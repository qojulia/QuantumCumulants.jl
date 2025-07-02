using QuantumCumulants
using QuantumOpticsBase
using ModelingToolkit: System
using OrdinaryDiffEq
using Test
using Random;
Random.seed!(0)

@testset "numeric-conversion" begin

    # Test fock basis conversion
    hfock = FockSpace(:fock)
    @qnumbers a::Destroy(hfock)
    bfock = FockBasis(7)
    @test to_numeric(a, bfock) == destroy(bfock)
    @test to_numeric(a', bfock) == create(bfock)

    # NLevelSpace conversion
    hnlevel = NLevelSpace(:nlevel, 3)
    σ(i, j) = Transition(hnlevel, :σ, i, j)
    bnlevel = NLevelBasis(3)
    for i = 1:3, j = 1:3
        op = σ(i, j)
        @test to_numeric(op, bnlevel) == QuantumOpticsBase.transition(bnlevel, i, j)
    end

    # with symbolic levels
    levels = (:g, :e, :a)
    hnlevel_sym = NLevelSpace(:nlevel_sym, levels)
    σ_sym(i, j) = Transition(hnlevel_sym, :σ, i, j)
    @test_throws ArgumentError to_numeric(σ_sym(:e, :g), bnlevel)
    level_map = Dict((levels .=> (1, 2, 3))...)
    for i = 1:3, j = 1:3
        lvl1 = levels[i]
        lvl2 = levels[j]
        op = σ_sym(lvl1, lvl2)
        @test to_numeric(op, bnlevel; level_map = level_map) ==
              QuantumOpticsBase.transition(bnlevel, i, j)
    end

    # On composite bases
    hprod = hfock ⊗ hnlevel
    a = Destroy(hprod, :a)
    σprod(i, j) = Transition(hprod, :σ, i, j)
    bprod = bfock ⊗ bnlevel
    for i = 1:3, j = 1:3
        i == j == 1 && continue  # rewritten as sum, see below
        op1 = a*σprod(i, j)
        op2 = a'*σprod(i, j)
        @test to_numeric(op1, bprod) == LazyTensor(
            bprod,
            [1, 2],
            (destroy(bfock), QuantumOpticsBase.transition(bnlevel, i, j)),
        )
        @test to_numeric(op2, bprod) == LazyTensor(
            bprod,
            [1, 2],
            (create(bfock), QuantumOpticsBase.transition(bnlevel, i, j)),
        )
    end

    op1_num = to_numeric(a*σprod(1, 1), bprod)
    @test op1_num isa LazySum
    @test sparse(op1_num) == destroy(bfock) ⊗ QuantumOpticsBase.transition(bnlevel, 1, 1)

    op2_num = to_numeric(a'*σprod(1, 1), bprod)
    @test op2_num isa LazySum
    @test sparse(op2_num) == create(bfock) ⊗ QuantumOpticsBase.transition(bnlevel, 1, 1)

    @test to_numeric(a'*a, bprod) ==
          LazyTensor(bprod, [1], (create(bfock) * destroy(bfock),))

    # Composite basis with symbolic levels
    σsym_prod(i, j) = Transition(hfock ⊗ hnlevel_sym, :σ, i, j)
    a = Destroy(hfock ⊗ hnlevel_sym, :a)
    @test_throws ArgumentError to_numeric(a*σsym_prod(:e, :g), bprod)
    for i = 1:3, j = 1:3
        i == j == 1 && continue  # see below
        op1 = a*σsym_prod(levels[i], levels[j])
        op2 = a'*σsym_prod(levels[i], levels[j])
        @test to_numeric(op1, bprod; level_map = level_map) == LazyTensor(
            bprod,
            [1, 2],
            (destroy(bfock), QuantumOpticsBase.transition(bnlevel, i, j)),
        )
        @test to_numeric(op2, bprod; level_map = level_map) == LazyTensor(
            bprod,
            [1, 2],
            (create(bfock), QuantumOpticsBase.transition(bnlevel, i, j)),
        )
    end

    op1_num = to_numeric(a*σsym_prod(:g, :g), bprod; level_map = level_map)
    @test op1_num isa LazySum
    @test sparse(op1_num) == destroy(bfock) ⊗ QuantumOpticsBase.transition(bnlevel, 1, 1)

    op2_num = to_numeric(a'*σsym_prod(:g, :g), bprod; level_map = level_map)
    @test op2_num isa LazySum
    @test sparse(op2_num) == create(bfock) ⊗ QuantumOpticsBase.transition(bnlevel, 1, 1)

    # composite basis with a "gap"
    hprod_gap = hfock ⊗ hnlevel ⊗ hnlevel
    bprod_gap = bprod ⊗ bnlevel
    a = Destroy(hprod_gap, :a)
    σprod_gap(i, j) = Transition(hprod_gap, :σ, i, j, 3)
    for i = 1:3, j = 1:3
        i == j == 1 && continue
        op1 = a*σprod_gap(i, j)
        op2 = a'*σprod_gap(i, j)
        @test to_numeric(op1, bprod_gap) == LazyTensor(
            bprod_gap,
            [1, 3],
            (destroy(bfock), QuantumOpticsBase.transition(bnlevel, i, j)),
        )
        @test to_numeric(op2, bprod_gap) == LazyTensor(
            bprod_gap,
            [1, 3],
            (create(bfock), QuantumOpticsBase.transition(bnlevel, i, j)),
        )
    end

    op1_num = to_numeric(a*σprod_gap(1, 1), bprod_gap)
    @test op1_num isa LazySum
    @test sparse(op1_num) ==
          destroy(bfock) ⊗ one(bnlevel) ⊗ QuantumOpticsBase.transition(bnlevel, 1, 1)

    op2_num = to_numeric(a'*σprod_gap(1, 1), bprod_gap)
    @test op2_num isa LazySum
    @test sparse(op2_num) ==
          create(bfock) ⊗ one(bnlevel) ⊗ QuantumOpticsBase.transition(bnlevel, 1, 1)


    # Numeric average values
    α = 0.1 + 0.2im
    ψ = coherentstate(bfock, α)
    a = Destroy(hfock, :a)
    @test numeric_average(a, ψ) ≈ numeric_average(average(a), ψ) ≈ α
    @test numeric_average(a'*a, ψ) ≈ numeric_average(average(a'*a), ψ) ≈ abs2(α)

    ψprod = ψ ⊗ nlevelstate(bnlevel, 1)
    @test_throws ArgumentError numeric_average(σsym_prod(:e, :g), ψprod)
    idfock = one(bfock)
    for i = 1:3, j = 1:3
        op = σprod(i, j)
        op_sym = σsym_prod(levels[i], levels[j])
        op_num = idfock ⊗ QuantumOpticsBase.transition(bnlevel, i, j)
        @test numeric_average(op, ψprod) ≈ expect(op_num, ψprod)
        @test numeric_average(op_sym, ψprod; level_map = level_map) ≈ expect(op_num, ψprod)
    end

    # LazyKet
    if isdefined(QuantumOpticsBase, :LazyKet)
        ψlazy = LazyKet(bprod, (ψ, nlevelstate(bnlevel, 1)))
        @test_throws ArgumentError numeric_average(σsym_prod(:e, :g), ψlazy)
        for i = 1:3, j = 1:3
            op = σprod(i, j)
            op_sym = σsym_prod(levels[i], levels[j])
            op_num = LazyTensor(bprod, [2], (QuantumOpticsBase.transition(bnlevel, i, j),))
            @test numeric_average(op, ψlazy) ≈ expect(op_num, ψlazy)
            @test numeric_average(op_sym, ψlazy; level_map = level_map) ≈
                  expect(op_num, ψlazy)
        end
    end

    # Initial values in actual equations
    levels = (:g, :e)
    h = FockSpace(:cavity) ⊗ NLevelSpace(:atom, levels)
    a = Destroy(h, :a)
    s(i, j) = Transition(h, :σ, i, j)

    @cnumbers Δ g κ γ η

    H = Δ*a'*a + g*(a'*s(:g, :e) + a*s(:e, :g)) + η*(a + a')
    ops = [a, s(:g, :e), a'*a, s(:e, :e), a'*s(:g, :e)]
    eqs = meanfield(ops, H, [a]; rates = [κ], order = 2)

    bcav = FockBasis(10)
    batom = NLevelBasis(2)
    b = bcav ⊗ batom
    ψ0 = randstate(b)

    level_map = Dict((levels .=> [1, 2])...)
    u0 = initial_values(eqs, ψ0; level_map = level_map)

    @test u0[1] ≈ expect(destroy(bcav) ⊗ one(batom), ψ0)
    @test u0[2] ≈ expect(one(bcav) ⊗ QuantumOpticsBase.transition(batom, 1, 2), ψ0)
    @test u0[3] ≈ expect(number(bcav) ⊗ one(batom), ψ0)
    @test u0[4] ≈ expect(one(bcav) ⊗ QuantumOpticsBase.transition(batom, 2, 2), ψ0)
    @test u0[5] ≈ expect(create(bcav) ⊗ QuantumOpticsBase.transition(batom, 1, 2), ψ0)


    if isdefined(QuantumOpticsBase, :LazyKet)
        ψlazy = LazyKet(b, (randstate(bcav), randstate(batom)))
        ψfull = Ket(ψlazy)
        u0 = initial_values(eqs, ψlazy; level_map = level_map)

        @test u0[1] ≈ expect(destroy(bcav) ⊗ one(batom), ψfull)
        @test u0[2] ≈ expect(one(bcav) ⊗ QuantumOpticsBase.transition(batom, 1, 2), ψfull)
        @test u0[3] ≈ expect(number(bcav) ⊗ one(batom), ψfull)
        @test u0[4] ≈ expect(one(bcav) ⊗ QuantumOpticsBase.transition(batom, 2, 2), ψfull)
        @test u0[5] ≈
              expect(create(bcav) ⊗ QuantumOpticsBase.transition(batom, 1, 2), ψfull)
    end

    # Test sufficiently large hilbert space; from issue #109
    hfock = FockSpace(:fock)
    @qnumbers a::Destroy(hfock)
    bfock = FockBasis(100)

    diff = (2*create(bfock)+2*destroy(bfock)) - to_numeric((2*(a)+2*(a')), bfock)
    @test isequal(2*create(bfock)+2*destroy(bfock), to_numeric((2*(a)+2*(a')), bfock))
    @test iszero(diff)

    @test isequal(to_numeric(2*a, bfock), 2*to_numeric(a, bfock))
    @test iszero(to_numeric(2*a, bfock) - 2*to_numeric(a, bfock))

    # Test indexed initial state (superradiant pulse)
    order = 2 #order of the cumulant expansion
    @cnumbers κ g Γ Δ N
    hc = FockSpace(:cavity)
    ha = NLevelSpace(:atom, 2)
    h = hc ⊗ ha
    i = Index(h, :i, N, ha)
    j = Index(h, :j, N, ha)
    k = Index(h, :k, N, ha)
    @qnumbers a::Destroy(h, 1)
    σ(α, β, k) = IndexedOperator(Transition(h, :σ, α, β, 2), k)
    bc = FockBasis(3)
    ba = NLevelBasis(2)
    b = tensor(bc, [ba for i = 1:order]...)
    ψc = fockstate(bc, 0)
    ψa = normalize(nlevelstate(ba, 1) + nlevelstate(ba, 2))
    ψ = tensor(ψc, [ψa for i = 1:order]...)
    a_ = LazyTensor(b, [1], (destroy(bc),))
    σ_(i, j, k) = LazyTensor(b, [1+k], (QuantumOpticsBase.transition(ba, i, j),))
    ranges=[1, 2]
    @test to_numeric(σ(1, 2, 1), b; ranges = ranges) == σ_(1, 2, 1)
    @test to_numeric(σ(2, 2, 2), b; ranges = ranges) == σ_(2, 2, 2)
    @test to_numeric(a, b; ranges = ranges) == a_
    @test to_numeric(a*σ(2, 2, 2), b; ranges = ranges) == σ_(2, 2, 2)*a_
    @test numeric_average(σ(2, 2, 2), ψ; ranges = ranges) ≈ 0.5
    @test numeric_average(average(σ(2, 2, 1)), ψ; ranges = ranges) ≈ 0.5
    @test numeric_average(average(a'a), ψ; ranges = ranges) ≈ 0.0
    @test numeric_average(average(a*σ(2, 2, 1)), ψ; ranges = ranges) ≈ 0.0
    @test_throws ArgumentError numeric_average(average(a'a), ψ)

    if isdefined(QuantumOpticsBase, :LazyKet)
        ψlazy = LazyKet(b, (ψc, (ψa for i = 1:order)...))
        @test numeric_average(σ(2, 2, 2), ψlazy; ranges = ranges) ≈ 0.5
        @test numeric_average(average(σ(2, 2, 1)), ψlazy; ranges = ranges) ≈ 0.5
        @test numeric_average(average(a'a), ψlazy; ranges = ranges) ≈ 0.0
        @test numeric_average(average(a*σ(2, 2, 1)), ψlazy; ranges = ranges) ≈ 0.0
        @test_throws ArgumentError numeric_average(average(a'a), ψlazy)
    end

    # Hamiltonian
    H = -Δ*a'a + g*(Σ(a'*σ(1, 2, i), i) + Σ(a*σ(2, 1, i), i))
    J = [a, σ(1, 2, i)]
    rates = [κ, Γ]
    ops = [a, σ(2, 2, j)]
    eqs = meanfield(ops, H, J; rates = rates, order = order)
    eqs_c = complete(eqs)
    eqs_sc = scale(eqs_c)
    @named sys = System(eqs_sc)
    @test_throws ArgumentError initial_values(eqs_sc, ψ)

    u0 = initial_values(eqs_sc, ψ; ranges = ranges)

    if isdefined(QuantumOpticsBase, :LazyKet)
        ψlazy = LazyKet(b, (ψc, (ψa for i = 1:order)...))
        @test u0 ≈ initial_values(eqs_sc, ψlazy; ranges = ranges)
    end

    N_ = 2e5
    Γ_ = 1.0 #Γ=1mHz
    Δ_ = 2500Γ_ #Δ=2.5Hz
    g_ = 1000Γ_ #g=1Hz
    κ_ = 5e6*Γ_ #κ=5kHz
    ps = [N, Δ, g, κ, Γ]
    p0 = [N_, Δ_, g_, κ_, Γ_]

    prob = ODEProblem(sys, u0, (0.0, 1e-4/Γ_), ps .=> p0)
    sol = solve(prob, Tsit5(); save_on = false, save_everystep = false)
    @test sol.retcode == ReturnCode.Success

    # Pauli and Spin in tes_spin.jl

end # testset
