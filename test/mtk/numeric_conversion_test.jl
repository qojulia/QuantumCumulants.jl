using QuantumCumulants
using QuantumOpticsBase
using Symbolics: @variables
using ModelingToolkitBase: @named, mtkcompile, ODEProblem, parameters, unknowns
using OrdinaryDiffEqTsit5: Tsit5, solve
using Random
using Test

Random.seed!(0)

# Values of an `initial_values` map in the order of `eqs.states`.
u0_vector(eqs, u0) = [u0[v] for v in values(moment_variable_map(eqs))]

@testset "initial_values: Fock + 2-level state to u0 map (integer levels)" begin
    h = FockSpace(:cavity) ⊗ NLevelSpace(:atom, 2)
    a = Destroy(h, :a)
    s(i, j) = Transition(h, :σ, i, j)

    @variables Δ::Real g::Real κ::Real η::Real
    H = Δ * a' * a + g * (a' * s(1, 2) + a * s(2, 1)) + η * (a + a')
    ops = [a, s(1, 2), a' * a, s(2, 2), a' * s(1, 2)]
    eqs = meanfield(ops, H, [a]; rates = [κ], order = 2)

    bcav = FockBasis(10)
    batom = NLevelBasis(2)
    b = bcav ⊗ batom
    ψ0 = randstate(b)

    u0 = u0_vector(eqs, initial_values(eqs, ψ0))

    @test u0[1] ≈ expect(destroy(bcav) ⊗ one(batom), ψ0)
    @test u0[2] ≈ expect(one(bcav) ⊗ QuantumOpticsBase.transition(batom, 1, 2), ψ0)
    @test u0[3] ≈ expect(number(bcav) ⊗ one(batom), ψ0)
    @test u0[4] ≈ expect(one(bcav) ⊗ QuantumOpticsBase.transition(batom, 2, 2), ψ0)
    @test u0[5] ≈ expect(create(bcav) ⊗ QuantumOpticsBase.transition(batom, 1, 2), ψ0)
end

@testset "System parameter order is deterministic" begin
    ha = NLevelSpace(:atom, 2)
    hf = FockSpace(:cavity)
    h = ha ⊗ hf

    a = Destroy(h, :a)
    σ(α, β, i) = IndexedOperator(Transition(h, :σ, α, β), i)

    @variables η::Real Γ::Real κ::Real Δc::Real N::Real
    g(i) = IndexedVariable(:g, i)
    i = Index(h, :i, N, ha)

    H = -Δc * a' * a +
        ∑(g(i) * (σ(2, 1, i) * a + σ(1, 2, i) * a'), i) +
        1im * η * (a' - a)
    J = [σ(1, 2, i), a]

    eqs = meanfield(a' * a, H, J; rates = [Γ, κ], order = [1, 2])
    complete!(eqs)
    evaled = evaluate(eqs; limits = (N => 3))

    sys = System(evaled; name = :deterministic_params)
    names = string.(parameters(sys))

    @test names == sort(names)
    @test Set(names) == Set(["Δc", "g", "η", "Γ", "κ"])
end

@testset "get_solution: untracked operator throws KeyError" begin
    hc = FockSpace(:cavity)
    a = Destroy(hc, :a)
    @variables ω::Real κ::Real
    H = ω * a' * a
    # Order-1 system tracks only ⟨a⟩; a second moment is not a state and has no conjugate
    # partner among the states, so the resolver exhausts every lookup and throws. Pass a
    # real solved `sol` (not `nothing`) so the test does not silently depend on
    # get_solution skipping `sol` on the lookup-miss path.
    eqs = meanfield([a], H, [a]; rates = [κ], order = 1)
    @named sys = System(eqs)
    sysc = mtkcompile(sys)
    u0 = Dict(unknowns(sysc) .=> ComplexF64[0.5])
    prob = ODEProblem(sysc, merge(u0, Dict(ω => 1.0, κ => 0.5)), (0.0, 1.0))
    sol = solve(prob, Tsit5())
    @test_throws KeyError get_solution(sol, a * a, eqs)
end

@testset "initial_values: LazyKet route matches Ket route" begin
    if isdefined(QuantumOpticsBase, :LazyKet)
        h = FockSpace(:cavity) ⊗ NLevelSpace(:atom, 2)
        a = Destroy(h, :a)
        s(i, j) = Transition(h, :σ, i, j)

        @variables Δ::Real g::Real κ::Real
        H = Δ * a' * a + g * (a' * s(1, 2) + a * s(2, 1))
        ops = [a, s(1, 2), a' * a, s(2, 2)]
        eqs = meanfield(ops, H, [a]; rates = [κ], order = 2)

        bcav = FockBasis(8)
        batom = NLevelBasis(2)
        b = bcav ⊗ batom
        ψlazy = LazyKet(b, (randstate(bcav), randstate(batom)))
        ψfull = Ket(ψlazy)

        u0_lazy = u0_vector(eqs, initial_values(eqs, ψlazy))
        u0_full = u0_vector(eqs, initial_values(eqs, ψfull))
        @test u0_lazy ≈ u0_full
    else
        @info "QuantumOpticsBase.LazyKet not available; skipping"
    end
end

@testset "initial_values from a state feeds ModelingToolkit and the kernel alike" begin
    h = FockSpace(:cavity) ⊗ NLevelSpace(:atom, 2)
    a = Destroy(h, :a)
    s(i, j) = Transition(h, :σ, i, j)
    @variables Δ::Real g::Real κ::Real
    H = Δ * a' * a + g * (a' * s(1, 2) + a * s(2, 1))
    eqs = complete(meanfield([a], H, [a]; rates = [κ], order = 2))
    ψ0 = coherentstate(FockBasis(8), 0.3) ⊗ nlevelstate(NLevelBasis(2), 2)
    ps = Dict(Δ => 0.1, g => 0.7, κ => 0.4)

    u0 = initial_values(eqs, ψ0)
    @test u0 isa AbstractDict
    @test u0_vector(eqs, u0) ≈ u0_vector(eqs, initial_values(eqs, u0_vector(eqs, u0)))

    sys = mtkcompile(System(eqs; name = :sys))
    sol_mtk = solve(ODEProblem(sys, merge(u0, ps), (0.0, 2.0)), Tsit5(); abstol = 1.0e-10, reltol = 1.0e-10)
    for u0_kernel in (ψ0, u0)
        prob = ODEProblem(eqs, u0_kernel, (0.0, 2.0), ps; backend = KernelBackend())
        sol = solve(prob, Tsit5(); abstol = 1.0e-10, reltol = 1.0e-10)
        @test get_solution(sol, a, eqs)(2.0) ≈ get_solution(sol_mtk, a, eqs)(2.0) rtol = 1.0e-6
    end
end
