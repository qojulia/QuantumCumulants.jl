using QuantumCumulants
using QuantumOpticsBase
using ModelingToolkit
using OrdinaryDiffEq
using Test
using Random
using SteadyStateDiffEq
const qc = QuantumCumulants

@testset "indexed_meanfield" begin

    order = 2
    @cnumbers Δc η Δa κ

    N = 2 #number of atoms
    hc = FockSpace(:cavity)
    ha = NLevelSpace(Symbol(:atom), 2)
    h = hc ⊗ ha

    #define indices
    i_ind = Index(h, :i, N, ha)
    j_ind = Index(h, :j, N, ha)
    k_ind = Index(h, :k, N, ha)

    #define indexed variables
    g(k) = IndexedVariable(:g, k)
    Γ_ij = DoubleIndexedVariable(:Γ, i_ind, j_ind)
    Ω_ij = DoubleIndexedVariable(:Ω, i_ind, j_ind; identical = false)

    @qnumbers a::Destroy(h)
    σ(i, j, k) = IndexedOperator(Transition(h, :σ, i, j), k)

    # Hamiltonian

    DSum = Σ(Ω_ij*σ(2, 1, i_ind)*σ(1, 2, j_ind), j_ind, i_ind; non_equal = true)

    @test DSum isa DoubleSum
    @test isequal(Σ(Σ(Ω_ij*σ(2, 1, i_ind)*σ(1, 2, j_ind), i_ind, [j_ind]), j_ind), DSum)

    Hc = Δc*a'a + η*(a' + a)
    Ha = Δa*Σ(σ(2, 2, i_ind), i_ind) + DSum
    Hi = Σ(g(i_ind)*(a'*σ(1, 2, i_ind) + a*σ(2, 1, i_ind)), i_ind)
    H = Hc + Ha + Hi

    J = [a, [σ(1, 2, i_ind), σ(1, 2, j_ind)]]
    rates = [κ, Γ_ij]

    J_2 = [a, σ(1, 2, i_ind)]
    rates_2 = [κ, Γ_ij]

    ops = [a, σ(2, 2, k_ind), σ(1, 2, k_ind)]
    eqs = meanfield(ops, H, J; rates = rates, order = order)

    eqs_2 = meanfield(ops, H, J_2; rates = rates_2, order = order)

    @test eqs.equations == eqs_2.equations

    @test isequal([i_ind, j_ind, k_ind], sort(qc.get_indices_equations(eqs)))
    @test isequal([:i, :j, :k], sort(qc.getIndName.(qc.get_indices_equations(eqs))))

    @test length(eqs) == 3

    ind1 = Index(h, :q, N, ha)
    ind2 = Index(h, :r, N, ha)
    ind3 = Index(h, :s, N, ha)

    eqs_comp = qc.complete(eqs; extra_indices = [ind1, ind2, ind3])
    eqs_comp2 = qc.complete(eqs)

    @test length(eqs_comp.equations) == length(eqs_comp2.equations)

    eqs_ = evaluate(eqs_comp)
    eqs_2 = evaluate(eqs_comp2)

    @test length(eqs_2) == length(eqs_)

    @test length(eqs_) == 18

    eqs_ord = meanfield(ops, H, J; rates = rates)
    @test eqs_ord.order === nothing
    eqs_ord_c = complete(eqs_ord)
    eqs_ord_c2 = complete(eqs_ord; order = 2)
    @test isequal(eqs_ord_c.states, eqs_ord_c2.states)

    @named sys = System(eqs_);

    u0 = zeros(ComplexF64, length(eqs_))
    # parameter
    Γ_ = 1.0
    d = 2π*0.8 #0.8λ
    θ = π/2

    Ωij_(i, j) =
        Γ_*(-3/4)*((1-(cos(θ))^2)*cos(d)/d-(1-3*(cos(θ))^2)*(sin(d)/(d^2)+(cos(d)/(d^3))))
    function Γij_(i, j)
        i==j ? Γ_ :
        Γ_*(3/2)*((1-(cos(θ))^2)*sin(d)/d+(1-3*(cos(θ))^2)*((cos(d)/(d^2))-(sin(d)/(d^3))))
    end

    ΓMatrix = [Γij_(i, j) for i = 1:2, j = 1:2]
    ΩMatrix = [Ωij_(i, j) for i = 1:2, j = 1:2]


    g_ = 2Γ_
    κ_ = 3Γ_
    Δa_ = 0Γ_
    Δc_ = 0Γ_
    η_ = κ_/10

    g_v = [g_*(-1)^j for j = 1:2]
    ps = [Δc, η, Δa, κ, g(i_ind), Γ_ij, Ω_ij];

    eqs_4 = meanfield(ops, H, J; rates = rates, order = 4)

    Δc_i = -3*Γ_
    Δa_i = Δc_i + Ωij_(1, 2) #cavity on resonace with the shifted collective emitter
    p0_i = [Δc_i, η_, Δa_i, κ_, g_v, ΓMatrix, ΩMatrix]

    ps_ = value_map(ps, p0_i) #Combine all the parameters + values to one list for solving
    prob = ODEProblem(sys, u0, (0.0, 10Γ_), ps_);

    sol_ss =
        solve(prob, Tsit5(), save_everystep = false, save_on = false, save_start = false)

    @test length(eqs_4) == length(eqs)

    @test get_solution(sol_ss, 2a + σ(1, 1, 1))[1] ≈
          (2*sol_ss[a][1] + 1 - sol_ss[σ(2, 2, 1)][1]) ≈
          get_solution(sol_ss, average(2a + σ(1, 1, 1)))[1]
    @test isequal(σ(1, 1, 2), 1-σ(2, 2, 2))

    order = 1
    @cnumbers g_ N κ

    # Hilbertspace
    hc = FockSpace(:cavity)
    hf = FockSpace(:filter)

    h = hc ⊗ hf

    i = Index(h, :i, N, hf)
    j = Index(h, :j, N, hf)
    k = Index(h, :k, N, hf)

    xij = IndexedVariable(:x, i, j)

    @qnumbers a_::Destroy(h, 1)
    b(k) = IndexedOperator(Destroy(h, :b, 2), k)

    H = g_*a_'a_ + Σ(xij*a_'a_*b(i)'b(j), i, j)
    J = [a_]
    rates = [κ]


    eqs1 = meanfield(a_, H, J; rates = rates, order = order)
    eqs2 = meanfield([a_], H, J; rates = rates, order = order)
    @test isequal(eqs1.equations, eqs2.equations)

    @test isequal(sort([i, j]), sort(qc.get_all_indices(eqs1)))


    #example for testing evaluation of individual hilbertspaces
    @cnumbers N N2 Δ κ Γ R ν M

    hc = FockSpace(:cavity)
    ha = NLevelSpace(:atom, 2)

    h = hc ⊗ ha

    k = Index(h, :k, N, ha)
    l = Index(h, :l, N, ha)

    m = Index(h, :m, N2, hc)
    n = Index(h, :n, N2, hc)

    order = 2

    σ(i, j, k) = IndexedOperator(Transition(h, :σ, i, j), k)
    ai(k) = IndexedOperator(Destroy(h, :a), k)

    H_2 =
        -Δ*∑(ai(m)'ai(m), m) +
        g_*(∑(Σ(ai(m)'*σ(1, 2, k), k), m) + ∑(Σ(ai(m)*σ(2, 1, k), k), m))

    J_2 = [ai(m), σ(1, 2, k), σ(2, 1, k), σ(2, 2, k)]
    rates_2 = [κ, Γ, R, ν]
    ops_2 = [ai(n)'*ai(n), σ(2, 2, l)]
    eqs_2 = meanfield(ops_2, H_2, J_2; rates = rates_2, order = order)

    q = Index(h, :q, N, ha)
    r = Index(h, :r, N2, hc)

    extra_indices = [q, r]

    eqs_com = qc.complete(eqs_2; extra_indices = extra_indices);
    eqs_com2 = qc.complete(eqs_2)

    @test length(eqs_com) == 15
    @test length(eqs_com) == length(eqs_com2)

    @test isequal(sort([k, m, n, l, q, r]), sort(qc.get_all_indices(eqs_com)))

    e_1 = evaluate(eqs_com; h = 2, limits = (N=>5))
    e_2 = evaluate(eqs_com; h = 1, limits = (N2=>6))

    @test length(e_1) != length(e_2)
    @test !(e_1.equations == e_2.equations)
    @test !(e_1.states == e_2.states)

    @test sort(qc.get_indices_equations(e_1)) == sort([m, n, r])
    @test sort(qc.get_indices_equations(e_2)) == sort([k, l, q])

    limits = Dict(N=>5, N2=>6)
    s1 = evaluate(eqs_com; h = [1, 2], limits = limits)
    s2 = evaluate(eqs_com; limits = limits)

    s1_ = evaluate(eqs_com; h = [ha], limits = (N=>5))
    s2_ = evaluate(eqs_com; h = [2], limits = (N=>5))

    @test s1_.equations == s2_.equations

    @test length(s1) == length(s2)
    @test s1.equations == s2.equations

    @test qc.get_indices_equations(s1) == []

    @test s1.states == s2.states

    # test for linear combinations in jump operators
    ha = NLevelSpace(:atoms, 2)
    hc = FockSpace(:cavity)
    h = hc ⊗ ha

    @cnumbers N Δ κ γ ν χ gg
    i = Index(h, :i, N, ha)
    j = Index(h, :j, N, ha)
    N_ = 3

    @qnumbers b_::Destroy(h)
    σ(x, y, z) = IndexedOperator(Transition(h, :σ, x, y), z)
    gi = IndexedVariable(:g, i)

    H = Δ*b_'*b_ + ∑(gi*(b_*σ(2, 1, i) + b_'*σ(1, 2, i)), i)
    ops = [b_'b_, σ(2, 2, j)]

    J = [b_, σ(1, 2, i), σ(2, 1, i), σ(2, 2, i)]
    rates = [κ, γ, ν, χ]

    eqs = meanfield(ops, H, J; rates = rates, order = 2)
    eqs_c = complete(eqs);

    #test with scale
    abstol=1e-8
    reltol=1e-8
    eqs_sc = scale(eqs_c);
    @named sys_sc = System(eqs_sc);
    u0 = zeros(ComplexF64, length(eqs_sc))
    u0[1] = 2 # initial value for ⟨b'*b⟩
    ps_sc = [IndexedVariable(:g, 1), N, Δ, κ, γ, ν, χ]
    p_sc = [0.5, N_, 0.0, 1.0, 1.0, 0.5, 2.0]
    P_sc = value_map(ps_sc, p_sc);
    prob_sc = ODEProblem(sys_sc, u0, (0.0, 10.0), P_sc);
    sol_sc = solve(prob_sc, Tsit5(); abstol, reltol);
    @test sol_sc isa ODESolution

    #test with evaluate
    eqs_eval = evaluate(eqs_c; limits = (N=>N_));
    @named sys = System(eqs_eval);
    u0 = zeros(ComplexF64, length(eqs_eval))
    u0[1] = 2
    ps = [gi, Δ, κ, γ, ν, χ]
    p = [[0.5, 0.5, 0.5], 0.0, 1.0, 1.0, 0.5, 2.0]
    P = value_map(ps, p; limits = (N=>N_));
    prob = ODEProblem(sys, u0, (0.0, 10.0), P);
    sol = solve(prob, Tsit5(); abstol, reltol);
    @test sol isa ODESolution

    #test comparison straight forward (without indexing)
    ha_(i) = NLevelSpace(Symbol(:atoms, i), 2)
    ha = tensor([ha_(i) for i = 1:N_]...)
    hc = FockSpace(:cavity)
    hh = hc ⊗ ha

    @qnumbers b2::Destroy(hh)
    s(x, y, i) = Transition(hh, Symbol(:s, i), x, y, i+1)
    s_ls(x, y) = [s(x, y, i) for i = 1:N_]

    H_ = Δ*b2'*b2 + sum(gg*(b2*s(2, 1, i) + b2'*s(1, 2, i)) for i = 1:N_)
    ops_ = [b2'b2]
    J_ = [b2; s_ls(1, 2); s_ls(2, 1); s_ls(2, 2)]
    rates_ = [κ; [γ for i = 1:N_]; [ν for i = 1:N_]; [χ for i = 1:N_]]
    eqs_ = meanfield(ops_, H_, J_; rates = rates_, order = 2)
    eqs_c_ = complete(eqs_);

    @named sys_ = System(eqs_c_)
    u0_ = zeros(ComplexF64, length(eqs_c_))
    u0_[1] = 2
    ps_ = [gg, Δ, κ, γ, ν, χ]
    p_ = [0.5, 0.0, 1.0, 1.0, 0.5, 2.0]
    prob_ = ODEProblem(sys_, u0_, (0, 10.0), ps_ .=> p_)
    sol_ = solve(prob_, Tsit5(); abstol, reltol)

    @test sol_sc[b_'b_][end] ≈ sol[b_'b_][end] ≈ sol_[b2'b2][end]
    @test sol_sc[σ(1, 2, 1)][end] ≈ sol[σ(1, 2, 1)][end] ≈ sol_[s(1, 2, 1)][end]
    @test sol_sc[b_'σ(1, 2, 1)][end] ≈ sol[b_'σ(1, 2, 1)][end] ≈ sol_[b2's(1, 2, 1)][end]
    @test sol_sc[σ(2, 2, 1)][end] ≈ sol[σ(2, 2, 1)][end] ≈ sol_[s(2, 2, 1)][end]

    ### another test case
    # Hilbertspace
    hc = FockSpace(:cavity)
    ha = NLevelSpace(:atom, 2)
    h = hc ⊗ ha
    # operators
    @qnumbers a::Destroy(h)
    σ(α, β, i) = IndexedOperator(Transition(h, :σ, α, β), i)
    @cnumbers N Δ κ Γ R ν
    g(i) = IndexedVariable(:g, i)
    i = Index(h, :i, N, ha)
    j = Index(h, :j, N, ha)
    # Hamiltonian
    H = -Δ*a'a + Σ(g(i)*(a'*σ(1, 2, i) + a*σ(2, 1, i)), i)
    # Jump operators with corresponding rates
    J = [a, σ(1, 2, i), σ(2, 1, i), σ(2, 2, i)]
    rates = [κ, Γ, R, ν]
    # Derive equations
    ops = [a, σ(2, 2, j)]
    eqs_i = meanfield(ops, H, J; rates = rates)

    # Hilbertspace
    hc = FockSpace(:cavity)
    ha = ⊗([NLevelSpace(:atom, 2) for i = 1:3]...)
    h = hc ⊗ ha
    # operators
    @qnumbers a::Destroy(h)
    σ(α, β, k) = Transition(h, Symbol("σ$(k)"), α, β, k+1)
    @cnumbers N Δ κ Γ R ν gg
    # Hamiltonian
    H = -Δ*a'a + sum(gg*(a'*σ(1, 2, i) + a*σ(2, 1, i)) for i = 1:3)
    # Jump operators with corresponding rates
    J = [a; [σ(1, 2, i) for i = 1:3]; [σ(2, 1, i) for i = 1:3]; [σ(2, 2, i) for i = 1:3]]
    rates = [κ; [Γ for i = 1:3]; [R for i = 1:3]; [ν for i = 1:3]]
    # Derive equations
    ops = [a, σ(2, 2, 1)]
    eqs_os = meanfield(ops, H, J; rates = rates)

    eqs_i_c = complete(eqs_i)
    eqs_i_c2 = complete(eqs_i; order = 2)
    eqs_os_c = complete(eqs_os)
    eqs_os_c2 = complete(eqs_os; order = 2)

    @test eqs_i_c.order === nothing
    @test eqs_i_c2.order === nothing
    @test eqs_os_c.order === nothing
    @test eqs_os_c2.order === nothing
    @test isequal(eqs_i_c.states, eqs_i_c2.states)
    @test isequal(eqs_os_c.states, eqs_os_c2.states)
    eqs_i_ev = evaluate(eqs_i_c; limits = (N=>3))
    eqs_i_ev2 = evaluate(eqs_i_c2; limits = (N=>3))
    @test isequal(length(eqs_i_ev.states), length(eqs_os_c.states))
    @test isequal(length(eqs_i_ev2.states), length(eqs_os_c2.states))

    @test_throws MethodError complete(eqs_i; order = 1)
    @test_throws MethodError complete(eqs_os; order = 1)

end
