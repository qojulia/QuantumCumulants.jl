using QuantumCumulants
using ModelingToolkitBase
using Symbolics: @variables
using OrdinaryDiffEqTsit5: ODEProblem, solve, Tsit5
using Test

@testset "heterodyne_detection" begin
    @variables N ωa γ η χ ωc κ g ξ ωl
    @register_symbolic _het_pulse(tt)
    hc = FockSpace(:resonator); ha = NLevelSpace(:atom, 2); h = hc ⊗ ha
    j = Index(h, :j, N, ha); k = Index(h, :k, N, ha)
    @qnumbers a::Destroy(h, 1)
    σ(α, β, idx) = IndexedOperator(Transition(h, :σ, α, β, 2), idx)
    eqs_seed = meanfield(a, ωc * a' * a, [a]; rates = [κ])
    tv = eqs_seed.iv
    H = ωc * a' * a + ωa * Σ(σ(2, 2, j), j) +
        g * a' * Σ(σ(1, 2, j), j) +
        g * a * Σ(σ(2, 1, j), j)
    J = [a * exp(1.0im * ωl * tv), σ(1, 2, j), σ(2, 1, j), σ(2, 2, j)]
    rates = [κ, γ, η * _het_pulse(tv), 2 * χ]
    efficiencies = [ξ, 0, 0, 0]
    ops = [a', a' * a, σ(2, 2, k), σ(1, 2, k), a * a]
    eqs = meanfield(
        ops, H, J;
        rates = rates, efficiencies = efficiencies,
        direction = Forward(), order = 2, iv = tv,
    )
    eqs_c = complete(eqs; get_adjoints = false)
    @test length(eqs_c.equations) == 13
    eqs_sc = scale(eqs_c)
    @test length(eqs_sc.equations) == 12
    sys_c = mtkcompile(System(eqs_sc; name = :het_sys))
    @test length(unknowns(sys_c)) == 12
end

@testset "heterodyne_detection (single-trajectory SDE bounded)" begin
    using StochasticDiffEq: SDEProblem, EM, RealWienerProcess
    using StochasticDiffEq.SciMLBase: ReturnCode
    import Random
    @variables N ωa γ η χ ωc κ g ξ ωl
    @register_symbolic _het_pulse2(tt)
    hc = FockSpace(:resonator); ha = NLevelSpace(:atom, 2); h = hc ⊗ ha
    j = Index(h, :j, N, ha); k = Index(h, :k, N, ha)
    @qnumbers a::Destroy(h, 1)
    σ(α, β, idx) = IndexedOperator(Transition(h, :σ, α, β, 2), idx)
    eqs_seed = meanfield(a, ωc * a' * a, [a]; rates = [κ])
    tv = eqs_seed.iv
    H = ωc * a' * a + ωa * Σ(σ(2, 2, j), j) +
        g * a' * Σ(σ(1, 2, j), j) +
        g * a * Σ(σ(2, 1, j), j)
    J = [a * exp(1.0im * ωl * tv), σ(1, 2, j), σ(2, 1, j), σ(2, 2, j)]
    rates = [κ, γ, η * _het_pulse2(tv), 2 * χ]
    efficiencies = [ξ, 0, 0, 0]
    ops = [a', a' * a, σ(2, 2, k), σ(1, 2, k), a * a]
    eqs = meanfield(
        ops, H, J;
        rates = rates, efficiencies = efficiencies,
        direction = Forward(), order = 2, iv = tv,
    )
    eqs_c = complete(eqs; get_adjoints = false)
    eqs_sc = scale(eqs_c)
    sys_st = mtkcompile(System(eqs_sc; name = :hetst_sys))
    ωc_ = 0.0; κ_ = 2π * 1.13e3; ξ_ = 0.12; N_ = 5.0e4
    ωa_ = 0.0; γ_ = 2π * 0.375; η_ = 2π * 20; χ_ = 0.016
    g_ = 2π * 0.73; ωl_ = 2π * 1.0e3
    t0 = 0.0; t1 = 20.0e-3
    _het_pulse2(tt) = (tt > t0 && tt < t0 + t1) * 1.0
    T_end = 0.1
    p_pairs = Dict(
        N => N_, ωa => ωa_, γ => γ_, η => η_, χ => χ_,
        ωc => ωc_, κ => κ_, g => g_, ξ => ξ_, ωl => ωl_,
    )
    u0_vec = zeros(ComplexF64, length(eqs_sc))
    dict_st = parameter_map(
        sys_st,
        merge(
            Dict{Any, Any}(unknowns(sys_st) .=> u0_vec),
            Dict{Any, Any}(p_pairs)
        ),
    )
    Random.seed!(2)
    noise = RealWienerProcess(0.0, 0.0)
    prob_st = SDEProblem(sys_st, dict_st, (0.0, T_end); noise = noise)
    sol = solve(prob_st, EM(); dt = T_end / 2.0e5)
    # The SDE must stay bounded; example trajectory peaks are ~600.
    adag_a_traj = real.(get_solution(sol, a' * a, eqs_sc).(sol.t))
    @test sol.retcode == ReturnCode.Success
    @test all(isfinite, adag_a_traj)
    @test maximum(abs, adag_a_traj) < 1.0e8
end

@testset "heterodyne_detection (drift-only numerical)" begin
    @variables N ωa γ η χ ωc κ g ωl
    @register_symbolic _het_pulse(tt)
    hc = FockSpace(:resonator); ha = NLevelSpace(:atom, 2); h = hc ⊗ ha
    j = Index(h, :j, N, ha); k = Index(h, :k, N, ha)
    @qnumbers a::Destroy(h, 1)
    σ(α, β, idx) = IndexedOperator(Transition(h, :σ, α, β, 2), idx)
    eqs_seed = meanfield(a, ωc * a' * a, [a]; rates = [κ])
    tv = eqs_seed.iv
    H = ωc * a' * a + ωa * Σ(σ(2, 2, j), j) +
        g * a' * Σ(σ(1, 2, j), j) +
        g * a * Σ(σ(2, 1, j), j)
    J = [a * exp(1.0im * ωl * tv), σ(1, 2, j), σ(2, 1, j), σ(2, 2, j)]
    rates = [κ, γ, η * _het_pulse(tv), 2 * χ]
    ops = [a', a' * a, σ(2, 2, k), σ(1, 2, k), a * a]
    eqs = meanfield(
        ops, H, J;
        rates = rates, direction = Forward(), order = 2, iv = tv,
    )
    eqs_c = complete(eqs; get_adjoints = false)
    eqs_sc = scale(eqs_c)
    sys_c = mtkcompile(System(eqs_sc; name = :het_det_sys))
    _het_pulse(_) = 1.0
    u0 = zeros(ComplexF64, length(eqs_sc.equations))
    init = initial_values(eqs_sc, u0)
    pmap = parameter_map(
        eqs_sc, Dict(
            ωc => 0.0, ωa => 0.0, ωl => 0.0,
            κ => 1.0, γ => 0.1, η => 1.0, χ => 0.1, g => 1.0, N => 4.0,
        )
    )
    prob = ODEProblem(sys_c, merge(init, pmap), (0.0, 5.0))
    sol = solve(prob, Tsit5(); reltol = 1.0e-8, abstol = 1.0e-8)
    @test isapprox(
        real(get_solution(sol, a' * a, eqs_sc).(sol.t[end])),
        1.906585813419867; rtol = 1.0e-4,
    )
    @test isapprox(
        real(get_solution(sol, eqs_sc.states[3], eqs_sc).(sol.t[end])),
        0.4702676694506166; rtol = 1.0e-4,
    )
end
