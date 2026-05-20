using QuantumCumulants
using Symbolics: Symbolics, @variables
using SymbolicUtils
using ModelingToolkitBase: ModelingToolkitBase, @named, mtkcompile, ODEProblem, unknowns
using OrdinaryDiffEq: Tsit5, solve, ReturnCode
using Test

@testset "CorrelationFunction: damped cavity construction" begin
    hc = FockSpace(:cavity)
    @qnumbers a::Destroy(hc)
    @variables ω κ
    H = ω * a' * a
    eqs = meanfield([a], H, [a]; rates = [κ])
    complete!(eqs)
    c = CorrelationFunction(a', a, eqs)
    @test c isa CorrelationFunction
    @test c.op1 === a'
    @test c.op2 === a
    @test length(c.eqs.equations) >= 1
end

@testset "Spectrum constructor" begin
    hc = FockSpace(:cavity)
    @qnumbers a::Destroy(hc)
    @variables ω κ
    H = ω * a' * a
    eqs = meanfield([a], H, [a]; rates = [κ])
    complete!(eqs)
    c = CorrelationFunction(a', a, eqs)
    s = Spectrum(c)
    @test s isa Spectrum
end

@testset "Damped cavity: Spectrum matches analytic Lorentzian" begin
    # Closed-form check on the full pipeline.
    # For H = ωc·a†a, J = [a] with rate κ, vacuum bath:
    #   ⟨a(τ)·a†(0)⟩_ss = exp(-(iωc + κ/2)τ),  initial value ⟨a·a†⟩_ss = 1.
    #   S(ω) = Re ∫₀^∞ ⟨a(τ)a†(0)⟩ exp(-iωτ) dτ = (κ/2)/((ω + ωc)² + (κ/2)²).
    hc = FockSpace(:cavity)
    a = Destroy(hc, :a)
    @variables ωc::Real κ::Real
    H = ωc * a' * a
    he = meanfield(a' * a, H, [a]; rates = [κ])
    complete!(he)

    # CorrelationFunction(a, a', he) ⇒ ⟨a(τ)·a†(0)⟩ correlator.
    c = CorrelationFunction(a, a', he)
    @test length(c.eqs.equations) == 1   # single state: ⟨a·a†_anc⟩

    # τ=0 initial value (vacuum steady state ⟨a†a⟩_ss = 0 ⇒ ⟨a·a†⟩_ss = 1).
    u_end = Dict(SymbolicUtils.unwrap(average(a' * a)) => 0.0 + 0im)
    u0 = correlation_u0(c, u_end)
    @test only(values(u0)) ≈ 1.0

    # Spectrum vs analytic Lorentzian.
    ωc_num = 1.0
    κ_num = 0.5
    ps = (ωc, κ)
    pn = ComplexF64[ωc_num, κ_num]
    ω = collect(range(-2π, 2π; length = 81))
    S = Spectrum(c, ps)
    s_num = S(ω, u_end, pn)
    s_an = @. (κ_num / 2) / ((ω + ωc_num)^2 + (κ_num / 2)^2)
    @test maximum(abs.(s_num .- s_an)) < 1e-10
end

@testset "Driven cavity: ODE conjugate substitution is non-folding" begin
    # Regression test for a subtle correctness bug: meanfield equations
    # frequently contain `⟨X⟩ - ⟨X†⟩` terms. The `_conj_substitution_dict`
    # pass in `to_system` maps `⟨X†⟩ → _qc_conj(avg_X(t))`. Using bare
    # `conj(...)` folds to `avg_X(t)` (because state vars carry `Real`
    # symtype), silently zeroing the entire driving term. Here:
    #   d⟨a'a⟩/dt = -κ⟨a'a⟩ + iΩ(⟨a⟩ - ⟨a'⟩)
    # The closed-form steady state is ⟨a'a⟩_ss = |⟨a⟩_ss|² where
    # ⟨a⟩_ss = -iΩ/(iωc + κ/2). Verifies the ODE actually integrates the
    # conjugate-paired drive.
    hc = FockSpace(:cavity)
    a = Destroy(hc, :a)
    @variables ωc::Real κ::Real Ω::Real
    H = ωc * a' * a + Ω * (a + a')
    he = meanfield([a, a'a], H, [a]; rates = [κ])
    complete!(he)

    @named sys = to_system(he)
    sys_c = mtkcompile(sys)
    ωc_n, κ_n, Ω_n = 1.0, 0.5, 0.3
    ps = [ωc, κ, Ω]
    pn = [ωc_n, κ_n, Ω_n]
    u0 = Dict(unknowns(sys_c) .=> zeros(ComplexF64, length(unknowns(sys_c))))
    sol = solve(
        ODEProblem(sys_c, merge(u0, Dict(ps .=> pn)), (0.0, 80.0)),
        Tsit5(); abstol = 1e-10, reltol = 1e-10,
    )
    @test sol.retcode == ReturnCode.Success

    a_ss_an = -im * Ω_n / (im * ωc_n + κ_n / 2)
    n_ss_an = abs2(a_ss_an)
    @test abs(get_solution(sol, a, he)(80.0) - a_ss_an) < 1e-7
    @test abs(get_solution(sol, a' * a, he)(80.0) - n_ss_an) < 1e-7
end

@testset "_lookup_avg resolves compound (normal-ordered) expressions" begin
    # Regression test for the bug found while writing the Lorentzian check:
    # `_undo_ancilla` rebuilds `a · a'_anc → a · a'` via SQA arithmetic, which
    # normal-orders to `a'·a + 1`. The internal `_lookup_avg` must resolve
    # such compound expressions against `u_end`, not return 0.
    hc = FockSpace(:cavity)
    a = Destroy(hc, :a)
    @variables ωc::Real κ::Real
    H = ωc * a' * a
    he = meanfield(a' * a, H, [a]; rates = [κ])
    complete!(he)

    c = CorrelationFunction(a, a', he)
    # Non-zero ⟨a'·a⟩ to ensure the +1 offset is detectable.
    u_end = Dict(SymbolicUtils.unwrap(average(a' * a)) => 2.5 + 0im)
    u0 = correlation_u0(c, u_end)
    # State is ⟨a · a'_anc⟩; undo_ancilla normal-orders to ⟨a'·a + 1⟩ = 3.5.
    @test only(values(u0)) ≈ 3.5
end

@testset "Damped cavity: τ-integration matches exp(-(iωc + κ/2)τ)" begin
    # End-to-end check of the τ-ODE pipeline: `to_system(c)`, `correlation_u0`,
    # `correlation_p0`, ODE solve. The closed form for a damped vacuum cavity
    # is ⟨a(τ)·a†(0)⟩_ss = exp(-(iωc + κ/2)τ).
    hc = FockSpace(:cavity)
    a = Destroy(hc, :a)
    @variables ωc::Real κ::Real
    H = ωc * a' * a
    he = meanfield(a' * a, H, [a]; rates = [κ])
    complete!(he)

    c = CorrelationFunction(a, a', he)
    ps = (ωc, κ)
    pn = ComplexF64[1.0, 0.5]

    # u_end at vacuum.
    u_end = Dict(SymbolicUtils.unwrap(average(a' * a)) => 0.0 + 0im)
    u0 = correlation_u0(c, u_end)
    p0_dict = correlation_p0(c, u_end, collect(ps) .=> pn)
    @named csys = to_system(c)
    csys_c = mtkcompile(csys)
    prob = ODEProblem(csys_c, merge(u0, p0_dict), (0.0, 30.0))
    sol = solve(prob, Tsit5(); abstol = 1e-10, reltol = 1e-10)
    @test sol.retcode == ReturnCode.Success

    # The first (and only) unknown of `csys_c` is ⟨a·a'_anc⟩(τ).
    g_var = first(unknowns(csys_c))
    τs = [0.5, 1.0, 2.0, 4.0]
    g_num = [sol(τ; idxs = g_var) for τ in τs]
    g_an = @. exp(-(im * pn[1] + pn[2] / 2) * τs)
    @test maximum(abs.(g_num .- g_an)) < 1e-7

    # g(τ → ∞) → 0 (vacuum, no asymptote). κ/2 = 0.25 ⇒ exp(-7.5) ≈ 5.5e-4.
    @test abs(sol(30.0; idxs = g_var)) < 1e-3
end

@testset "Mollow correlation: setup + ODE solve + CorrelationFunction shapes" begin
    h = NLevelSpace(:atom, (:g, :e))
    @variables Δ::Real Ω::Real γ::Real
    σ(i, j) = Transition(h, :σ, i, j)
    H = Δ * σ(:e, :e) + Ω * (σ(:g, :e) + σ(:e, :g))
    J = [σ(:g, :e)]
    eqs = meanfield([σ(:e, :g), σ(:e, :e)], H, J; rates = [γ])
    complete!(eqs)

    ps = (Δ, Ω, γ)
    p0 = ComplexF64[20.0, 5.0, 1.0]
    u0 = zeros(ComplexF64, length(eqs.equations))
    @named sys = to_system(eqs)
    sys_c = mtkcompile(sys)
    dict = merge(Dict(unknowns(sys_c) .=> u0), Dict(ps .=> p0))
    prob = ODEProblem(sys_c, dict, (0.0, 20.0))
    sol = solve(prob, Tsit5())
    @test sol.retcode == ReturnCode.Success

    # Physicality: excited-state population is Hermitian (real) and lives in [0, 1].
    assert_real(sol, σ(:e, :e), eqs)
    assert_population(sol, σ(:e, :e), eqs)

    # CorrelationFunction in both regimes; non-SS needs the seed-state
    # coupling so it has at least as many equations.
    c_ss = CorrelationFunction(σ(:e, :g), σ(:g, :e), eqs; steady_state = true)
    @test c_ss isa CorrelationFunction
    @test length(c_ss.eqs.equations) >= 1

    c_nss = CorrelationFunction(σ(:e, :g), σ(:g, :e), eqs; steady_state = false)
    @test c_nss isa CorrelationFunction
    @test length(c_nss.eqs.equations) >= length(c_ss.eqs.equations)
end
