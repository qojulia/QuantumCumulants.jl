using QuantumCumulants
using Symbolics: Symbolics, @variables
using SymbolicUtils
using ModelingToolkitBase: ModelingToolkitBase, @named, mtkcompile, ODEProblem, unknowns
using OrdinaryDiffEq: Tsit5, solve, ReturnCode
using LinearAlgebra: I
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

@testset "correlation_u0: invalid u_end type throws" begin
    hc = FockSpace(:cavity)
    a = Destroy(hc, :a)
    @variables ωc::Real κ::Real
    H = ωc * a' * a
    he = meanfield(a' * a, H, [a]; rates = [κ])
    complete!(he)
    c = CorrelationFunction(a, a', he)
    # `u_end` must be a Dict (avg => value) or a Vector aligned with the states.
    @test_throws ArgumentError correlation_u0(c, 5)
end

@testset "Damped cavity: Spectrum matches analytic Lorentzian" begin
    # For H = ωc·a†a, J = [a] with rate κ, vacuum bath:
    #   ⟨a(τ)·a†(0)⟩_ss = exp(-(iωc + κ/2)τ), ⟨a·a†⟩_ss = 1.
    #   S(ω) = 2 Re ∫₀^∞ ⟨a(τ)a†(0)⟩ exp(-iωτ) dτ = κ / ((ω + ωc)² + (κ/2)²).
    hc = FockSpace(:cavity)
    a = Destroy(hc, :a)
    @variables ωc::Real κ::Real
    H = ωc * a' * a
    he = meanfield(a' * a, H, [a]; rates = [κ])
    complete!(he)

    c = CorrelationFunction(a, a', he)
    @test length(c.eqs.equations) == 1

    # Vacuum steady state: ⟨a†a⟩_ss = 0 ⇒ ⟨a·a†⟩_ss = 1.
    u_end = Dict(SymbolicUtils.unwrap(average(a' * a)) => 0.0 + 0im)
    u0 = correlation_u0(c, u_end)
    @test only(values(u0)) ≈ 1.0

    ωc_num = 1.0
    κ_num = 0.5
    ps = (ωc, κ)
    pn = ComplexF64[ωc_num, κ_num]
    ω = collect(range(-2π, 2π; length = 81))
    S = Spectrum(c, ps)
    s_num = S(ω, u_end, pn)
    s_an = @. κ_num / ((ω + ωc_num)^2 + (κ_num / 2)^2)
    @test maximum(abs.(s_num .- s_an)) < 1.0e-10
end

@testset "Spectrum: anti-linear (conjugate) augmentation, Hermitian op2" begin
    # The cross-correlation ⟨a(τ)·σ₂₂(0)⟩ pairs the cavity field with the atomic
    # population σ₂₂, which is Hermitian. Its conjugate ⟨a†(τ)·σ₂₂⟩ folds onto the
    # state (the correlation closes with get_adjoints=false), so the spectrum kernel
    # classifies that leaf as anti-linear and assembles the 2n augmented system
    #   M = [A_lin A_alin; conj(A_alin) conj(A_lin)]
    # instead of the bare n-dimensional A_lin. A plain Jaynes-Cummings model suffices:
    # the σ₁₂↔σ₂₁ coupling under a Hermitian op2 is what triggers the conjugate column,
    # no squeezing required.
    hc = FockSpace(:cavity)
    ha = NLevelSpace(:atom, 2)
    h = hc ⊗ ha
    a = Destroy(h, :a)
    σ(i, j) = Transition(h, :σ, i, j)
    @variables Δ::Real ωa::Real g::Real κ::Real γ::Real
    H = Δ * a' * a + ωa * σ(2, 2) + g * (a' * σ(1, 2) + a * σ(2, 1))
    he = complete(meanfield([a' * a, σ(2, 2)], H, [a, σ(1, 2)]; rates = [κ, γ], order = 2))

    c = CorrelationFunction(a, σ(2, 2), he)
    n = length(c.eqs.equations)
    ps = (Δ, ωa, g, κ, γ)
    p0 = ComplexF64[0.0, 0.0, 1.0, 1.0, 0.5]
    u_end = Dict(SymbolicUtils.unwrap(s) => ComplexF64(0.05k) for (k, s) in enumerate(he.states))
    S = Spectrum(c, ps)

    # The kernel takes the anti-linear branch: the system is augmented to 2n, and the
    # anti-linear block is actually populated (a no-op augmentation would leave it zero).
    M, rhs_b, dim = QuantumCumulants._spectrum_kernel(S, u_end, p0)
    @test dim == 2n
    @test size(M) == (2n, 2n)
    @test any(!iszero, M[1:n, (n + 1):(2n)])

    # The conjugate coupling is load-bearing: zeroing the off-diagonal (anti-linear)
    # blocks of M changes the spectrum. A dropped/no-op augmentation, or an A_alin that
    # never reaches the solve, would leave the two spectra equal.
    ωp = [0.4, 0.9, 1.5]
    s_full = [2 * real(((im * ω) * I(2n) - M) \ rhs_b)[1] for ω in ωp]
    M0 = copy(M)
    M0[1:n, (n + 1):(2n)] .= 0
    M0[(n + 1):(2n), 1:n] .= 0
    s_lin = [2 * real(((im * ω) * I(2n) - M0) \ rhs_b)[1] for ω in ωp]
    @test maximum(abs.(s_full .- s_lin)) > 1.0e-3

    # The resulting symmetric power spectrum is real and finite.
    s = S(collect(range(-3.0, 3.0; length = 41)), u_end, p0)
    @test eltype(s) == Float64
    @test all(isfinite, s)
end

@testset "Driven cavity: ODE conjugate substitution is non-folding" begin
    # d⟨a'a⟩/dt = -κ⟨a'a⟩ + iΩ(⟨a⟩ - ⟨a'⟩); steady state
    # ⟨a⟩_ss = -iΩ/(iωc + κ/2), ⟨a'a⟩_ss = |⟨a⟩_ss|².
    hc = FockSpace(:cavity)
    a = Destroy(hc, :a)
    @variables ωc::Real κ::Real Ω::Real
    H = ωc * a' * a + Ω * (a + a')
    he = meanfield([a, a'a], H, [a]; rates = [κ])
    complete!(he)

    @named sys = System(he)
    sys_c = mtkcompile(sys)
    ωc_n, κ_n, Ω_n = 1.0, 0.5, 0.3
    ps = [ωc, κ, Ω]
    pn = [ωc_n, κ_n, Ω_n]
    u0 = Dict(unknowns(sys_c) .=> zeros(ComplexF64, length(unknowns(sys_c))))
    sol = solve(
        ODEProblem(sys_c, merge(u0, Dict(ps .=> pn)), (0.0, 80.0)),
        Tsit5(); abstol = 1.0e-10, reltol = 1.0e-10,
    )
    @test sol.retcode == ReturnCode.Success

    a_ss_an = -im * Ω_n / (im * ωc_n + κ_n / 2)
    n_ss_an = abs2(a_ss_an)
    @test abs(get_solution(sol, a, he)(80.0) - a_ss_an) < 1.0e-7
    @test abs(get_solution(sol, a' * a, he)(80.0) - n_ss_an) < 1.0e-7
end

@testset "_lookup_avg resolves compound (normal-ordered) expressions" begin
    hc = FockSpace(:cavity)
    a = Destroy(hc, :a)
    @variables ωc::Real κ::Real
    H = ωc * a' * a
    he = meanfield(a' * a, H, [a]; rates = [κ])
    complete!(he)

    c = CorrelationFunction(a, a', he)
    u_end = Dict(SymbolicUtils.unwrap(average(a' * a)) => 2.5 + 0im)
    u0 = correlation_u0(c, u_end)
    # ⟨a · a'⟩ normal-orders to ⟨a'·a + 1⟩ = 3.5.
    @test only(values(u0)) ≈ 3.5
end

@testset "Damped cavity: τ-integration matches exp(-(iωc + κ/2)τ)" begin
    # Damped vacuum cavity: ⟨a(τ)·a†(0)⟩_ss = exp(-(iωc + κ/2)τ).
    hc = FockSpace(:cavity)
    a = Destroy(hc, :a)
    @variables ωc::Real κ::Real
    H = ωc * a' * a
    he = meanfield(a' * a, H, [a]; rates = [κ])
    complete!(he)

    c = CorrelationFunction(a, a', he)
    ps = (ωc, κ)
    pn = ComplexF64[1.0, 0.5]

    u_end = Dict(SymbolicUtils.unwrap(average(a' * a)) => 0.0 + 0im)
    u0 = correlation_u0(c, u_end)
    p0_dict = correlation_p0(c, u_end, collect(ps) .=> pn)
    @named csys = System(c)
    csys_c = mtkcompile(csys)
    prob = ODEProblem(csys_c, merge(u0, p0_dict), (0.0, 30.0))
    sol = solve(prob, Tsit5(); abstol = 1.0e-10, reltol = 1.0e-10)
    @test sol.retcode == ReturnCode.Success

    g_var = first(unknowns(csys_c))
    τs = [0.5, 1.0, 2.0, 4.0]
    g_num = [sol(τ; idxs = g_var) for τ in τs]
    g_an = @. exp(-(im * pn[1] + pn[2] / 2) * τs)
    @test maximum(abs.(g_num .- g_an)) < 1.0e-7

    # g(τ → ∞) → 0: κ/2 = 0.25 ⇒ exp(-7.5) ≈ 5.5e-4.
    @test abs(sol(30.0; idxs = g_var)) < 1.0e-3
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
    @named sys = System(eqs)
    sys_c = mtkcompile(sys)
    dict = merge(Dict(unknowns(sys_c) .=> u0), Dict(ps .=> p0))
    prob = ODEProblem(sys_c, dict, (0.0, 20.0))
    sol = solve(prob, Tsit5())
    @test sol.retcode == ReturnCode.Success

    # Excited-state population is real and lies in [0, 1].
    assert_real(sol, σ(:e, :e), eqs)
    assert_population(sol, σ(:e, :e), eqs)

    c_ss = CorrelationFunction(σ(:e, :g), σ(:g, :e), eqs; steady_state = true)
    @test c_ss isa CorrelationFunction
    @test length(c_ss.eqs.equations) >= 1

    c_nss = CorrelationFunction(σ(:e, :g), σ(:g, :e), eqs; steady_state = false)
    @test c_nss isa CorrelationFunction
    @test length(c_nss.eqs.equations) >= length(c_ss.eqs.equations)
end
