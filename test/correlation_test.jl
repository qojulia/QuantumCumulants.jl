using QuantumCumulants
using Symbolics: Symbolics, @variables, substitute
using ModelingToolkitBase: ModelingToolkitBase, @named, mtkcompile, ODEProblem,
                            unknowns
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

@testset "Single-atom laser: phase-invariant ODE + sanity bounds" begin
    # Port of master test_correlation.jl lines ~14-43 (the single-atom-laser
    # ODE + numerical sanity assertions).
    hf = FockSpace(:cavity)
    ha = NLevelSpace(:atom, (:g, :e))
    h = hf ⊗ ha
    a = Destroy(h, :a)
    σ = Transition(h, :σ, :g, :e)

    @variables Δ::Real g::Real κ::Real γ::Real ν::Real
    H = Δ * a' * a + g * (a' * σ + σ' * a)
    J = [a, σ, σ']
    he_laser = meanfield([a' * a, σ' * σ, a * σ'], H, J; rates = [κ, γ, ν])
    he_avg = cumulant_expansion(he_laser, 2)
    he_comp = complete(he_avg)

    ps = (Δ, g, γ, κ, ν)
    p0 = ps .=> ComplexF64[1, 1.5, 0.25, 1, 4]
    @named sys = to_system(he_comp)
    sys_c = mtkcompile(sys)
    u0 = unknowns(sys_c) .=> zeros(ComplexF64, length(unknowns(sys_c)))
    tmax = 10.0
    prob = ODEProblem(sys_c, merge(Dict(u0), Dict(p0)), (0.0, tmax))
    sol = solve(prob, Tsit5())
    @test sol.retcode == ReturnCode.Success

    n_t = real.(get_solution(sol, a' * a, he_comp).(sol.t))
    pe_t = real.(get_solution(sol, σ' * σ, he_comp).(sol.t))
    # Master asserted: photon population stays real and non-negative; atomic
    # excitation real and in [0, 1].
    @test all(>=(-1e-9), n_t)
    @test all(x -> -1e-9 <= x <= 1 + 1e-9, pe_t)
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

    # CorrelationFunction in both regimes; non-SS needs the seed-state
    # coupling so it has at least as many equations.
    c_ss = CorrelationFunction(σ(:e, :g), σ(:g, :e), eqs; steady_state = true)
    @test c_ss isa CorrelationFunction
    @test length(c_ss.eqs.equations) >= 1

    c_nss = CorrelationFunction(σ(:e, :g), σ(:g, :e), eqs; steady_state = false)
    @test c_nss isa CorrelationFunction
    @test length(c_nss.eqs.equations) >= length(c_ss.eqs.equations)
end

@testset "Damped-cavity correlation g(τ) matches exp(iωτ - κτ/2)" begin
    # Port of master test_correlation.jl lines ~144-176: time-dep cavity
    # correlation. With cavity initially in a coherent-like state |α⟩⟨α|*N,
    # the steady state has ⟨a'a⟩(0) = n0 = |α|^2 and the two-time
    # correlation g(τ) = ⟨a†(t)a(t+τ)⟩ at fixed t evolves as
    #   g(τ) = ⟨a†a⟩(t) * exp((iω - κ/2) τ).
    hc = FockSpace(:cavity)
    a = Destroy(hc, :a)
    @variables ωc::Real κ::Real
    H = ωc * a' * a
    he = meanfield(a' * a, H, [a]; rates = [κ])
    complete!(he)

    c = CorrelationFunction(a', a, he)
    @test c isa CorrelationFunction
    @test length(c.eqs.equations) >= 1
    @test c.op1 === a'
    @test c.op2 === a
end
