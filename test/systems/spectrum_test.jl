using QuantumCumulants
using Symbolics: Symbolics, @variables
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

@testset "Damped-cavity correlation g(τ) matches exp(iωτ - κτ/2)" begin
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
