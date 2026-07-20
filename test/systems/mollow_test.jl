using QuantumCumulants
using ModelingToolkitBase
using Symbolics: @variables
using OrdinaryDiffEqTsit5: ODEProblem, solve, Tsit5
using Test

@testset "mollow" begin
    @variables Δ Ω γ
    h = NLevelSpace(:atom, (:g, :e))
    σ(α, β) = Transition(h, :σ, α, β)
    H = Δ * σ(:e, :e) + Ω * (σ(:e, :g) + σ(:g, :e))
    eqs = meanfield([σ(:e, :g), σ(:e, :e)], H, [σ(:g, :e)]; rates = [γ])
    complete!(eqs)
    @test length(eqs.equations) == 2
    sys_c = mtkcompile(System(eqs; name = :mollow_sys))
    @test length(unknowns(sys_c)) == 2
    u0 = zeros(ComplexF64, length(eqs.equations))
    init = initial_values(eqs, u0)
    ps = Dict(Δ => 1.0, Ω => 0.5, γ => 0.1)
    prob = ODEProblem(sys_c, merge(init, ps), (0.0, 80.0))
    sol = solve(prob, Tsit5(); reltol = 1.0e-10, abstol = 1.0e-10)
    @test isapprox(
        real(get_solution(sol, σ(:e, :e), eqs).(sol.t[end])),
        0.16491193049237177; rtol = 1.0e-5,
    )
    @test isapprox(
        get_solution(sol, σ(:e, :g), eqs).(sol.t[end]),
        -0.33088938061453155 + 0.016559911920571824im; rtol = 1.0e-5,
    )
    c = CorrelationFunction(σ(:e, :g), σ(:g, :e), eqs; steady_state = true)
    @test length(c.eqs.equations) == 3
    ps_tup = (Δ, Ω, γ)
    p0 = (0.0, 2.0, 1.0)
    sys = mtkcompile(System(eqs; name = :mollow_steady))
    u0_ss = zeros(ComplexF64, length(unknowns(sys)))
    dict_ss = merge(Dict(unknowns(sys) .=> u0_ss), Dict(ps_tup .=> p0))
    prob_ss = ODEProblem(sys, dict_ss, (0.0, 20.0))
    sol_ss = solve(prob_ss, Tsit5(); reltol = 1.0e-10, abstol = 1.0e-10)
    S = Spectrum(c, ps_tup)
    s_vals = S([0.0, 2.0, -2.0], sol_ss.u[end], p0)
    @test isapprox(s_vals[1], 1.0257952485186328; rtol = 1.0e-5)
    @test isapprox(s_vals[2], 0.14359486828829415; rtol = 1.0e-5)
    # spectrum symmetry S(ω) = S(-ω)
    @test isapprox(s_vals[2], s_vals[3]; rtol = 1.0e-10)
end
