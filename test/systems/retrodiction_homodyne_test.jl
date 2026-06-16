using QuantumCumulants
using ModelingToolkitBase
using Symbolics: @variables
using OrdinaryDiffEqTsit5: ODEProblem, solve, Tsit5
using Test

@testset "retrodiction_homodyne (forward, noise)" begin
    using StochasticDiffEqLowOrder: SDEProblem, EM, ReturnCode
    using DiffEqNoiseProcess: RealWienerProcess
    import Random
    h = PhaseSpace(:motion)
    @qnumbers x::Position(h)
    p = Momentum(h, :p)
    @variables Ω Γ η s
    m = 1
    a = (x + 1im * p) * s
    H = p^2 / (2m) + 0.5m * Ω^2 * x^2
    eqs = meanfield(
        [x, p, x * x, x * p, p * p], H, [a];
        rates = [Γ], efficiencies = [η], order = 2,
    )
    @test length(eqs.equations) == 5
    sys_c = mtkcompile(System(eqs; name = :ret_sys))
    @test length(unknowns(sys_c)) == 5
    u0 = ComplexF64[5.0, 0.0, 5 + 25.0, 0.5im, 5.0]
    init = initial_values(eqs, u0)
    pmap = Dict(Ω => 1.0, Γ => 1 / 6, η => 0.5, s => 1 / √2)
    Random.seed!(11)
    noise = RealWienerProcess(0.0, 0.0)
    prob_st = SDEProblem(sys_c, merge(init, pmap), (0.0, 2.0); noise = noise)
    sol = solve(prob_st, EM(); dt = 1.0e-4)
    x_traj = real.(get_solution(sol, x, eqs).(sol.t))
    @test all(isfinite, x_traj)
    @test maximum(abs, x_traj) < 1.0e6
end

# Forward drift only (no measurement noise), solved as an ODE over the
# PhaseSpace / Position / Momentum operator path.
@testset "retrodiction_homodyne (forward drift only)" begin
    h = PhaseSpace(:motion)
    @qnumbers x::Position(h)
    p = Momentum(h, :p)
    @variables Ω Γ s
    m = 1
    a = (x + 1im * p) * s
    H = p^2 / (2m) + 0.5m * Ω^2 * x^2
    eqs = meanfield(
        [x, p, x * x, x * p, p * p], H, [a];
        rates = [Γ], order = 2,
    )
    sys_c = mtkcompile(System(eqs; name = :ret_drift_sys))
    u0 = ComplexF64[5.0, 0.0, 5 + 25.0, 0.5im, 5.0]
    init = initial_values(eqs, u0)
    pmap = Dict(Ω => 1.0, Γ => 1 / 6, s => 1 / √2)
    prob = ODEProblem(sys_c, merge(init, pmap), (0.0, 5.0))
    sol = solve(prob, Tsit5(); reltol = 1.0e-9, abstol = 1.0e-9)
    @test isapprox(
        real(get_solution(sol, x, eqs).(sol.t[end])),
        0.935008189535323; rtol = 1.0e-5,
    )
    @test isapprox(
        real(get_solution(sol, p, eqs).(sol.t[end])),
        3.1608092157050347; rtol = 1.0e-5,
    )
    @test isapprox(
        real(get_solution(sol, x * x, eqs).(sol.t[end])),
        3.3299322539872325; rtol = 1.0e-5,
    )
end
