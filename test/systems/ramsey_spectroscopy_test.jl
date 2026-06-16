using QuantumCumulants
using ModelingToolkitBase
using Symbolics: @variables
using OrdinaryDiffEqTsit5: ODEProblem, solve, Tsit5
using Test

@testset "ramsey_spectroscopy" begin
    @variables Δ Ω Γ ν
    @register_symbolic _ramsey_f(tt)
    h = NLevelSpace(:atom, 2)
    σ(i, j) = Transition(h, :σ, i, j)
    eqs_seed = meanfield(
        [σ(2, 2), σ(1, 2)], -Δ * σ(2, 2),
        [σ(1, 2), σ(2, 2)]; rates = [Γ, ν],
    )
    tv = eqs_seed.iv
    H = -Δ * σ(2, 2) + _ramsey_f(tv) * Ω / 2 * (σ(1, 2) + σ(2, 1))
    eqs = meanfield(
        [σ(2, 2), σ(1, 2)], H, [σ(1, 2), σ(2, 2)];
        rates = [Γ, ν], iv = tv,
    )
    complete!(eqs)
    @test length(eqs.equations) == 2
    sys_c = mtkcompile(System(eqs; name = :ramsey_sys))
    @test length(unknowns(sys_c)) == 2
    _ramsey_f(_) = 1.0
    u0 = zeros(ComplexF64, length(eqs.equations))
    init = initial_values(eqs, u0)
    ps = Dict(Δ => 0.0, Ω => 1.0, Γ => 0.1, ν => 0.05)
    prob = ODEProblem(sys_c, merge(init, ps), (0.0, 20.0))
    sol = solve(prob, Tsit5(); reltol = 1.0e-10, abstol = 1.0e-10)
    @test isapprox(
        real(get_solution(sol, σ(2, 2), eqs).(sol.t[end])),
        0.45407705143934646; rtol = 1.0e-5,
    )
    @test isapprox(
        get_solution(sol, σ(1, 2), eqs).(sol.t[end]),
        0.0 - 0.12468144623895208im; rtol = 1.0e-5,
    )
end
