using QuantumCumulants
using ModelingToolkitBase
using Symbolics: @variables
using OrdinaryDiffEqTsit5: ODEProblem, solve, Tsit5
using Test

@testset "optomechanical-cooling" begin
    hc = FockSpace(:cavity); hm = FockSpace(:motion); h = hc ⊗ hm
    @qnumbers a::Destroy(h, 1) b::Destroy(h, 2)
    @variables Δ ωm E G κ
    H = -Δ * a' * a + ωm * b' * b + G * a' * a * (b + b') + E * (a + a')
    eqs = meanfield([a' * a, b' * b], H, [a]; rates = [κ], order = 2)
    eqs_c = complete!(deepcopy(eqs))
    @test length(eqs_c.equations) == 14
    sys_c = mtkcompile(System(eqs_c; name = :om_sys))
    @test length(unknowns(sys_c)) == 14
    u0 = zeros(ComplexF64, length(eqs_c.equations))
    init = initial_values(
        eqs_c, u0; defaults = Dict(average(b' * b) => 100.0 + 0im),
    )
    ps = Dict(Δ => -1.0, ωm => 1.0, E => 0.5, G => 0.05, κ => 0.1)
    prob = ODEProblem(sys_c, merge(init, ps), (0.0, 50.0))
    sol = solve(prob, Tsit5(); reltol = 1.0e-9, abstol = 1.0e-9)
    @test isapprox(
        real(get_solution(sol, a' * a, eqs_c).(sol.t[end])),
        0.28986377256026735; rtol = 1.0e-4,
    )
    @test isapprox(
        real(get_solution(sol, b' * b, eqs_c).(sol.t[end])),
        0.03417729886018334; rtol = 1.0e-3,
    )
end
