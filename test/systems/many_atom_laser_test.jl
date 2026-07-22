using QuantumCumulants
using ModelingToolkitBase
using Symbolics: @variables
using OrdinaryDiffEqTsit5: ODEProblem, solve, Tsit5
using Test

@testset "many-atom-laser" begin
    N = 2
    @variables κ g Γ23 Γ13 Γ12 Ω Δc Δ3
    hf = FockSpace(:cavity)
    ha = ⊗([NLevelSpace(Symbol(:atom, i), 3) for i in 1:N]...)
    h = hf ⊗ ha
    a = Destroy(h, :a)
    σ(i, j, k) = Transition(h, Symbol("σ_{$k}"), i, j, k + 1)
    H = -Δc * a' * a +
        sum(g * (a' * σ(1, 2, i) + a * σ(2, 1, i)) for i in 1:N) +
        sum(Ω * (σ(3, 1, i) + σ(1, 3, i)) for i in 1:N) -
        sum(Δ3 * σ(3, 3, i) for i in 1:N)
    J = [a; [σ(1, 2, i) for i in 1:N]; [σ(1, 3, i) for i in 1:N]; [σ(2, 3, i) for i in 1:N]]
    rates = [κ; [Γ12 for _ in 1:N]; [Γ13 for _ in 1:N]; [Γ23 for _ in 1:N]]
    ops = [a' * a, σ(2, 2, 1), σ(3, 3, 1)]
    eqs = meanfield(ops, H, J; rates = rates, order = 2)
    complete!(eqs)
    @test length(eqs.equations) == 63
    sys_c = mtkcompile(System(eqs; name = :mal_sys))
    @test length(unknowns(sys_c)) == 63
    u0 = zeros(ComplexF64, length(eqs.equations))
    init = initial_values(eqs, u0)
    ps = Dict(
        κ => 0.5, g => 2.0, Γ23 => 20.0, Γ13 => 2.0, Γ12 => 1.0,
        Ω => 10.0, Δc => 0.0, Δ3 => 0.0,
    )
    prob = ODEProblem(sys_c, merge(init, ps), (0.0, 10.0))
    sol = solve(prob, Tsit5(); reltol = 1.0e-9, abstol = 1.0e-9)
    @test isapprox(
        real(get_solution(sol, a' * a, eqs).(sol.t[end])),
        12.908870121878037; rtol = 1.0e-4,
    )
    @test isapprox(
        real(get_solution(sol, σ(2, 2, 1), eqs).(sol.t[end])),
        0.4479907519944381; rtol = 1.0e-4,
    )
end
