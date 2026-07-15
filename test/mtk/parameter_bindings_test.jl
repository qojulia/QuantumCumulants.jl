using QuantumCumulants
using Symbolics: @variables
using ModelingToolkitBase: mtkcompile, ODEProblem, parameters, unknowns
using OrdinaryDiffEqTsit5: Tsit5, solve
using Test

# `System(eqs; name, kwargs...)` forwards extra keyword arguments to the underlying
# `ModelingToolkitBase.System`. The `bindings` kwarg marks a parameter as *derived*
# from others, so it is computed from the supplied values rather than set independently.
# Here the collective drive is `N * gΩ`; binding `gΩ` to `2.5 / N` makes that product
# intensive (independent of `N`), so both the compiled parameter set and the dynamics
# are `N`-independent. This mirrors the effective coupling in the unique-squeezing example.
@testset "System forwards bindings: derived parameter is intensive in N" begin
    h = FockSpace(:cavity)
    a = Destroy(h, :a)
    @variables Δ::Real κ::Real N::Real gΩ::Real

    H = Δ * a' * a + N * gΩ * (a + a')
    eqs = meanfield([a], H, [a]; rates = [κ], order = 1)

    sys = mtkcompile(System(eqs; name = :bound, bindings = [gΩ => 2.5 / N]))

    # gΩ is derived, so it is not a tunable parameter of the compiled system.
    pnames = string.(parameters(sys))
    @test !("gΩ" in pnames)
    @test "N" in pnames

    Δv, κv = 1.0, 0.5
    aend = ComplexF64[]
    for Nv in (1.0, 10.0, 137.0)
        u0 = Dict(unknowns(sys) .=> ComplexF64[0.0])
        prob = ODEProblem(sys, merge(u0, Dict(Δ => Δv, κ => κv, N => Nv)), (0.0, 200.0))
        sol = solve(prob, Tsit5(); reltol = 1.0e-10, abstol = 1.0e-10)
        push!(aend, get_solution(sol, a, eqs)(200.0))
    end

    # Intensive N*gΩ ⇒ steady state independent of N.
    @test aend[1] ≈ aend[2] atol = 1.0e-8
    @test aend[2] ≈ aend[3] atol = 1.0e-8

    # Magnitude matches the analytic steady state with N*gΩ = 2.5 (sign-convention free).
    @test abs(aend[1]) ≈ 2.5 / sqrt(Δv^2 + (κv / 2)^2) atol = 1.0e-6
end
