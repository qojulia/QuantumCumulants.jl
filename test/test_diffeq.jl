using QuantumCumulants
using OrdinaryDiffEq
using SteadyStateDiffEq
using ModelingToolkit
using Test

@testset "diffeq" begin

    hf = FockSpace(:cavity)
    ha = NLevelSpace(:atom, (:g, :e))
    h = tensor(hf, ha)

    @qnumbers a::Destroy(h) σ::Transition(h, :g, :e)

    # Single-atom laser
    @cnumbers Δ g κ γ ν

    H = Δ*a'*a + g*(a'*σ + σ'*a)
    J = [a, σ, σ']
    he_avg = meanfield([a'*a, σ'*σ, a*σ'], H, J; rates = [κ, γ, ν])

    he_exp = cumulant_expansion(he_avg, 2)
    @test isequal(
        he_exp.equations,
        meanfield([a'*a, σ'*σ, a*σ'], H, J; rates = [κ, γ, ν], order = 2).equations,
    )

    ps = [Δ, g, κ, γ, ν]
    missed = find_missing(he_exp)
    @test !any(p ∈ Set(missed) for p in ps)

    # Exploit phase invariance
    subs = Dict([missed; QuantumCumulants._conj.(missed)] .=> 0)
    he_nophase = substitute(he_exp, subs)
    @test isempty(find_missing(he_nophase))

    @named sys = System(he_nophase)

    # Numerical solution
    p0 = ps .=> [0.0, 0.5, 1.0, 0.1, 0.9]
    u0 = zeros(ComplexF64, 3)
    tmax = 10.0

    prob = ODEProblem(sys, u0, (0.0, tmax), p0)
    sol = solve(prob, RK4())
    n = sol[average(a'*a)]
    pe = getindex.(sol.u, 2)

    @test all(iszero.(imag.(n)))
    @test all(iszero.(imag.(pe)))
    @test all(real.(n) .>= 0.0)
    @test all(1.0 .>= real.(pe) .>= 0.0)

    prob_ss = SteadyStateProblem(prob)
    sol_ss = solve(prob_ss, DynamicSS(Tsit5()), abstol = 1e-6, reltol = 1e-6)

    @test isapprox(sol_ss[a'*a], sol[a'a][end], atol = 1e-4)

end # testset
