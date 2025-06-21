using QuantumCumulants
using Test
using OrdinaryDiffEq, ModelingToolkit

@testset "higher-order" begin

    # Parameters
    tspan = range(0.0, 10.0, length = 101)
    ω = range(-2pi, 2pi, length = 201)

    # Define parameters
    @cnumbers Δ g γ κ ν

    # Define hilbert space
    hf = FockSpace(:cavity)
    ha = NLevelSpace(:atom, (:g, :e))
    h = hf ⊗ ha

    # Define the fundamental operators
    a = Destroy(h, :a)
    s = Transition(h, :σ, :g, :e)

    # Hamiltonian
    H = Δ*a'*a + g*(a'*s + a*s')

    # Collapse operators
    J = [a, s, s']
    rates = [κ, γ, ν]


    # Custom filter function -- include only phase-invaraint terms
    ϕ(x) = 0
    ϕ(::Destroy) = -1
    ϕ(::Create) = 1
    function ϕ(t::Transition)
        if (t.i==:e && t.j==:g)
            1
        elseif (t.i==:g && t.j==:e)
            -1
        else
            0
        end
    end
    ϕ(avg::Average) = ϕ(avg.arguments[1])
    function ϕ(t::QuantumCumulants.QMul)
        p = 0
        for arg in t.args_nc
            p += ϕ(arg)
        end
        return p
    end
    phase_invariant(x) = iszero(ϕ(x))


    # Derive equations
    he_n = meanfield(a'*a, H, J; rates = rates)

    ## Fourth order
    he4 = complete(he_n; order = 4, filter_func = phase_invariant)


    # Numerical solution
    ps = (Δ, g, γ, κ, ν)
    p0 = (0.5, 1.5, 0.25, 1, 4)
    @named sys4 = ODESystem(he4)
    u0 = zeros(ComplexF64, length(he4))
    prob4 = ODEProblem(sys4, u0, (0.0, tspan[end]), ps .=> p0)
    sol4 = solve(prob4, RK4())

    n4 = real.(getindex.(sol4.(tspan), 1))

    # Correlation function and spectrum
    c4 = CorrelationFunction(a', a, he4; steady_state = true, filter_func = phase_invariant)
    S4 = Spectrum(c4, ps)
    s4 = S4(ω, sol4.u[end], p0)

    ## Sixth order
    he6 = complete(he_n; order = 6, filter_func = phase_invariant)


    # Numerical solution
    ps = (Δ, g, γ, κ, ν)
    p0 = (0.5, 1.5, 0.25, 1, 4)
    @named sys6 = ODESystem(he6)
    u0 = zeros(ComplexF64, length(he6))
    prob6 = ODEProblem(sys6, u0, (0.0, tspan[end]), ps .=> p0)
    sol6 = solve(prob6, RK4())

    n6 = real.(getindex.(sol6.(tspan), 1))

    # Correlation function and spectrum
    c6 = CorrelationFunction(a', a, he6; steady_state = true, filter_func = phase_invariant)
    S6 = Spectrum(c6, ps)
    s6 = S6(ω, sol6.u[end], p0)

    ## Test
    @test maximum(abs.(n6 .- n4)) < 0.05
    @test maximum(abs.(s6 .- s4)) < 0.2


end # testset
