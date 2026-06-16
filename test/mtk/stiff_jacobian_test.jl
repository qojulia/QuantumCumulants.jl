using QuantumCumulants
using ModelingToolkitBase: System, mtkcompile, ODEProblem
using OrdinaryDiffEqRosenbrock: OrdinaryDiffEqRosenbrock, Rosenbrock23, solve
using Symbolics: @variables
using Test

# A stiff solver builds a Jacobian ∂f/∂u; with the Number-symtype complex unknowns the
# state vector is ComplexF64. ForwardDiff (the default AD) cannot differentiate complex
# states, a pre-existing limitation shared by the old Real-symtype + complex-u0 approach,
# so a stiff solve must use a complex-compatible AD (finite differencing). This guards
# that Number unknowns integrate correctly through a stiff (Jacobian-based) solver.
@testset "stiff solver, complex state, Number unknowns" begin
    h = FockSpace(:c)
    a = Destroy(h, :a)
    @variables Δ κ
    eqs = complete(meanfield([a], Δ * a' * a, [a]; rates = [κ]))
    sys = mtkcompile(System(eqs; name = :s))
    u0 = initial_values(eqs; defaults = Dict(average(a) => 1.0 + 0.0im))
    prob = ODEProblem(sys, vcat(collect(u0), [Δ => 2.0, κ => 1.0]), (0.0, 3.0))
    sol = solve(prob, Rosenbrock23(autodiff = OrdinaryDiffEqRosenbrock.AutoFiniteDiff()); reltol = 1.0e-9, abstol = 1.0e-9)
    g = get_solution(sol, average(a), eqs)
    exact = exp((-1im * 2.0 - 0.5) * 3.0)
    @test isapprox(g(3.0), exact; atol = 1.0e-4)
end
