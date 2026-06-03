using QuantumCumulants
using Symbolics: Symbolics, @variables
using ModelingToolkitBase: ModelingToolkitBase, @named, mtkcompile, ODEProblem, unknowns
using OrdinaryDiffEq: Tsit5, solve, ReturnCode
using Test

@testset "Dicke model (Pauli): meanfield + complete + ODE" begin
    @variables Δ_ g κ η
    hf = FockSpace(:cavity)
    hs1 = PauliSpace(:spin1)
    hs2 = PauliSpace(:spin2)
    h = hf ⊗ hs1 ⊗ hs2
    a = Destroy(h, :a)
    σ(s, axis) = Pauli(h, Symbol(:σ, s), axis, s + 1)

    H = Δ_ * a' * a + g * (a' + a) * (σ(1, 1) + σ(2, 1)) + η * (a' + a)
    J = [a]
    rates = [κ]
    eq = meanfield(
        [σ(1, 3), σ(2, 3), σ(1, 3) * σ(2, 3)], H, J;
        rates = rates, order = 2
    )
    eqs = complete(eq)
    @test isempty(find_missing(eqs))

    ps = [Δ_, g, κ, η]
    @named sys = System(eqs)
    sys_c = mtkcompile(sys)
    u0 = Dict(unknowns(sys_c) .=> zeros(ComplexF64, length(unknowns(sys_c))))
    p0 = [0.5, 1.0, 1.25, 0.85]
    dict = merge(u0, Dict(ps .=> p0))
    prob = ODEProblem(sys_c, dict, (0.0, 0.5))
    sol = solve(prob, Tsit5())
    @test sol.retcode == ReturnCode.Success

    # Pauli-Z expectation is real and bounded in [-1, 1] for each spin.
    assert_real(sol, σ(1, 3), eqs)
    assert_bounded(sol, σ(1, 3), eqs, -1.0, 1.0)
    assert_real(sol, σ(2, 3), eqs)
    assert_bounded(sol, σ(2, 3), eqs, -1.0, 1.0)
end
