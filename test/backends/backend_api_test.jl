using QuantumCumulants
using Symbolics: Symbolics, @variables
using SciMLBase: SciMLBase, ODEProblem, ODEFunction
using OrdinaryDiffEqLowOrderRK: RK4, solve
using Test

# Two-spin order-1 model (small; the API surface is what is under test here).
Ns = 2
hs = ⊗([PauliSpace(Symbol(:s, i)) for i in 1:Ns]...)
sz(i) = Pauli(hs, :σ, 3, i)
sx(i) = Pauli(hs, :σ, 1, i)
sy(i) = Pauli(hs, :σ, 2, i)
sm(i) = (sx(i) - 1im * sy(i)) / 2
@variables J hx γ
Hs = -J * sz(1) * sz(2) - hx * (sx(1) + sx(2))
eqs_s = meanfield([sz(i) for i in 1:Ns], Hs, [sm(i) for i in 1:Ns]; rates = [γ, γ], order = 1)
complete!(eqs_s)
ps = Dict(J => 1.0, hx => 1.0, γ => 0.2)
u0s = zeros(ComplexF64, length(eqs_s.states))

@testset "AutoBackend routes kernel first, sharded on NonPolynomialDriftError" begin
    prob = ODEProblem(eqs_s, u0s, (0.0, 1.0), ps)                     # default AutoBackend
    @test prob.p isa QuantumCumulants.KernelParameters
    eqs_np = modify_equations(eqs_s, (op, d) -> d + exp(average(op)))
    prob_np = @test_logs (:info, r"falling back") match_mode = :any ODEProblem(
        eqs_np, u0s, (0.0, 0.1), ps
    )
    @test !(prob_np.p isa QuantumCumulants.KernelParameters)
    @test solve(prob_np, RK4()).retcode == SciMLBase.ReturnCode.Success
    # every OTHER lowering error stays hard through AutoBackend
    @test_throws TimeDependentCoefficientError ODEProblem(
        modify_equations(eqs_s, (op, d) -> cos(eqs_s.iv) * d), u0s, (0.0, 1.0), ps
    )
end

@testset "guard methods" begin
    # parameters are required at construction
    @test_throws ArgumentError ODEProblem(eqs_s, u0s, (0.0, 1.0))
    @test_throws ArgumentError ODEFunction(eqs_s)
end
