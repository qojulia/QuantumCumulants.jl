using QuantumCumulants
using SymbolicUtils: SymbolicUtils
using Symbolics: Symbolics, @variables
using SciMLBase: SciMLBase, ODEProblem
using Test

# Two-spin transverse-field Ising model, order 1: small fixture for the error taxonomy.
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
u0 = zeros(ComplexF64, length(eqs_s.states))
function mkprob(e; kw...)
    return ODEProblem(
        e, zeros(ComplexF64, length(e.states)), (0.0, 1.0), ps;
        backend = KernelBackend(), kw...,
    )
end

@testset "kernel lowering error taxonomy" begin
    @test_throws NonPolynomialDriftError mkprob(
        modify_equations(eqs_s, (op, d) -> d + exp(average(op)))
    )
    @test_throws TimeDependentCoefficientError mkprob(
        modify_equations(eqs_s, (op, d) -> cos(eqs_s.iv) * d)
    )
    imvar = Symbolics.variable(:im)
    @test_throws ImParameterCollisionError mkprob(
        modify_equations(eqs_s, (op, d) -> imvar * d)
    )
    @test_throws ImParameterCollisionError ODEProblem(
        eqs_s, u0, (0.0, 1.0), merge(ps, Dict(imvar => 2.0)); backend = KernelBackend()
    )
    # order-2 hierarchy left uncompleted: RHS moments missing from the states
    eqs_open = meanfield([sz(1)], Hs, [sm(i) for i in 1:Ns]; rates = [γ, γ], order = 2)
    @test_throws UnresolvedMomentError mkprob(eqs_open)
end

@testset "parameter and eltype contract" begin
    prob = mkprob(eqs_s)
    # a plain vector in prob.p is a typed error at the first RHS call
    prob_bad = SciMLBase.remake(prob, p = zeros(3))
    @test_throws ArgumentError prob_bad.f(u0, u0, prob_bad.p, 0.0)
    # dual/eltype guard: the kernel is compiled for ComplexF64 states
    @test_throws ArgumentError prob.f(u0, zeros(Float64, length(u0)), prob.p, 0.0)
end
