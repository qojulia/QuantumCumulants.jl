using QuantumCumulants
using SymbolicUtils: SymbolicUtils
using Symbolics: Symbolics, @variables
using SciMLBase: SciMLBase, ODEFunction, ODEProblem
using Test

const QC = QuantumCumulants

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

    # Error rendering is part of the typed diagnostic contract, not merely an implementation
    # detail of @test_throws.
    @test occursin("non-polynomial", sprint(showerror, NonPolynomialDriftError(2, :residual)))
    @test occursin("independent variable", sprint(showerror, TimeDependentCoefficientError(:c)))
    @test occursin("parameter named `im`", sprint(showerror, ImParameterCollisionError()))
    @test occursin("holomorphic-only", sprint(showerror, HolomorphicJacobianError()))
    @test occursin("does not resolve", sprint(showerror, UnresolvedMomentError(:m)))
end

@testset "parameter and eltype contract" begin
    prob = mkprob(eqs_s)
    # a plain vector in prob.p is a typed error at the first RHS call
    prob_bad = SciMLBase.remake(prob, p = zeros(3))
    @test_throws ArgumentError prob_bad.f(u0, u0, prob_bad.p, 0.0)
    @test_throws ArgumentError update_parameters!(prob_bad, Dict(J => 2.0))
    # dual/eltype guard: the kernel is compiled for ComplexF64 states
    @test_throws ArgumentError prob.f(u0, zeros(Float64, length(u0)), prob.p, 0.0)
    @test_throws ArgumentError update_parameters!(prob.f, Dict(J => 2.0))

    ir = QC._lower_moment_ir(eqs_s)
    @test_throws ArgumentError QC.kernel_pdict(ir.params, Dict())
    @test isempty(QC.kernel_pdict(ir.params, Dict(); strict = false))
end

@testset "direct construction guards" begin
    @test_throws ArgumentError KernelBackend(parallel = :threads)
    @test_throws ArgumentError ODEProblem(eqs_s, u0, (0.0, 1.0))
    @test_throws ArgumentError ODEFunction(eqs_s)
    @test_throws DimensionMismatch ODEProblem(
        eqs_s,
        zeros(ComplexF64, length(u0) + 1),
        (0.0, 1.0),
        ps;
        backend = KernelBackend(),
    )

    # The standalone ODEFunction constructor is supported even though parameter state is kept
    # separately by ODEProblem; constructing it exercises the public direct-function surface.
    f = ODEFunction(eqs_s, ps; backend = KernelBackend())
    @test f isa SciMLBase.ODEFunction

    @test !QC._resolve_kernel_parallel(:auto, 1)
    @test QC._resolve_kernel_parallel(false, 10_000) === false
    @test QC._resolve_kernel_parallel(true, 1) === true
end
