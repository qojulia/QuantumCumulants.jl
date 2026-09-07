using QuantumCumulants
using Symbolics: Symbolics, @variables
using Test

const QC = QuantumCumulants

@testset "MomentIR unresolved and time-dependent errors" begin
    h = PauliSpace(:errors)
    sx = Pauli(h, :σ, 1)
    sy = Pauli(h, :σ, 2)
    sz = Pauli(h, :σ, 3)
    sm = (sx - 1im * sy) / 2
    @variables J hx γ
    H = -J * sz - hx * sx

    open = meanfield([sz], H, [sm]; rates = [γ], order = 2)
    @test_throws QC.UnresolvedMomentError QC._lower_moment_ir(open)

    closed = complete(open; get_adjoints = true)
    timed = modify_equations(closed, (op, drift) -> cos(closed.iv) * drift)
    @test_throws QC.TimeDependentCoefficientError QC._lower_moment_ir(timed)

    nonpoly = modify_equations(closed, (op, drift) -> drift + exp(average(op)))
    @test_throws QC.NonPolynomialDriftError QC._lower_moment_ir(nonpoly)
end
