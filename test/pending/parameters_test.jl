# PENDING: port of master test/test_parameters.jl
#
# Status: blocked on two v1 issues —
#   1. `eqs.operator_equations[i].rhs` is stored as a SymbolicUtils SymReal
#      symbolic instead of a QAdd, so the master assertions
#      `iszero(simplify(he.operator_equations[i].rhs - da))` cannot subtract
#      a QAdd `da` from the SymReal-typed stored rhs.
#   2. SQA v0.5 dropped the legacy `QTerm` exported alias path used in
#      `@test p*q*a isa QuantumCumulants.QTerm`. v1 reexports `QTerm` but
#      `p*q*a` is now a `QMul`-style internal form, so the `isa QTerm`
#      assertion needs reformulation.
#
# To enable: unwrap the `#= ... =#` block below, address (1) by exposing a
# QAdd view of operator_equations OR rewriting the assertions to use the
# averaged `eqs.equations[i].rhs` (a SymbolicUtils symbolic, comparable to
# `average(...)`-built targets after `_rewrite_complex_literals`).

#=
using QuantumCumulants
using Symbolics: Symbolics, @variables, simplify
using SymbolicUtils
using Test

@testset "parameters: JC commutators with @variables" begin
    hf = FockSpace(:cavity)
    ha = NLevelSpace(:atom, (:g, :e))
    h = hf ⊗ ha
    a = Destroy(h, :a)

    @variables p::Real q::Real
    @test p * q * a isa QuantumCumulants.SecondQuantizedAlgebra.QTerm  # see (2) above

    σ = Transition(h, :σ, :g, :e)
    σee = Transition(h, :σ, :e, :e)

    @variables ωc::Real ωa::Real g::Real
    H = ωc * a' * a + ωa * σ' * σ + g * (a' * σ + σ' * a)

    da = commutator(1im * H, a)
    @test iszero(simplify(da - ((-1im * g) * σ + (-1im * ωc) * a)))
    ds = commutator(1im * H, σ)
    @test iszero(simplify(ds - ((-1im * g) * a + (-1im * ωa) * σ + (2im * g) * a * σee)))
    dn = commutator(1im * H, a' * a)
    @test iszero(simplify(dn - ((-1im * g) * a' * σ + (1im * g) * a * σ')))

    he = meanfield([a, σ, a' * a], H)
    @test iszero(simplify(he.operator_equations[1].rhs - da))
    @test iszero(simplify(he.operator_equations[2].rhs - ds))
    @test iszero(simplify(he.operator_equations[3].rhs - dn))

    @variables κ::Real γ::Real
    J = [a, σ]
    he_diss = meanfield([a, σ, σ' * σ], H, J; rates = [κ, γ])

    @test iszero(simplify(he_diss.operator_equations[1].rhs -
                          ((-1im * ωc - 0.5κ) * a + (-1im * g) * σ); expand = true))
    @test iszero(simplify(he_diss.operator_equations[2].rhs -
                          ((-1im * g) * a + (-1im * ωa - 0.5γ) * σ +
                           (2im * g) * a * σee); expand = true))
    @test iszero(simplify(he_diss.operator_equations[3].rhs -
                          ((-γ) * σee + (1im * g) * a' * σ + (-1im * g) * a * σ')))

    @variables ν::Real
    J = [a, σ, σ']
    he_laser = meanfield([a' * a, σ' * σ, a' * σ], H, J; rates = [κ, γ, ν])

    @test iszero(simplify(he_laser.operator_equations[1].rhs -
                          ((-κ) * a' * a + (-1im * g) * a' * σ + (1im * g) * a * σ')))
    @test iszero(simplify(he_laser.operator_equations[2].rhs -
                          (ν + (-ν - γ) * σee + (1im * g) * a' * σ +
                           (-1im * g) * a * σ'); expand = true))

    @variables Δ::Real
    h2 = NLevelSpace(:a, 2)
    s(i, j) = Transition(h2, :s, i, j)
    H = Δ * s(2, 2) - Δ * s(1, 1)
    H = simplify(H)
    @test iszero(simplify(meanfield(s(1, 2), H).operator_equations[1].rhs +
                          2im * Δ * s(1, 2)))
end
=#
