using QuantumCumulants
using Symbolics: Symbolics, @variables, simplify
using SymbolicUtils
using Test

# v1 surface for parameter handling. Master's test_parameters.jl asserts
# operator-level QAdd equality on `eqs.operator_equations[i].rhs`; in v1
# that field is stored as a `SymbolicUtils.BasicSymbolic{SymReal}` (Symbolics
# unwraps the QAdd through `~`), so the assertions compare via `average(...)`
# on the symbolic side and use `_is_symbolic_zero` to bridge symbolic-vs-Bool
# zero detection (`_is_zero(simplify(x))` on a BasicSymbolic returns a symbolic
# equation, not a Bool).
#
# Ground-state projectors `σ_gg` stay atomic in SQA's canonical form (see
# SQA devdocs "Simplification vs. normal ordering"). Comparing against
# master's expanded forms requires `expand_completeness` on the expected
# side.

SQA = QuantumCumulants.SecondQuantizedAlgebra

# Symbolic / numeric zero check that handles both QAdd and BasicSymbolic.
_is_zero(x::SQA.QAdd) = iszero(x)
_is_zero(x::Number) = iszero(x)
function _is_zero(x::SymbolicUtils.BasicSymbolic)
    s = simplify(x; expand = true)
    s isa Number && return iszero(s)
    SymbolicUtils.isconst(s) && return iszero(s.val)
    return isequal(s, 0)
end
_is_zero(x) = isequal(x, 0)

@testset "parameters: JC commutators with @variables" begin
    hf = FockSpace(:cavity)
    ha = NLevelSpace(:atom, (:g, :e))
    h = hf ⊗ ha
    a = Destroy(h, :a)

    σ = Transition(h, :σ, :g, :e)
    σee = Transition(h, :σ, :e, :e)

    @variables ωc::Real ωa::Real g::Real
    H = ωc * a' * a + ωa * σ' * σ + g * (a' * σ + σ' * a)

    da = commutator(1im * H, a)
    @test _is_zero(simplify(da - ((-1im * g) * σ + (-1im * ωc) * a)))
    # `ds` keeps σ_gg atomic (SQA canonical form); master's expected form
    # had σ_gg expanded via completeness. Apply `expand_completeness` to ds
    # before comparing.
    ds_raw = commutator(1im * H, σ)
    ds = SQA.expand_completeness(ds_raw)
    @test _is_zero(simplify(ds - ((-1im * g) * a + (-1im * ωa) * σ + (2im * g) * a * σee)))
    dn = commutator(1im * H, a' * a)
    @test _is_zero(simplify(dn - ((-1im * g) * a' * σ + (1im * g) * a * σ')))
end

@testset "parameters: meanfield averaged RHS matches commutator average" begin
    hf = FockSpace(:cavity)
    ha = NLevelSpace(:atom, (:g, :e))
    h = hf ⊗ ha
    a = Destroy(h, :a)
    σ = Transition(h, :σ, :g, :e)

    @variables ωc::Real ωa::Real g::Real
    H = ωc * a' * a + ωa * σ' * σ + g * (a' * σ + σ' * a)

    he = meanfield([a, σ, a' * a], H)
    da = commutator(1im * H, a)
    ds = SQA.expand_completeness(commutator(1im * H, σ))
    dn = commutator(1im * H, a' * a)

    @test _is_zero(simplify(he.equations[1].rhs - average(da); expand = true))
    @test _is_zero(simplify(he.equations[2].rhs - average(ds); expand = true))
    @test _is_zero(simplify(he.equations[3].rhs - average(dn); expand = true))
end

@testset "parameters: dissipative meanfield" begin
    hf = FockSpace(:cavity)
    ha = NLevelSpace(:atom, (:g, :e))
    h = hf ⊗ ha
    a = Destroy(h, :a)
    σ = Transition(h, :σ, :g, :e)
    σee = Transition(h, :σ, :e, :e)

    @variables ωc::Real ωa::Real g::Real κ::Real γ::Real
    H = ωc * a' * a + ωa * σ' * σ + g * (a' * σ + σ' * a)
    J = [a, σ]
    he_diss = meanfield([a, σ, σ' * σ], H, J; rates = [κ, γ])

    expected_1 = average((-1im * ωc - 0.5κ) * a + (-1im * g) * σ)
    expected_2 = average((-1im * g) * a + (-1im * ωa - 0.5γ) * σ + (2im * g) * a * σee)
    expected_3 = average((-γ) * σee + (1im * g) * a' * σ + (-1im * g) * a * σ')

    @test _is_zero(simplify(he_diss.equations[1].rhs - expected_1; expand = true))
    @test _is_zero(simplify(he_diss.equations[2].rhs - expected_2; expand = true))
    @test _is_zero(simplify(he_diss.equations[3].rhs - expected_3; expand = true))
end

@testset "parameters: laser-mode equations with pumping" begin
    hf = FockSpace(:cavity)
    ha = NLevelSpace(:atom, (:g, :e))
    h = hf ⊗ ha
    a = Destroy(h, :a)
    σ = Transition(h, :σ, :g, :e)
    σee = Transition(h, :σ, :e, :e)

    @variables ωc::Real ωa::Real g::Real κ::Real γ::Real ν::Real
    H = ωc * a' * a + ωa * σ' * σ + g * (a' * σ + σ' * a)
    J = [a, σ, σ']
    he_laser = meanfield([a' * a, σ' * σ, a' * σ], H, J; rates = [κ, γ, ν])

    expected_n = average((-κ) * a' * a + (-1im * g) * a' * σ + (1im * g) * a * σ')
    expected_pop = average(ν + (-ν - γ) * σee + (1im * g) * a' * σ + (-1im * g) * a * σ')

    @test _is_zero(simplify(he_laser.equations[1].rhs - expected_n; expand = true))
    @test _is_zero(simplify(he_laser.equations[2].rhs - expected_pop; expand = true))
end

@testset "parameters: detuned two-level commutator sign" begin
    @variables Δ::Real
    h2 = NLevelSpace(:a, 2)
    s(i, j) = Transition(h2, :s, i, j)
    H = Δ * s(2, 2) - Δ * s(1, 1)
    eqs = meanfield(s(1, 2), H)
    # Use `Symbolics.IM` (symbolic imaginary unit) so the literal matches the
    # form the meanfield codepath emits; Julia's `im` becomes a `complex(0, …)`
    # literal that SymbolicUtils does not unify with the factored form.
    diff = eqs.equations[1].rhs + 2 * Symbolics.IM * Δ * average(s(1, 2))
    @test _is_zero(simplify(diff; expand = true))
end
