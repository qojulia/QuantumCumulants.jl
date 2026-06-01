using QuantumCumulants
using Test

import QuantumCumulants.SecondQuantizedAlgebra as SQA

# Defensive invariant for TODO §4: no completeness-redundant ground-state
# projector `σ^{gg}` (a Transition with i == j == ground_state) may survive
# into a post-cumulant closure. `expand_completeness` folds `σ^{gg} = 1 - Σ
# σ^{kk}`; the mean-field/noise drift builders call it, but the fold is
# point-local ("every caller remembers to fold"). This test catches a future
# operator-product site that forgets the fold and leaks a redundant `⟨σ^{gg}⟩`
# leaf — the failure mode behind the det/stoch closure-overcount bugs.
#
# Detection works on GENUINE leaf averages only (`⟨single operator string⟩`).
# A *product* of averages such as `⟨σ¹²⟩·⟨a σ²¹⟩` also reports
# `SQA.is_average`, and `undo_average` would reassemble it into a same-site
# product that collapses to a spurious `σ^{gg}` — a false positive. Guarding on
# `_is_leaf_average` avoids that.

_qadd_has_gs(::Any) = false
function _qadd_has_gs(op::SQA.QAdd)
    for (term, _) in op.arguments, o in term.ops
        o isa SQA.Transition || continue
        o.i == o.j == o.ground_state && return true
    end
    return false
end

function _leaf_has_gs(x)
    x isa SymbolicUtils.BasicSymbolic || return false
    if QuantumCumulants._is_leaf_average(x)
        op = SQA.undo_average(x)
        return op isa SQA.QAdd && _qadd_has_gs(op)
    end
    SymbolicUtils.iscall(x) || return false
    return any(_leaf_has_gs, SymbolicUtils.arguments(x))
end

# Scan every state and every (drift + noise) RHS of an equation set.
function leaks_ground_projector(eqs)
    any(_leaf_has_gs, eqs.states) && return true
    any(eq -> _leaf_has_gs(eq.rhs), eqs.equations) && return true
    if hasproperty(eqs, :noise_equations)
        any(eq -> _leaf_has_gs(eq.rhs), eqs.noise_equations) && return true
    end
    return false
end

@testset "completeness invariant (no σ^gg leak)" begin
    # Indexed two-level ensemble + cavity with the full raising/lowering/
    # dephasing channel set (heterodyne shape). Same-site `σ²¹_k · σ¹²_k`
    # collapses inside cross-atom cumulant blocks are exactly what can
    # surface a redundant `σ^{gg}`. Check both `get_adjoints` settings and
    # the scaled closure.
    @testset "indexed 2-level (raising+lowering+dephasing)" begin
        @variables N Δc Δa g κ γ χ
        hc = FockSpace(:cavity); ha = NLevelSpace(:atom, 2); h = hc ⊗ ha
        k = Index(h, :k, N, ha)
        @qnumbers a::Destroy(h, 1)
        σ(α, β, i) = IndexedOperator(Transition(h, :σ, α, β, 2), i)
        H = Δc * a' * a + Δa * Σ(σ(2, 2, k), k) +
            g * (a' * Σ(σ(1, 2, k), k) + a * Σ(σ(2, 1, k), k))
        J = [a, σ(1, 2, k), σ(2, 1, k), σ(2, 2, k)]
        rates = [κ, γ, γ, 2 * χ]
        ops = [a' * a, σ(2, 2, k), σ(1, 2, k)]
        for ga in (true, false)
            eqs = meanfield(ops, H, J; rates = rates, order = 2)
            eqs_c = complete(eqs; get_adjoints = ga)
            @test !leaks_ground_projector(eqs_c)
            @test !leaks_ground_projector(scale(eqs_c))
        end
    end

    # Indexed Dicke (cavity + N two-level atoms, single decay channel).
    @testset "indexed Dicke" begin
        @variables N Δc Δa g κ Γ
        hc = FockSpace(:cavity); ha = NLevelSpace(:atom, 2); h = hc ⊗ ha
        i = Index(h, :i, N, ha)
        @qnumbers a::Destroy(h, 1)
        σ(α, β, x) = IndexedOperator(Transition(h, :σ, α, β, 2), x)
        H = Δc * a' * a + Δa * Σ(σ(2, 2, i), i) +
            g * (a' * Σ(σ(1, 2, i), i) + a * Σ(σ(2, 1, i), i))
        eqs_c = complete(meanfield([a' * a, σ(2, 2, i)], H, [a, σ(1, 2, i)]; rates = [κ, Γ], order = 2))
        @test !leaks_ground_projector(eqs_c)
        @test !leaks_ground_projector(scale(eqs_c))
    end

    # Multilevel (three-level) indexed system: σ^{gg} folds onto level 1, and
    # the level-2/3 cross terms still must not resurface a ground projector.
    @testset "indexed 3-level" begin
        @variables N Δ Ω g κ Γ
        hc = FockSpace(:cavity); ha = NLevelSpace(:atom, 3); h = hc ⊗ ha
        i = Index(h, :i, N, ha)
        @qnumbers a::Destroy(h, 1)
        σ(α, β, x) = IndexedOperator(Transition(h, :σ, α, β, 2), x)
        H = Δ * a' * a + Ω * Σ(σ(3, 1, i) + σ(1, 3, i), i) +
            g * (a' * Σ(σ(1, 2, i), i) + a * Σ(σ(2, 1, i), i))
        eqs_c = complete(meanfield([a' * a, σ(2, 2, i), σ(3, 3, i)], H, [a, σ(1, 2, i), σ(1, 3, i)]; rates = [κ, Γ, Γ], order = 2))
        @test !leaks_ground_projector(eqs_c)
        @test !leaks_ground_projector(scale(eqs_c))
    end

    # Scalar (non-indexed) two-level + cavity, for coverage of the
    # non-summed product path.
    @testset "scalar two-level" begin
        @variables Δ g κ γ
        hc = FockSpace(:cavity); ha = NLevelSpace(:atom, 2); h = hc ⊗ ha
        a = Destroy(h, :a); σ(α, β) = Transition(h, :σ, α, β)
        H = Δ * a' * a + g * (a' * σ(1, 2) + a * σ(2, 1))
        eqs_c = complete(meanfield([a' * a, σ(2, 2)], H, [a, σ(1, 2)]; rates = [κ, γ], order = 2))
        @test !leaks_ground_projector(eqs_c)
    end
end
