using QuantumCumulants
using Test

import QuantumCumulants.SecondQuantizedAlgebra as SQA

# A redundant ground-state projector σ^{gg} must not survive into a post-cumulant
# closure; the basis change σ^{gg} = 1 - Σ σ^{kk} should already have folded it away.

_qadd_has_gs(::Any) = false
function _qadd_has_gs(op::SQA.QAdd)
    for (term, _) in op.arguments, o in term.ops
        SQA.is_transition(o) || continue
        o.l1 == o.l2 == o.g && return true
    end
    return false
end

function _state_has_gs(s)
    op = SQA.undo_average(s)
    return _qadd_has_gs(op)
end

function leaks_ground_projector(eqs)
    return any(_state_has_gs, eqs.states)
end

@testset "completeness invariant (no σ^gg leak)" begin
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

    @testset "scalar two-level" begin
        @variables Δ g κ γ
        hc = FockSpace(:cavity); ha = NLevelSpace(:atom, 2); h = hc ⊗ ha
        a = Destroy(h, :a); σ(α, β) = Transition(h, :σ, α, β)
        H = Δ * a' * a + g * (a' * σ(1, 2) + a * σ(2, 1))
        eqs_c = complete(meanfield([a' * a, σ(2, 2)], H, [a, σ(1, 2)]; rates = [κ, γ], order = 2))
        @test !leaks_ground_projector(eqs_c)
    end
end
