using QuantumCumulants
using Symbolics: Symbolics, @variables
using SymbolicUtils
using Test

# 6th-order cumulant closure of the damped JC system.

@testset "higher-order: 6th-order closure runs (no spectrum assertion)" begin
    @variables Δ::Real g::Real γ::Real κ::Real ν::Real
    hf = FockSpace(:cavity)
    ha = NLevelSpace(:atom, 2)
    h = hf ⊗ ha
    a = Destroy(h, :a)
    s = Transition(h, :σ, 1, 2)

    H = Δ * a' * a + g * (a' * s + a * s')
    J = [a, s, s']
    rates = [κ, γ, ν]

    SQA = QuantumCumulants.SecondQuantizedAlgebra
    ϕ(::Number) = 0
    ϕ(::SQA.Destroy) = -1
    ϕ(::SQA.Create) = 1
    function ϕ(t::SQA.Transition)
        t.i == t.j && return 0
        return t.i == 2 ? 1 : -1
    end
    ϕ(q::SQA.QAdd) = isempty(q.arguments) ? 0 : sum(ϕ(term) for (term, _) in q.arguments)
    ϕ(t::SQA.QTerm) = isempty(t.ops) ? 0 : sum(ϕ(op) for op in t.ops)
    function ϕ(avg)
        avg isa SymbolicUtils.BasicSymbolic && SQA.is_average(avg) || return 0
        return ϕ(SQA.undo_average(avg))
    end
    phase_invariant(x) = iszero(ϕ(x))

    he6 = complete(
        meanfield(a' * a, H, J; rates = rates);
        order = 6, filter_func = phase_invariant
    )
    @test length(he6.equations) > 0
    @test isempty(find_missing(he6; filter_func = phase_invariant))
end
