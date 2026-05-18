# PENDING: port of master test/test_indexed_correlation.jl
#
# Status: master's 418-line test exercises a large surface of removed /
# unimplemented v1 APIs —
#   * `IndexedCorrelationFunction` (collapsed into `CorrelationFunction`
#     for indexed inputs; the indexed semantics behind it still need
#     porting — the v1 unified `CorrelationFunction` does not yet handle
#     `IndexedOperator` op1/op2 the same way).
#   * `evaluate(eqs_c; limits = (N => k))` to unroll Σ-sums.
#   * `evaluate(corr, 1, 2; limits = ...)` to evaluate at fixed indices.
#   * `scale(corr)` per-correlation scaling.
#   * `split_sums(avg_sum, idx, k)` to partition Σ-sums into chunks.
#   * `value_map`, `DoubleNumberedVariable`, `SingleNumberedVariable`,
#     `scale_term` — SQA-internal helpers no longer in v0.5's surface.
#
# Faithful porting of every assertion is not useful while these APIs are
# missing. The skeleton below covers the core scenarios that *should* keep
# working once the features land: closure of the indexed JC laser at order
# 2, ODE roundtrip, CorrelationFunction construction with phase filter.

#=
using QuantumCumulants
using Symbolics: Symbolics, @variables
using SymbolicUtils
using ModelingToolkitBase: @named, mtkcompile, ODEProblem, unknowns
using OrdinaryDiffEq: Tsit5, solve
using Test

@testset "indexed_correlation: laser closure + correlation construction" begin
    order = 2
    @variables N::Real Δ::Real g::Real κ::Real Γ::Real R::Real ν::Real M::Real

    hc = FockSpace(:cavity)
    ha = NLevelSpace(:atom, 2)
    h = hc ⊗ ha

    k = Index(h, :k, N, ha)
    l = Index(h, :l, N, ha)
    m = Index(h, :m, N, ha)

    @qnumbers a::Destroy(h)
    σ(i, j, k) = IndexedOperator(Transition(h, :σ, i, j), k)

    H = -Δ * a' * a + g * (Σ(a' * σ(1, 2, k), k) + Σ(a * σ(2, 1, k), k))
    J = [a, σ(1, 2, l), σ(2, 1, l), σ(2, 2, l)]
    rates = [κ, Γ, R, ν]
    ops = [a' * a, σ(2, 2, m)]
    eqs = meanfield(ops, H, J; rates = rates, order = order)
    @test length(eqs.equations) == 2

    SQA = QuantumCumulants.SecondQuantizedAlgebra
    ϕ(x) = 0
    ϕ(::Destroy) = -1
    ϕ(::Create) = 1
    function ϕ(t::Transition)
        t.i == t.j && return 0
        return t.i == 2 ? 1 : -1
    end
    ϕ(q::SQA.QAdd) = sum(ϕ(arg) for (arg, _) in q.arguments)
    ϕ(q::SQA.QTerm) = sum(ϕ(op) for op in q.ops)
    ϕ(op::IndexedOperator) = ϕ(op.op)
    function ϕ(avg)
        avg isa SymbolicUtils.BasicSymbolic && SQA.is_average(avg) || return 0
        return ϕ(SQA.undo_average(avg))
    end
    phase_invariant(x) = iszero(ϕ(x))

    eqs_c = complete(eqs; filter_func = phase_invariant)
    @test isempty(find_missing(eqs_c; filter_func = phase_invariant))

    eqs_sc1 = scale(eqs_c)
    @test length(eqs_sc1.equations) >= 1

    op1 = σ(1, 2, m)
    op2 = σ(2, 1, Index(h, :n, N, ha))
    corr = CorrelationFunction(op1, op2, eqs_c; filter_func = phase_invariant)
    @test corr isa CorrelationFunction
    corr_sc = scale(corr)
    @test corr_sc isa CorrelationFunction
end
=#
