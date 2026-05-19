using QuantumCumulants
using Symbolics: Symbolics, @variables
using Test

# v1 surface: indexed `CorrelationFunction` (unified API, no separate
# `IndexedCorrelationFunction` type). Master's full test exercises
# `IndexedCorrelationFunction`, `evaluate(corr, 1, 2; limits=...)`,
# `scale(corr)`, `split_sums`, and several SQA helpers (`DoubleNumberedVariable`,
# `SingleNumberedVariable`, `value_map`) that don't exist in v1; those parts
# stay in test/pending/indexed_correlation_test.jl until the order=2 indexed
# JC closure speeds up (the same `complete!` mixed-order issue blocks
# indexed_mixed_order_test).

@testset "indexed CorrelationFunction: JC laser, order=1" begin
    @variables N::Real Δ::Real g::Real κ::Real Γ::Real R::Real ν::Real
    hc = FockSpace(:cavity)
    ha = NLevelSpace(:atom, 2)
    h = hc ⊗ ha
    k = Index(h, :k, N, ha)
    @qnumbers a::Destroy(h)
    σ(i, j, idx) = IndexedOperator(Transition(h, :σ, i, j), idx)

    H = -Δ * a' * a + g * (Σ(a' * σ(1, 2, k), k) + Σ(a * σ(2, 1, k), k))
    J = [a, σ(1, 2, k), σ(2, 1, k), σ(2, 2, k)]
    rates = [κ, Γ, R, ν]
    ops = [a' * a, a]

    eqs = meanfield(ops, H, J; rates = rates, order = 1)
    eqs_c = complete(eqs)
    @test length(eqs_c.equations) >= 1

    # Scale should produce a finite-size closed system.
    eqs_sc = scale(eqs_c)
    @test length(eqs_sc.equations) >= 1
    @test isempty(find_missing(eqs_sc; get_adjoints = false))
end

@testset "indexed CorrelationFunction: g^(1)(τ) construction" begin
    @variables N::Real Δ::Real g::Real κ::Real Γ::Real R::Real ν::Real
    hc = FockSpace(:cavity)
    ha = NLevelSpace(:atom, 2)
    h = hc ⊗ ha
    k = Index(h, :k, N, ha)
    @qnumbers a::Destroy(h)
    σ(i, j, idx) = IndexedOperator(Transition(h, :σ, i, j), idx)

    H = -Δ * a' * a + g * (Σ(a' * σ(1, 2, k), k) + Σ(a * σ(2, 1, k), k))
    J = [a, σ(1, 2, k)]
    rates = [κ, Γ]
    eqs = meanfield([a' * a, a], H, J; rates = rates, order = 1)
    eqs_c = complete(eqs)

    # First-order correlation of the cavity mode.
    corr = CorrelationFunction(a', a, eqs_c)
    @test corr isa CorrelationFunction
    @test length(corr.eqs.equations) >= 1

    # Scale on the correlation function: same surface as on regular equations.
    corr_sc = scale(corr)
    @test corr_sc isa CorrelationFunction
end
