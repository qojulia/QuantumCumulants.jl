using QuantumCumulants
using Symbolics: Symbolics, @variables
using SymbolicUtils
using Test

# v1 surface: comparison between indexed noise meanfield + scale and the
# explicit per-atom expansion path. Master compares equation-by-equation
# at order=4; here we assert that both pipelines (indexed + scale, vs
# `indexed_complete` of the stochastic form) reach a closed system on the
# same input. Subset of master's `test_measurement_backaction_indices_comparison`
# that doesn't need `pulse(t)`/`@syms` or evaluate-at-fixed-indices.

@testset "measurement_backaction_indices_comparison: pipeline closure" begin
    @variables N::Real ωa::Real γ::Real η::Real χ::Real ωc::Real κ::Real g::Real ξ::Real

    hc = FockSpace(:resonator)
    ha = NLevelSpace(:atom, 2)
    h = hc ⊗ ha

    j = Index(h, :j, N, ha)
    k = Index(h, :k, N, ha)

    @qnumbers a::Destroy(h, 1)
    σ(α, β, k) = IndexedOperator(Transition(h, :σ, α, β, 2), k)

    H = ωc * a' * a + ωa * Σ(σ(2, 2, j), j) +
        g * a' * Σ(σ(1, 2, j), j) + g * a * Σ(σ(2, 1, j), j)
    J = [a, σ(1, 2, j), σ(2, 1, j), σ(2, 2, j)]
    rates = [κ, γ, η, 2 * χ]
    efficiencies = [ξ, 0, 0, 0]
    ops = [a, a' * a, σ(2, 2, k), σ(1, 2, k)]

    eqs = meanfield(ops, H, J; rates = rates, efficiencies = efficiencies, order = 2)
    @test eqs isa NoiseMeanFieldEquations
    eqs_c = complete(eqs)
    @test length(eqs_c.equations) >= length(ops)
    @test isempty(find_missing(eqs_c; get_adjoints = false))

    # Scale the closed indexed noise system.
    scaled_eqs = scale(eqs_c)
    @test scaled_eqs isa NoiseMeanFieldEquations
    @test length(scaled_eqs.equations) >= 1
    # Noise channel survives scaling. Counts may differ from deterministic
    # part by one if a state collapses under scale while its noise term
    # does not; both lists should be non-empty.
    @test length(scaled_eqs.noise_equations) >= 1
end

@testset "measurement_backaction_indices_comparison: deterministic vs stochastic LHS match" begin
    @variables N::Real ωa::Real γ::Real η::Real χ::Real ωc::Real κ::Real g::Real ξ::Real

    hc = FockSpace(:resonator)
    ha = NLevelSpace(:atom, 2)
    h = hc ⊗ ha
    j = Index(h, :j, N, ha)
    k = Index(h, :k, N, ha)

    @qnumbers a::Destroy(h, 1)
    σ(α, β, k) = IndexedOperator(Transition(h, :σ, α, β, 2), k)

    H = ωc * a' * a + ωa * Σ(σ(2, 2, j), j) +
        g * a' * Σ(σ(1, 2, j), j) + g * a * Σ(σ(2, 1, j), j)
    J = [a, σ(1, 2, j), σ(2, 1, j), σ(2, 2, j)]
    rates = [κ, γ, η, 2 * χ]
    efficiencies = [ξ, 0, 0, 0]
    ops = [a, a' * a, σ(2, 2, k), σ(1, 2, k)]

    stoch_eqs = meanfield(
        ops, H, J;
        rates = rates, efficiencies = efficiencies, order = 2
    )
    det_eqs = meanfield(ops, H, J; rates = rates, order = 2)
    stoch_c = complete(stoch_eqs)
    det_c = complete(det_eqs)

    # det and stoch now select ONE structural NE policy (population iff the
    # channel set contains a dephasing jump `σ^{αα}`), injected per-derived-op
    # (independent of completion-iteration order), applied identically to both
    # pipelines. So they reach the SAME closure with identical drifts.
    det_lhs = Set(e.lhs for e in det_c.equations)
    stoch_lhs = Set(e.lhs for e in stoch_c.equations)
    @test det_lhs == stoch_lhs
    @test !isempty(det_lhs)

    # Every shared deterministic drift agrees symbolically. The noise
    # meanfield's `equations` are the deterministic drift terms; only
    # `noise_equations` carry the stochastic part. So for each shared LHS the
    # stoch and det rhs must differ by zero.
    det_by_lhs = Dict(e.lhs => e.rhs for e in det_c.equations)
    n_checked = 0
    for eq in stoch_c.equations
        haskey(det_by_lhs, eq.lhs) || continue
        @test _is_zero(eq.rhs - det_by_lhs[eq.lhs])
        n_checked += 1
    end
    @test n_checked == length(stoch_c.equations)
end
