using QuantumCumulants
using Symbolics: Symbolics, @variables
using SymbolicUtils
using Test

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
    @test eqs isa NoiseMeanfieldEquations
    eqs_c = complete(eqs)
    @test length(eqs_c.equations) >= length(ops)
    @test isempty(find_missing(eqs_c))

    scaled_eqs = scale(eqs_c)
    @test scaled_eqs isa NoiseMeanfieldEquations
    @test length(scaled_eqs.equations) >= 1
    # The noise channel survives scaling.
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

    # The deterministic and stochastic pipelines reach the same closure with
    # identical LHS sets.
    det_lhs = Set(e.lhs for e in det_c.equations)
    stoch_lhs = Set(e.lhs for e in stoch_c.equations)
    @test det_lhs == stoch_lhs
    @test !isempty(det_lhs)

    # Every shared deterministic drift agrees symbolically. The noise meanfield's
    # `equations` carry only the deterministic drift; the stochastic part lives in
    # `noise_equations`.
    det_by_lhs = Dict(e.lhs => e.rhs for e in det_c.equations)
    n_checked = 0
    for eq in stoch_c.equations
        haskey(det_by_lhs, eq.lhs) || continue
        @test _is_zero(eq.rhs - det_by_lhs[eq.lhs])
        n_checked += 1
    end
    @test n_checked == length(stoch_c.equations)
end
