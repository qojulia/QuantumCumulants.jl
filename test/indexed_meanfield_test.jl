using QuantumCumulants
using SecondQuantizedAlgebra: SecondQuantizedAlgebra
using Symbolics: Symbolics, @variables
using Test

@testset "indexed meanfield: collective atomic emission" begin
    hc = FockSpace(:cavity); ha = NLevelSpace(:atom, 2); h = hc ⊗ ha
    @qnumbers a::Destroy(h, 1)
    σ(x, y, k) = IndexedOperator(Transition(h, :σ, x, y, 2), k)
    @variables N Δ g κ γ
    i = Index(h, :i, N, ha)
    H = Δ * a' * a + Σ(g * (a * σ(2, 1, i) + a' * σ(1, 2, i)), i)
    eqs = meanfield([a], H, [a, σ(1, 2, i)]; rates = [κ, γ])
    @test eqs isa MeanFieldEquations
    @test length(eqs.equations) == 1
end

@testset "indexed average appears in RHS" begin
    hc = FockSpace(:cavity); ha = NLevelSpace(:atom, 2); h = hc ⊗ ha
    @qnumbers a::Destroy(h, 1)
    σ(x, y, k) = IndexedOperator(Transition(h, :σ, x, y, 2), k)
    @variables N Δ g κ γ
    i = Index(h, :i, N, ha)
    H = Δ * a' * a + Σ(g * (a * σ(2, 1, i) + a' * σ(1, 2, i)), i)
    eqs = meanfield([a], H, [a, σ(1, 2, i)]; rates = [κ, γ])
    @test occursin("σ", repr(eqs.equations[1].rhs))
end

@testset "indexed meanfield: two-atom collective JC drift" begin
    @variables Δc::Real η::Real Δa::Real κ::Real Γ::Real
    N = 2
    hc = FockSpace(:cavity)
    ha = NLevelSpace(Symbol(:atom), 2)
    h = hc ⊗ ha

    i_ind = Index(h, :i, N, ha)
    k_ind = Index(h, :k, N, ha)

    g(k) = IndexedVariable(:g, k)

    @qnumbers a::Destroy(h)
    σ(i, j, k) = IndexedOperator(Transition(h, :σ, i, j), k)

    Hc = Δc * a' * a + η * (a' + a)
    Ha = Δa * Σ(σ(2, 2, i_ind), i_ind)
    Hi = Σ(g(i_ind) * (a' * σ(1, 2, i_ind) + a * σ(2, 1, i_ind)), i_ind)
    H = Hc + Ha + Hi

    J = [a, σ(1, 2, i_ind)]
    rates = [κ, Γ]

    ops = [a, σ(2, 2, k_ind), σ(1, 2, k_ind)]
    eqs = meanfield(ops, H, J; rates = rates, order = 2)
    @test length(eqs.equations) == 3
    @test eqs isa MeanFieldEquations
    ceqs = complete(eqs)
    @test ceqs isa MeanFieldEquations
    @test length(ceqs.equations) >= 3

    # `complete` is idempotent on a closed system: a second pass should
    # not discover any further missing states.
    ceqs_again = complete(ceqs)
    @test length(ceqs_again.equations) == length(ceqs.equations)

    # Every index appearing in the completed system must trace back to the
    # user's declared vocabulary. The completion path canonicalizes leaf
    # averages against the user-provided indices and never invents new
    # ones (names like `i1`, `_c1`, etc. must not appear).
    declared = Set([i_ind.name, k_ind.name])
    for op in ceqs.operators
        for idx in SecondQuantizedAlgebra.get_indices(op)
            @test idx.name in declared
        end
    end
end

@testset "indexed meanfield: alpha-equivalent leaves are one state" begin
    # When `H` is summed over `i_ind` and the user's `ops` use `k_ind`,
    # the leaf `⟨σ_{i_ind}^{12}⟩` (which appears in RHSes from H's sum)
    # and the user-declared state `⟨σ_{k_ind}^{12}⟩` denote the same
    # physical observable. `find_missing` must recognize the
    # alpha-equivalence and not introduce a duplicate state.
    @variables Δc::Real η::Real Δa::Real κ::Real Γ::Real
    N = 2
    hc = FockSpace(:cavity); ha = NLevelSpace(:atom, 2); h = hc ⊗ ha
    i_ind = Index(h, :i, N, ha); k_ind = Index(h, :k, N, ha)
    @qnumbers a::Destroy(h)
    σ(i, j, k) = IndexedOperator(Transition(h, :σ, i, j), k)
    g(k) = IndexedVariable(:g, k)

    H = Δc * a' * a + η * (a' + a) + Δa * Σ(σ(2, 2, i_ind), i_ind) +
        Σ(g(i_ind) * (a' * σ(1, 2, i_ind) + a * σ(2, 1, i_ind)), i_ind)
    ops = [a, σ(2, 2, k_ind), σ(1, 2, k_ind)]
    J = [a, σ(1, 2, i_ind)]
    eqs = meanfield(ops, H, J; rates = [κ, Γ], order = 2)

    # `find_missing` should not report `⟨σ_i^{12}⟩` (alpha-equivalent
    # to the already-declared `⟨σ_k^{12}⟩`) or `⟨σ_i^{22}⟩` (alpha-
    # equivalent to `⟨σ_k^{22}⟩`) as missing.
    missing_after_first = QuantumCumulants.find_missing(eqs)
    for m in missing_after_first
        s = repr(m)
        # any missing single-σ state must already be one of the
        # declared k-indexed ones; an i-indexed σ-single would be an
        # alpha-duplicate that we explicitly want to deduplicate.
        @test !occursin(r"⟨σ_i₁₂⟩", s)
        @test !occursin(r"⟨σ_i₂₂⟩", s)
        @test !occursin(r"⟨σ_i₂₁⟩", s)
    end
end
