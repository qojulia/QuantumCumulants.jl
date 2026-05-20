using QuantumCumulants
using SymbolicUtils: SymbolicUtils
using Test

@testset "find_missing on JC" begin
    hc = FockSpace(:cavity); ha = NLevelSpace(:atom, 2)
    h = hc ⊗ ha
    @qnumbers a::Destroy(h, 1)
    σ(i, j) = Transition(h, :σ, i, j, 2)
    @variables Δ g κ γ
    H = Δ * a' * a + g * (a * σ(2, 1) + a' * σ(1, 2))
    eqs = meanfield([a, σ(2, 2)], H, [a, σ(1, 2)]; rates = [κ, γ], order = 2)
    missing_states = find_missing(eqs)
    @test !isempty(missing_states)
    @test all(
        SymbolicUtils.iscall(m) && QuantumCumulants.SecondQuantizedAlgebra.is_average(m)
            for m in missing_states
    )
end

@testset "complete! closes JC at order 2" begin
    hc = FockSpace(:cavity); ha = NLevelSpace(:atom, 2)
    h = hc ⊗ ha
    @qnumbers a::Destroy(h, 1)
    σ(i, j) = Transition(h, :σ, i, j, 2)
    @variables Δ g κ γ
    H = Δ * a' * a + g * (a * σ(2, 1) + a' * σ(1, 2))
    eqs = meanfield([a, σ(2, 2)], H, [a, σ(1, 2)]; rates = [κ, γ], order = 2)
    complete!(eqs)
    @test isempty(find_missing(eqs))
    @test length(eqs.equations) >= 2
end

@testset "complete is non-mutating" begin
    hc = FockSpace(:cavity)
    @qnumbers a::Destroy(hc)
    @variables ω κ U
    H = ω * a' * a + U * a' * a' * a * a
    eqs = meanfield([a], H, [a]; rates = [κ], order = 2)
    n_before = length(eqs.equations)
    eqs2 = complete(eqs)
    @test length(eqs.equations) == n_before
    @test isempty(find_missing(eqs2))
end

@testset "find_missing: closure size is invariant of user index naming" begin
    # The indexed JC laser closure under `complete!` has a canonical-class
    # cardinality that depends only on the algebra and the cumulant order,
    # not on the symbol the user picked for the atom index. Two systems
    # that differ only in `k` vs `l` must close to the same equation count.
    # Previously the sequential `_canonical_key` rename let det/stoch
    # passes pick different concrete representatives mid-iteration, which
    # is the same family of failure mode.
    hc = FockSpace(:cavity); ha = NLevelSpace(:atom, 2); h = hc ⊗ ha
    @qnumbers a::Destroy(h, 1)
    σ(α, β, idx) = IndexedOperator(Transition(h, :σ, α, β, 2), idx)
    @variables N::Real Δ::Real g::Real κ::Real Γ::Real R::Real ν::Real

    function closure_size(idx_name::Symbol)
        idx = Index(h, idx_name, N, ha)
        H = -Δ * a' * a + g * (Σ(a' * σ(1, 2, idx), idx) + Σ(a * σ(2, 1, idx), idx))
        J = [a, σ(1, 2, idx), σ(2, 1, idx), σ(2, 2, idx)]
        eqs = complete(
            meanfield(
                [a' * a, σ(2, 2, idx)], H, J;
                rates = [κ, Γ, R, ν], order = 1
            )
        )
        return length(eqs.equations)
    end

    @test closure_size(:k) == closure_size(:l)
end

@testset "find_missing tracks conjugate pairs deterministically" begin
    # Driven cavity with a Fock-coherent observable: the drift of `a` puts
    # `⟨a⟩` on the LHS but the RHS never references `⟨a'⟩`, so the closure
    # has a single state. Switching to `a'` as the observable should yield
    # the conjugate scenario with a single state whose canonical key is
    # the conjugate canonical key of the first. Asserts that
    # `_collect_missing!`'s conjugate handling is symmetric: starting from
    # either rep, the canonical-class identity of the state set is the
    # same.
    hc = FockSpace(:cavity)
    @qnumbers a::Destroy(hc)
    @variables ω κ η
    H = ω * a' * a + 1im * η * (a' - a)

    eqs_a = complete(meanfield(a, H, [a]; rates = [κ]))
    eqs_ad = complete(meanfield(a', H, [a]; rates = [κ]))

    canon_a = QuantumCumulants._build_canonical_indices(eqs_a)
    canon_ad = QuantumCumulants._build_canonical_indices(eqs_ad)
    key_a = QuantumCumulants._canonical_key(only(eqs_a.states), canon_a)
    key_ad = QuantumCumulants._canonical_key(only(eqs_ad.states), canon_ad)
    # The state of one is the conjugate of the other's; the canonical key
    # of one matches the conjugate canonical key of the other.
    conj_key_a = QuantumCumulants._canonical_key(
        QuantumCumulants._avg_conj(only(eqs_a.states)), canon_a,
    )
    @test key_ad == conj_key_a
end
