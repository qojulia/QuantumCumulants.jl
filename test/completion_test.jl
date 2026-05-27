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

@testset "_build_canonical_indices filters H/J-bound indices" begin
    # Regression: a free LHS atom index must not appear in the canonical
    # pool if the same name is also bound by an H sum or a collective
    # jump. Pre-fix, `complete!` would happily map a LHS `i` onto H's
    # bound `i`, silently dropping the cross-atom commutator on the
    # derived equation.
    @variables Δ::Real g::Real κ::Real Γ::Real N::Real
    hc = FockSpace(:cavity); ha = NLevelSpace(:atom, 2); h = hc ⊗ ha
    @qnumbers a::Destroy(h, 1)
    σ(α, β, k) = IndexedOperator(Transition(h, :σ, α, β, 2), k)
    i_idx = Index(h, :i, N, ha)
    j_idx = Index(h, :j, N, ha)

    # `i` is bound by H's sum. `j` is the only free atom index.
    H = Δ * Σ(σ(2, 2, i_idx), i_idx) +
        g * Σ(a' * σ(1, 2, i_idx) + a * σ(2, 1, i_idx), i_idx)
    J = [a, σ(1, 2, i_idx)]
    eqs = meanfield([a' * a, σ(2, 2, j_idx)], H, J; rates = [κ, Γ], order = 2)
    canon = QuantumCumulants._build_canonical_indices(eqs)

    bound = QuantumCumulants._bound_indices(eqs)
    @test i_idx in bound
    @test !(j_idx in bound)

    atom_space = i_idx.space_index
    @test haskey(canon, atom_space)
    @test j_idx in canon[atom_space]
    @test !(i_idx in canon[atom_space])
end

@testset "_build_canonical_indices fallback when every user index is bound" begin
    # When the user declares only indices that all get bound (e.g. on a
    # cavity subspace where every index is consumed by an H sum), the
    # filtered pool is empty. The fallback mints `original[1](2)` so
    # downstream `_canonical_key` lookups produce a deterministic
    # suffix like `:i_2` instead of colliding with the bound name.
    # Setup mirrors the filter-cavity example: composite space, single
    # index symbol `i` that is bound on the filter subspace by every
    # source it appears in.
    @variables κ::Real ω::Real N::Real M::Real
    hc = FockSpace(:cavity); hf = FockSpace(:filter)
    ha = NLevelSpace(:atom, 2)
    h = hc ⊗ hf ⊗ ha
    @qnumbers a::Destroy(h, 1)
    b(k) = IndexedOperator(Destroy(h, :b, 2), k)
    i_idx = Index(h, :i, M, hf)
    H = ω * Σ(b(i_idx)' * b(i_idx), i_idx)
    J = [a, b(i_idx)]
    eqs = meanfield([a' * a], H, J; rates = [κ, κ], order = 1)
    canon = QuantumCumulants._build_canonical_indices(eqs)
    filter_space = i_idx.space_index
    @test haskey(canon, filter_space)
    fallback = only(canon[filter_space])
    @test fallback.name == :i_2
end

@testset "_canonical_key strips irrelevant NE pairs" begin
    # Two leaves that differ only in an NE pair referencing an index
    # absent from the operator atoms should canonicalise to the same
    # key. Pre-fix, the inherited NE pair survived in the key and
    # blocked dedup of equivalent states pulled in from a parent
    # derivation.
    @variables N::Real
    hc = FockSpace(:cavity); ha = NLevelSpace(:atom, 2); h = hc ⊗ ha
    @qnumbers a::Destroy(h, 1)
    σ(α, β, k) = IndexedOperator(Transition(h, :σ, α, β, 2), k)
    j_idx = Index(h, :j, N, ha)
    k_idx = Index(h, :k, N, ha)
    ℓ_idx = Index(h, :ℓ, N, ha)

    SQA = QuantumCumulants.SecondQuantizedAlgebra
    plain_qadd = σ(1, 2, j_idx) * 1  # QAdd form
    plain = average(plain_qadd)
    # Inject an NE-pair (k, ℓ) onto the QAdd. k and ℓ are not in the
    # operator atoms, so the canonical key should treat the leaf as
    # equivalent to the bare one.
    qadd_with_ne = SQA.assume_distinct_index(plain_qadd, [(k_idx, ℓ_idx)])
    with_ne = average(qadd_with_ne)

    eqs = meanfield([σ(1, 2, j_idx)], 0 * a; order = 1)
    canon = QuantumCumulants._build_canonical_indices(eqs)
    key_plain = QuantumCumulants._canonical_key(plain, canon)
    key_with_ne = QuantumCumulants._canonical_key(with_ne, canon)
    @test key_plain == key_with_ne
end

@testset "complete!: free LHS atom index does not collide with H-bound i" begin
    # End-to-end check that `_derive_for`'s alpha-rename + `_apply_undo`
    # round trip preserves the cross-atom commutator that was being
    # silently dropped pre-fix. We use the indexed-laser closure: with
    # H summed over `i` and the user-declared free index `i` (clashing
    # name), the σ(2,2,i) equation's RHS must include a
    # `(N-1) i g ⟨σ_{i_1}^{12} σ_{i_2}^{21}⟩` style cross-atom term
    # after `complete!`. We assert by checking that the dimension of
    # the closure grows past the single-atom case.
    @variables Δ::Real g::Real κ::Real Γ::Real N::Real
    hc = FockSpace(:cavity); ha = NLevelSpace(:atom, 2); h = hc ⊗ ha
    @qnumbers a::Destroy(h, 1)
    σ(α, β, k) = IndexedOperator(Transition(h, :σ, α, β, 2), k)
    i_idx = Index(h, :i, N, ha)  # deliberately reuse `i` for LHS

    H = Δ * Σ(σ(2, 2, i_idx), i_idx) +
        g * Σ(a' * σ(1, 2, i_idx) + a * σ(2, 1, i_idx), i_idx)
    J = [a, σ(1, 2, i_idx)]
    # Note: free LHS uses the SAME index symbol as H's bound sum.
    eqs = meanfield([a' * a, σ(2, 2, i_idx)], H, J; rates = [κ, Γ], order = 2)
    eqs_c = complete(eqs)

    # Pre-fix: closure collapsed to 4-5 equations because cross-atom
    # commutators were dropped. Post-fix: the cumulant order-2 indexed
    # JC laser closes at >= 7 equations (atom-atom cross states present).
    @test length(eqs_c.equations) >= 7
end
