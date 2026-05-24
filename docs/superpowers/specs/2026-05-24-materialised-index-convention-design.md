# Materialised-Index Convention (QC v1)

Date: 2026-05-24
Branch: rewrite
Scope: QuantumCumulants.jl (no new SQA changes; depends on already-committed
SQA edits, see "SQA prerequisites")

## Problem

QC's `scale` (cluster-symmetric canonical form) and `evaluate` (per-atom
enumeration) both mint new symbolic `Index` objects when they materialise
a free index into concrete slots. They independently use the same naming
pattern `Symbol(base.name, "_", k)`, so the symbols collide while the
*semantics* differ:

- `scale`'s slot-k rep is "the k-th canonical-form atom of an exchangeable
  cluster"
- `evaluate`'s atom-k is "the concrete atom at position k in a fully
  enumerated sum"

Concretely:

- [src/completion.jl:272-277](../../src/completion.jl#L272-L277)
  `_slot_rep(base, k)` returns `base` for k=1 and
  `Index(Symbol(base.name, "_", k), …)` for k>1.
- [src/evaluate.jl:303-308](../../src/evaluate.jl#L303-L308)
  `_fresh_index(b, k, canon)` always returns
  `Index(Symbol(canon_first.name, "_", k), …)`.

When the two passes are composed (test
`indexed_filter_cavity_test.jl:113-122` exercises both orderings), the
resulting equations have inconsistent literal symbols on LHS vs RHS, so
`find_missing` reports leaves like `⟨σ_j_1₁₂⟩` as missing and MTK
codegen can't bind `AvgFunc` literals to unknowns. Current state of the
`rewrite` branch: 690 pass / 4 fail / 1 error, all failures rooted here.

## SQA already ships the right primitive

`SecondQuantizedAlgebra.jl/src/expressions/index_types.jl:85-89` defines

```julia
function (i::Index)(k::Integer)
    name = Symbol(i.name, "_", k)
    sym_var = SymbolicUtils.Sym{SymbolicUtils.SymReal}(name; type = Int)
    return Index(name, i.range, i.space_index, Num(sym_var))
end
```

with docstring stating it "matches the indices that `evaluate` mints
internally so the resulting operators dedup-equal `evaluate`'s output".
This primitive was added for exactly this collision class. QC just isn't
using it.

## SQA prerequisites (already landed)

Two edits in SQA (`redesign-v2` branch, freshly committed) are
load-bearing for this design:

1. **[src/expressions/qterm.jl](../../../SecondQuantizedAlgebra.jl/src/expressions/qterm.jl)**:
   add `_ne_becomes_contradictory(ne, from, to)` and the batched
   `_ne_becomes_contradictory(ne, pairs::Dict)` variant.

2. **[src/expressions/index.jl](../../../SecondQuantizedAlgebra.jl/src/expressions/index.jl)**:
   wire those checks into both `change_index(QAdd, from, to)` and
   `change_index(QAdd, pairs::Dict)`. When a rename collapses an NE pair
   `(α, β)` into `(x, x)`, drop the term instead of silently relaxing
   the NE.

Without these, `scale`'s slot-folding silently strips NE constraints
between user-distinct atom indices, and cluster correlations
`⟨σ_i^{12} σ_j^{21}⟩` (with `i ≠ j` declared) mis-fold into single-atom
averages. This is the same bug class as the σᵢᵢ leak the prior commit
fixed, just at a deeper layer.

The QC plan adds a regression test that composes scale (slot-folding
under same-index assignment) with an NE-bearing source so the
contradictory-NE drop is exercised end-to-end.

## Design

### Convention

| Form              | Meaning                                       |
| ----------------- | --------------------------------------------- |
| Bare base (`j`)   | Symbolic, unmaterialised, free index          |
| Suffixed (`j_k`)  | Concrete materialised position, *always* with `_k` suffix even for k=1 |

The boundary between symbolic and materialised is a single SQA primitive:
`(i::Index)(k::Integer)`. Every site in QC that today mints a
"concrete atom-k" `Index` calls this primitive. There is no other
materialised-index spelling.

### Single source of truth (design invariant)

`_build_canonical_indices(eqs)` in
[src/completion.jl:75](../../src/completion.jl#L75) is the **sole
authority** for the per-subspace base-name mapping. All three passes
(`completion`, `scale`, `evaluate`) call it and consume the result
identically. There is no alternative way to derive "the canonical first
index for subspace s". If a subspace has no declared index, the
fallback (the first free index found in any operator on that subspace)
is also produced by `_build_canonical_indices` and shared across passes.

### Components

1. **SQA**: no new changes. `(i::Index)(k::Integer)` already exists.
   The `_ne_becomes_contradictory` edits already landed in SQA
   `redesign-v2`.

2. **[src/completion.jl](../../src/completion.jl)**, `_slot_rep`:
   - Delete the special-case for k=1.
   - Body becomes `_slot_rep(base, k) = base(k)` (or inline the SQA
     primitive at the call site and delete `_slot_rep` entirely).
   - `_build_canonical_indices` populates the slot table using
     `first_idx(k)` for every k, including k=1.

3. **[src/evaluate.jl](../../src/evaluate.jl)**, `_fresh_index`:
   - Replace the manual `Symbol(base_name, "_", k)` construction with
     `canon[b.space_index][1](k)` (apply the SQA primitive to the
     canonical first declared index for `b`'s subspace, where `canon`
     comes from `_build_canonical_indices`).
   - Side effect (intentional): a user who declares `j` first and
     threads `k_ind` through ops materialises *both* via the same
     primitive on `canon[…][1]`, so the resulting `Index` objects are
     `==` and dedup naturally.

4. **[src/scaling.jl](../../src/scaling.jl)**, `_min_slot_assignment`:
   - Continues to call `_slot_rep`, which now routes through SQA.
   - **Verification required (resolved):** the permutation-min
     comparator orders slot assignments by `(space_index, name)` via
     `Base.isless(::Index, ::Index)` ([SQA
     index_types.jl:62](../../../SecondQuantizedAlgebra.jl/src/expressions/index_types.jl#L62)).
     With every slot uniformly suffixed, the ordering is total
     (`j_1 < j_2 < … < j_K` by lex on suffixed name) and the
     permutation-min representative is unique. No code change required
     beyond routing `_slot_rep` through `base(k)`.

5. **`_aux_slot_rep` (Σ-prefixed)** and the `_alpha_rename_sources`
   workaround machinery: delete and revert callers to use `_slot_rep`
   directly. Specifically:
   - [completion.jl:438](../../src/completion.jl#L438) `_alpha_rename_sources`
   - [completion.jl:450](../../src/completion.jl#L450) `_observable_indices_by_space`
   - [completion.jl:460](../../src/completion.jl#L460) `_build_alpha_rename_map`
   - [completion.jl:491](../../src/completion.jl#L491) `_apply_rename`
   - [completion.jl:497-503](../../src/completion.jl#L497-L503) `_apply_rename_each`
   - [completion.jl:541](../../src/completion.jl#L541) `_aux_slot_rep`
   - [completion.jl:412 and 422](../../src/completion.jl#L412)
     callers in `_derive_for` revert to passing `eqs.hamiltonian`,
     `eqs.jumps`, `eqs.jumps_dagger` directly.
   Frees about 70 lines.

6. **`_canonicalise_avg_leaves`** at
   [evaluate.jl:664](../../src/evaluate.jl#L664), invoked at
   [evaluate.jl:131](../../src/evaluate.jl#L131):
   - The design closes the *materialised-state* mismatch case.
   - The post-pass may still handle the alpha-equivalence at completion
     time (see the "two-atom JC drift" testset and the
     "alpha-equivalent leaves are one state" testset at
     `indexed_meanfield_test.jl:176-209`).
   - **Action:** keep the post-pass in place for this design's
     implementation. File a follow-up to attempt deletion after the
     materialised path is unified and the alpha-equivalence path is
     reconfirmed to be subsumed (or to depend on it explicitly).
   - The follow-up does NOT block this design.

### MTK and codegen surface

- [mtk.jl:55](../../src/mtk.jl#L55) builds MTK Sym names from
  `idx.name`: `"_" * string(idx.name)`. After this design, MTK unknown
  names go from `_j` to `_j_1` for materialised states. This is
  user-visible in `unknowns(sys)`, `parameter_map(sys)`,
  `get_solution(sol, op, eqs)`, and printed error messages.
- [mtk.jl:20](../../src/mtk.jl#L20) concatenates NE pair names as
  `"$(a.name)neq$(b.name)"`. Cosmetic change: NE keys read
  `j_1neqj_2` instead of `jneqk`. No semantic impact.
- The state-literal rename happens entirely in the equations object
  before MTK sees them, so consistency is preserved: every leaf in
  every RHS uses the same `_<k>`-suffixed Index, MTK's identity-based
  binding matches by SQA equality, and codegen produces the same number
  of unknowns as before with renamed labels only.
- `correlation.jl` and `latexify.jl` do not reference `Index.name` at
  all (verified by grep). No render changes there.

### Test impact

This is a breaking change in user-visible state symbol names *for
materialised states only*. The known affected test sites:

- **`test/indexed_filter_cavity_test.jl:113-122`** (the 4 currently
  failing tests + 1 error): these are exactly the failures this design
  closes. They should pass after the rename.

- **`test/indexed_meanfield_test.jl:126-174`** ("indexed jumps: σ_jj
  equation has no σ²/σ³ leak after scale"): the test builds an
  expected symbolic RHS using `σ(2, 2, i)` where `i` is the user's
  declared index, then asserts
  `SymbolicUtils._iszero(SymbolicUtils.simplify(rhs - expected))`.
  After this design the scaled RHS uses `i(1)`. The fix is to update
  `expected` to use `σ(2, 2, i(1))` and similarly for `σ(1, 2, i(1))`,
  `σ(2, 1, i(1))`. **This is a convention update, not a test
  weakening.** The σ²/σ³-leak invariant is unchanged; only the symbol
  the test compares against has been renamed.

- **`test/indexed_meanfield_test.jl:176-209`** ("alpha-equivalent
  leaves are one state"): uses `find_missing` and `repr(m)` regex
  matches. The regexes already match `σ_i` patterns; verify they
  continue to behave correctly when matching against materialised
  states (likely already fine since the user's pre-materialisation
  indices stay bare).

- **`test/indexed_meanfield_test.jl:31-91`** ("two-atom collective JC
  drift"): contains `_is_user_derived(name)` accepting
  `Symbol(base, "_", k)` slot names. Already compatible.

- **Regression coverage for the NE-contradictory drop:** rather than
  add a new dedicated test, the regression is already exercised by:
  - The σ²/σ³-leak-after-scale test (lines 126-180), which uses the
    superradiant-laser setup with `J = [a, σ(1,2,i), σ(2,1,i),
    σ(2,2,i)]` at order=2 and asserts the dissipator's symbolic form
    via `Symbolics.IM` equality. Any incorrect NE relaxation under
    slot-folding would produce spurious σ²/σ³ terms here.
  - The filter-cavity per-Hilbert-space evaluate/scale commute tests,
    which previously failed at exactly the NE-relaxation bug class.
  A standalone test that constructs a minimal cluster with `ops =
  [a, σ(2,2,j)]` was attempted during implementation and turned out
  to be moot: at order=2 with this observable set the system closes
  on single-atom states only and never produces a 2-atom cluster
  correlation in `eqs_sc.operators`, so there is no fold-candidate
  to drop. The two existing tests reach the path; a third would not
  add coverage.

### What stays the same

- Operator-level semantics: scaled and evaluated equations describe the
  same physics as before.
- The cluster correlation work shipped earlier remains. The two-atom
  correlation that previously read `⟨σ_i^{12} σ_{i_2}^{21}⟩` (mixed
  bare and suffixed) now reads `⟨σ_{i_1}^{12} σ_{i_2}^{21}⟩`
  (uniformly suffixed) and is recognised as the same physical state by
  SQA equality.
- `MeanFieldEquations` struct and the user-facing `meanfield` /
  `complete` / `scale` / `evaluate` / `System` APIs.

## Acceptance criteria

1. `make test` reports zero failures and zero errors. The total test
   count must not regress versus the pre-change baseline minus the 4
   known filter-cavity failures (i.e., at least 695 passing, possibly
   more if regression tests are added).
2. No regression in `make jet` (run separately).
3. Only one site mints suffixed indices: SQA's
   `(i::Index)(k::Integer)`. Verify with
   `grep -rn 'Symbol(.*\.name.*"_".*,' src/` in QC `src/` returning no
   hits, and with all `_<k>`-suffix construction in QC going through
   `base(k)` or `_slot_rep(base, k) = base(k)`.
4. The six workaround helpers
   (`_alpha_rename_sources`, `_observable_indices_by_space`,
   `_build_alpha_rename_map`, `_apply_rename`, `_apply_rename_each`,
   `_aux_slot_rep`) are deleted. `_canonicalise_avg_leaves` is kept
   for this design's implementation with a comment explaining why
   (alpha-equivalence at completion time, see component 6).
5. The "indexed jumps: σ_jj equation has no σ²/σ³ leak after scale"
   regression continues to pass with `Symbolics.IM` symbolic-equality
   assertion (no string fallback), with the expected RHS updated to
   use `i(1)` per this design's convention.
6. Regression coverage for the NE-contradictory drop is verified by
   the existing σ²/σ³ leak test and the filter-cavity commute tests
   passing under the new convention. No new test is added; see Test
   impact for rationale.

## Out of scope

- New SQA changes. The `(i)(k)` primitive and the NE-contradictory
  edits are already in SQA `redesign-v2`.
- The cluster-correlation logic in `_canonical_key` and
  `_min_slot_assignment`. Unchanged. This design just lets those
  produce consistent symbols.
- Deletion of `_canonicalise_avg_leaves`. Deferred to a follow-up.
- Master-style `insert_index(., ., k::Integer)` substitution.
  Considered and rejected: SQA's existing `(i)(k)` primitive plus
  `change_index(., from::Index, to::Index)` covers the use case
  without adding a new SQA signature.
