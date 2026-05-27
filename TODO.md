# TODO

Open work items for the v1 rewrite. Each entry names the failure mode and where
it shows up so the fix can be verified end to end.

## Examples status

All 13 examples now run end to end on the rewrite branch. Side-by-side
comparison run via `tools/run_examples.jl` against the master worktree
(Julia 1.11, QC v0.4.3) and the rewrite branch (Julia 1.12, QC v0.5.0).
Report at `tmp/report.html`.

### Identical to master

- `single-atom-laser-spectrum.jl`
- `mollow.jl`
- `many-atom-laser.jl`
- `optomechanical-cooling.jl`
- `ramsey_spectroscopy.jl`
- `waveguide.jl`
- `excitation-transport-chain.jl` (deterministic + ensemble; tiny ~10%
  steady-state drift on the dashed end-of-chain trace, see below)
- `retrodiction_homodyne.jl` (plots 1, 2, 4 identical; plot 3 amplitude
  difference is SDE RNG)
- `superradiant_laser_indexed.jl` (time evolution AND spectrum correct)
- `unique_squeezing.jl` (after `i(1) -> j(1)` lookup fix)
- `heterodyne_detection.jl` (deterministic + SDE; closes at 12 unknowns
  matching master after the get_adjoints / distinct-atom fixes below)
- `cavity_antiresonance_indexed.jl` (antiresonance dip restored; the
  literal-key fix in `_canonicalise_avg_leaves` resolved the apparent
  dipole-dipole bug, which was actually an alpha-rename collapse of
  distinct per-atom states)
- `filter-cavity_indexed.jl` (filter populations stable at ~0.01,
  intensity spectrum is a clean Lorentzian; after canon-slot lookup
  update in the example file)

### Outstanding test failures (regressions from the dedup fix)

- `measurement_backaction_indices_test.jl::indexed measurement backaction:
  filter cavity drift agrees det vs stoch` and
  `measurement_backaction_indices_comparison_test.jl::measurement_backaction_indices_comparison:
  deterministic vs stochastic LHS match`: assert `det_lhs ⊆ stoch_lhs`.
  After the `_dedup_key_strip_free_ne` change in `find_missing`, det
  derives an extra `⟨σ_j_2₁₁⟩` equation that stoch does not. Root cause:
  the operator-level derivation for `⟨σ_j_2₂₂ * σ_j_2_2₂₁⟩` produces a
  factor-3-different Γ coefficient and uncollapsed 3-op cross-atom
  terms between det and stoch paths (stoch consolidates, det does not).
  Some asymmetry in how SQA's `_canonicalize!` / `assume_distinct_index`
  fires between `_meanfield_deterministic` and `_meanfield_noise`. The
  symbolic discrepancy then propagates through cumulant truncation and
  find_missing into one extra (extraneous) σ¹¹ state on the det side.
  This is not a numerical bug in the example outputs, only the test's
  subset assertion. Needs a focused SQA-level investigation.

### Minor numerical divergence

- `unique_squeezing.jl` plot 5: Full model curves are correct (X ≈ 2.29,
  P ≈ 0.44, matching master). Effective model curves are off by ~20%:
  X_a ≈ 1.88 (vs master 2.3), P_a ≈ 0.62 (vs master 0.45). The Full and
  Effective models should overlap exactly (the effective adiabatic
  Hamiltonian `H_a = Hf - N*gΩ*(a+a')²` is supposed to reproduce the
  Dicke model in the low-excitation limit). Verified to be PRE-EXISTING
  on the rewrite branch at clean HEAD (commit 116752a), not introduced
  by the current Spectrum/dedup work. The Effective model uses
  `meanfield([a, a'a, a*a], H_a, [b]; rates=[κ], order=2)` with no
  indexed completion or scaling, so the bug is upstream of all the
  cross-atom completion logic. Suspects: (a) the squeezed-bath jump
  handling with cosh/sinh ξ rates, (b) cumulant-order-2 truncation of
  cross-cavity products under the squeezed bath, or (c) the mtk codegen
  substitution of `⟨a'·a'⟩` as `conj(state_aa)`. Needs a focused
  derivation-vs-symbolic-equation comparison (Mathematica check).
- `excitation-transport-chain.jl`: dashed "end of chain" trace settles at
  ~0.10 on rewrite vs ~0.11 on master. Smaller than pre-fix gap and well
  within the noise of a JumpProblem ensemble; could be a remaining
  algebraic difference. Low priority.
- `retrodiction_homodyne.jl`: docs build entry in `docs/make.jl` is still
  commented out to match master (PR #266 history; suspected 5 x 10^4 step
  SDE solve being too heavy for docs).

## Fix landed: heterodyne SDE blowup (conjugate-pair + cross-atom consolidation)

The heterodyne SDE blowup had two intertwined root causes, both now fixed.

**1. Doubled state count from conjugate-pair handling.** Master's
`find_missing` filtered conjugates of explicit states out of the missing
pool (effectively `get_adjoints=false` semantics, despite the surface
default of `true`). The rewrite's `find_missing` literally pushed both
members of every conjugate pair into `missing_states`, so the laser
closure produced 18 unknowns instead of 12. The 6 extras were tracked as
separate ODE/SDE state vars whose drifts preserved the
`u_i(t) = conj(u_j(t))` invariant exactly but whose diffusion columns
evaluated `g[i]` and `conj(g[j])` on independent floating-point paths;
roundoff drifted the invariant and `Ng = 2.3e5` amplified the slip into
exponential SDE blowup within ~100 steps at high N.

Fix:

- `src/completion.jl::find_missing` defaults to `get_adjoints=false` and
  emits one representative per conjugate pair. `complete!` propagates the
  new default.
- `src/mtk.jl::_conj_substitution_dict` walks RHS leaves directly and
  matches conjugates via a **permutation-canonical key**
  (`_perm_canon_key`): stable-sort QTerm ops by `(acts_on, string(op))`
  so cross-atom commuting ops on the same Hilbert subspace collapse to a
  deterministic order. Without this the raw adjoint reverses operator
  order (e.g. state `⟨σ_k₁₁₂ · σ_k₂₂₂⟩` has conjugate
  `⟨σ_k₂₂₂ · σ_k₁₂₁⟩` while the RHS literal is
  `⟨σ_k₁₂₁ · σ_k₂₂₂⟩`) and `_safe_substitute`'s literal-key dict misses
  the leaf, leading `mtkcompile` to refuse with "Brownian appears
  non-linearly".

**2. Implicit `σ^11` leaked into the closure via cross-atom products.**
Pre-fix, `expand_completeness` only fired for atomic ground-state
projectors. A Lindblad commutator like
`[σ_j_21, σ_k_22 · σ_k_2_22] · σ_j_12` produced
`σ_k_12 · σ_k_2_22 · σ_k_21` (3 ops, two on atom k, one on atom k_2)
which SQA refused to reorder per the "Undetermined free-index pairs
stay put" invariant. The cumulant truncation then factorised this
across atoms and materialised the implicit `σ_k_12 · σ_k_21 = σ_k_11`
as a separate state. Master applies the same algebra but its
`*`/`commutator` knew to reduce same-atom-cross-atom products because
the algebra layer was simpler; the rewrite, building on SQA's stricter
distinct-index policy, does not.

Fix:

- `src/meanfield.jl::_assume_distinct_atom_indices` (called from
  `src/completion.jl::_derive_for`, plumbed via a `distinct_indices`
  kwarg through `_meanfield_deterministic` / `_meanfield_noise` /
  `_build_op_drift`) injects NE pairs over the auto-minted slot indices
  that `_canonical_key` creates for cross-atom moments (e.g. `k_1`,
  `k_2`), then re-canonicalises via `SQA.assume_distinct_index`. This
  unlocks SQA's cross-atom commutation + same-site collapse, so
  `σ_k_1₁₂ · σ_k_2₂₂ · σ_k_1₂₁` reduces to
  `σ_k_2₂₂ - σ_k_1₂₂ · σ_k_2₂₂` via `σ_12 · σ_21 = σ_11 = 1 - σ_22`
  applied through `expand_completeness`. Critically, this is only
  applied inside `_derive_for` (when minted slot indices are introduced
  by completion), NOT in the top-level `meanfield` call where the user
  may legitimately want diagonal (`i = j`) contributions in their
  H-sums.

Additional supporting changes:

- `src/completion.jl::_canonical_key` now stable-sorts `encountered`
  free indices by `(space_index, string)` so the slot assignment is
  independent of operator-product order. Without this, adjoint-reversed
  cross-atom products got swapped index renames vs the originals.
- `_canonical_key` ALSO strips NE pairs on present indices in the
  returned QAdd, so the LHS state representation is canonical and
  doesn't double-up via NE metadata. `_strip_present_ne` does the work.
- `src/completion.jl::_canonical_dedup_key` builds a permutation-
  canonical, NE-blind string key for `seen_keys` dedup, matching the
  signature used by `src/mtk.jl::_perm_canon_key` so RHS leaves and
  states with permuted same-Hilbert-subspace ops collapse to the same
  key.
- `_collect_missing!` deterministically tie-breaks rep choice by
  lex-smaller `_canonical_dedup_key` (with a fall-through to the
  conjugate rep if `filter_func` rejects the lex-smaller one). Without
  this, det-vs-noise paths over the same physical system picked
  different conjugates as the tracked state, breaking subset invariants.
- `src/mtk.jl::_collect_conj_subs!` falls back to direct
  `_perm_canon_key` match (substituting the state's var) before trying
  the conjugate match, so RHS leaves carrying NE metadata that the
  matching state lhs doesn't carry still get resolved at codegen time.

Net effect on the heterodyne laser: rewrite now closes at 12 unknowns
exactly matching master, the SDE solver returns `ReturnCode.Success`,
and max photon number (`~591` for the rewrite, `~442` for master) is
within the RNG-trajectory variance.

Regression tests (in `test/completion_test.jl`):

- `find_missing default get_adjoints=false: one rep per conjugate pair`
  asserts the laser closure is ≤12 unknowns and `mtkcompile` succeeds.
- `_build_op_drift consolidates cross-atom products` asserts no σ^11
  state survives in the closure and the size matches master's 12.

Additional fix to make `Spectrum` survive the new dedup policy:

- `src/correlation.jl::_spectrum_kernel` walks every RHS leaf and resolves
  conjugate-of-state averages via `_perm_canon_key` (matching the codegen
  policy in `src/mtk.jl::_collect_conj_subs!`). Without this the
  steady-state Jacobian dropped to rank 1 (out of 6) at order 4 with
  the phase-invariant filter, because conj leaves like
  `⟨a' * a * a (ancilla) * σ_21⟩` (touching original + atom + ancilla)
  were skipped by `_ambient_avgs` and silently zeroed in `zero_sub`,
  collapsing every column whose state's conjugate appeared on the RHS.

## Fix landed: indexed atom-cavity coupling regression

The catastrophic regression (superradiant flat output, unique_squeezing
blowup, heterodyne 4-orders-of-magnitude photon scale, cavity_antiresonance
missing dip, filter-cavity exploding populations) had one root cause in
`complete!`: `_build_canonical_indices` included Hamiltonian-bound sum
indices in the canonical-name pool. When `find_missing` then renamed a
missing state's free index onto an H-bound name (e.g. user-declared `j`
mapped to canon `i`, where `i` was H's bound sum index), the subsequent
`meanfield` call silently dropped cross-atom commutator terms
(`Sigma_{k != i} g_k <sigma_k^{21} sigma_i^{12}>`) because the free `i`
collided with H's bound scope.

Fix in `src/completion.jl`:

1. `_build_canonical_indices` filters out H/J-bound indices from the
   canonical slot pool (`_bound_indices` helper).
2. When the pool is empty after filtering (every user-declared index on
   a space is bound), mint a single deterministic fallback as the
   lex-first declared index's successor (`first[1](2)`, e.g. `:i`
   becomes `i_2`).
3. `_canonical_key` mints slot `k > len(canon)` as `canon[1](k)` so
   higher-order correlations always get disjoint canonical names.
4. `_canonical_key` strips per-term NE pairs that reference an index not
   present in the operator atoms (these are constraints inherited from
   a deeper derivation and prevent dedup).
5. `_derive_for` alpha-renames LHS indices that still happen to be bound
   to fresh successors before calling meanfield, then undoes the rename
   on the derived equations. The undo (`_apply_undo`) handles sum-bound
   collisions by first relabelling any bound name in the result that
   would clash with the undo target.
6. `src/correlation.jl`'s `_complete_ancilla!` accepts a `parent_canon`
   kwarg and threads `_build_canonical_indices(eqs0)` so the correlation
   meanfield reuses the parent's atom-index slot. Without this the
   ancilla-completion would mint a different canonical name and ambient
   steady-state lookups against `eqs0.states` failed, collapsing the
   `Spectrum` to a constant.
7. `find_missing` accepts an optional `canon = nothing` kwarg that
   defaults to its own `_build_canonical_indices(eqs)` call, enabling
   parent-canon threading from `_complete_ancilla!`.

Fix in `src/evaluate.jl`:

8. `_canonicalise_avg_leaves` tries a literal (no-alpha-rename) key
   FIRST, then falls back to the canon-equivalent key. The literal key
   prevents distinct per-atom states `<sigma_{j_1}>` and `<sigma_{j_2}>`
   from being alpha-collapsed onto canon slot 1 (which would silently
   pick the wrong atom's value when atoms are distinguishable).
9. When both keys miss, `_canonicalise_avg_leaves` rewrites the leaf to
   its NE-stripped form `average(literal)` so downstream MTK
   substitution (which builds keys from a non-NE operator structure)
   can resolve the leaf via `_conj_substitution_dict`. This handles the
   case where a cross-atom sum inherits an irrelevant NE pair from its
   parent derivation.

Example files updated to match the new canonical slot convention:

- `examples/superradiant_laser_indexed.jl` (`sigma(2, 2, i(1))` becomes `sigma(2, 2, j(1))`)
- `examples/unique_squeezing.jl` (same)
- `examples/filter-cavity_indexed.jl` (`b(k::Integer)` uses `i(2)(k)`
  instead of `i(k)` so per-filter lookup hits the canon-fallback slot
  `i_2`)

Test updated to match:

- `test/indexed_meanfield_test.jl` "sigma_jj equation has no sigma^2 leak
  after scale" (canonical slot is `j(1)` instead of `i(1)`).

All package tests pass after these changes.
