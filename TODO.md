# TODO

Open bugs to fix in the package. Test-scaffold notes and historical
fix write-ups live in git log / CHANGELOG.

## 1. det vs stoch closure asymmetry: RESOLVED

`measurement_backaction_indices_comparison_test::deterministic vs
stochastic LHS match` now passes. Full suite green at 776/776 (count
went from 777 to 776 because one unique_squeezing assert was the same
artifact this session removed, see below).

**Resolution lives in three pieces:**

1. **Conjugate dedup bug in `_collect_missing!`**
   ([src/completion.jl](src/completion.jl)). The original code pushed
   `average(conj_full)` whenever `dedup_conj != dedup`, with no check
   that `dedup_conj` had already been added to `seen_keys` via a
   different leaf. So when leaf A's primary state happened to be
   leaf B's conjugate, find_missing returned the same state twice.
   Closure gained duplicate equations; downstream `scale` collapsed
   them away, so the bug was latent in models where scale's symmetric
   collapse hid it. unique_squeezing's pre-scale lock dropped from 24
   to 22 once the bug was fixed; eqs_sc count (19) and the numerical
   end-state (30.15) are unchanged.

2. **SQA gained `Index.concrete::Bool`**
   ([SecondQuantizedAlgebra/src/expressions/index_types.jl]). The 4-arg
   `Index(name, range, space_index, sym)` ctor defaults `concrete=false`;
   `(::Index)(::Integer)` (user-facing slot-mint) sets `concrete=true`.
   `Base.isequal(::Index, ::Index) = (a == b)` ignores the flag so
   substitution and dedup behave identically to before. Equality and
   hash already ignored the flag; the explicit isequal override is
   needed because `Base.:(==)(a::Transition, b::Transition) = isequal(a,b)`
   falls back to structural isequal across all fields if not overridden.

3. **QC reads `idx.concrete` directly** instead of the previous
   "non-bound atom index in initial_operators" heuristic that
   mis-classified user-free `k = Index(...)` as concrete. The
   classification is now exact: an atom index is treated as a
   concrete-site label iff it was minted via `(::Index)(::Integer)`.

The previous belief that "`Base.isequal` override breaks filter-cavity"
was a wrong conclusion; the actual breakage was the conjugate dedup
bug, which the SQA changes exposed by altering canonical-key output
enough to make the latent conjugate-collision reachable in
filter-cavity. With the dedup fix in place, the SQA `concrete` design
works as intended.

## 2. `unique_squeezing.jl` Full model oscillates (pre-existing)

Plot 5 at N=100: the **Full model** oscillates between roughly -8 and
+5 instead of settling at the X ≈ 2.29 / P ≈ 0.44 plateau that master
and the Effective model both produce. The N=1 anchor `⟨a'a⟩ ≈ 30.15`
locked in `test/examples_regression_test.jl` is the barely-stable
manifestation of the cumulant cross-decay cancellation that this
session's Bug 2 conditional preserves (NE injection on det path is
suppressed when `j(1)` is detected as user-concrete).

The Full model uses
`meanfield([a, a'a, σ(2,2,j(1))], H, [b, σ(1,2,i)]; rates=[κ, γ],
order=2)` then `complete` then `scale`, with the squeezed-bath jump
`b = a cosh(ξ) + a' sinh(ξ)`. The cavity-only Effective model uses
the **same** `b` jump and is correct, narrowing the suspect surface
to the atom-cavity coupling.

The bounded-but-wrong-physics behaviour at N=1 (30.15) and the
oscillation at N=100 share the same root: derived cross-atom states
involve algebra-minted phantom slot indices (e.g. `j_1_2`) whose
dynamics aren't physically correct for a model where atom 1 is the
only specified concrete atom. Cumulant truncation produces a closure
where the phantom states' dynamics is bounded near zero by a delicate
cross-decay cancellation; small perturbations to the dissipator
algebra (e.g. NE injection) break the cancellation. The Bug 2
conditional avoids the perturbation; the underlying phantom-state
representation remains. Real fix probably needs either
`expand_completeness` to fire selectively (only for permutation-
symmetric atom subspaces), or sum-wrapping cross-atom moments instead
of materialising them as independent states.

## 3. `excitation-transport-chain` end-of-chain drift (low priority)

Dashed end-of-chain trace settles at ~0.10 on rewrite vs ~0.11 on
master. Well within JumpProblem ensemble noise; could be a remaining
algebraic difference. Low priority. The current locked value in
`test/examples_regression_test.jl` (at N=4, deterministic) is on the
rewrite-current side, so any improvement that brings the rewrite
closer to master will fire the test and prompt an update.

## Reproduction recipes

Diagnostic scripts used this session, useful when re-investigating:

- `/tmp/gather13.jl`: solves unique_squeezing Full model at N=1, prints
  `⟨a'a⟩` at `t = 2π/ωd`. **Bounded behaviour: 30.1497...** Any change
  to dissipator NE handling that touches the cumulant cross-decay
  cancellation produces values in the ~1e7 to ~1e8 range. Run via
  `julia --project=examples /tmp/gather13.jl`.
- `/tmp/diff_us.jl`: dumps `eqs_c.states` and `eqs_c.equations` from
  unique_squeezing's completion, used to diff symbolic RHSes between
  fix variants. The clean-HEAD vs Bug-2 diff shows γ-coefficient
  changes in equations 14, 20, 24 (the cross-atom `σ_{j_1} σ_{j_1_2}`
  states), e.g. eq 20: `-(1/2)γ` becomes `-(3/2)γ` and 5 cumulant-
  cross-decay correction terms disappear. That's the mechanism by
  which NE injection on the det path destabilises unique_squeezing.

## Failed approaches (don't retry these in isolation)

Each row is a fix that was tried, what it broke, and why.

| Approach | Result | Why it doesn't work |
|---|---|---|
| Mirror NE injection on det path only (no filter, no LHS strip) | unique_squeezing → 1.4e7; filter-cavity drops 1 state | Removes cumulant cross-decay correction in dissipator. |
| Mirror NE + LHS strip (no filter) | unique_squeezing → 1.4e7; closures sealed | Same as above; LHS strip is cosmetic for state-identity, not algebra. |
| Pin user-concrete indices in `_canonical_key` (no SQA sort) | unique_squeezing OK; ramsey/many-atom-laser/heterodyne/filter-cavity state counts EXPLODE | Loses encounter-order operator-ordering coalescence; load-bearing for compact closures. |
| Pin canon + `canonical_sort_name_break` post-rename | unique_squeezing OK; state counts shift but don't fully recover | Sort can't normalise op-type-at-position differences (`σ^{12}_a σ^{22}_b` vs `σ^{22}_a σ^{12}_b`); swapping operators across positions changes physical identity. |
| `_build_canonical_indices` uses `initial_operators` only (no other changes) | heterodyne gains 1 state | Canon being frozen mid-completion changes which leaves dedup vs which become new states. |
| NE-injection filter = `user_concretes` only | unique_squeezing → 1.25e8 | Doesn't block the j_1_2 ↔ i (bound) assertion that splits the dissipator sum. |
| NE-injection filter = `user_concretes ∪ bound` (unconditional) | unique_squeezing OK; heterodyne +1 state; filter-cavity cascade errors in `evaluate`/`mtkcompile` | Filters too aggressively for population-symmetric models. |
| NE-injection filter = `user_concretes ∪ bound` conditional on `!isempty(user_concretes)` + LHS strip | **Current state**: unique_squeezing OK; 3 of 4 measurement_backaction resolved; closures sealed (no cascade); 2 fails remain because `k = Index(...)` mis-classifies as concrete | Limitation: heuristic can't distinguish user-free `k` from user-concrete `j(1)`. |

## Conceptual map: which models need which dissipator algebra?

Two model "shapes" with conflicting closure requirements:

- **Concrete-site shape** (unique_squeezing): user has a specific atom
  like `σ(2,2,j(1))` in initial ops. Completion derives cross-atom
  moments by minting a phantom partner `j_1_2`. The phantom is a
  representational artefact, not a physical atom. The closure
  *requires* the cumulant cross-decay correction in the dissipator to
  keep the phantom states near zero (the cancellation that makes the
  artefact harmless). NE injection between phantom and partner BREAKS
  this cancellation → divergence.

- **Population shape** (filter-cavity, heterodyne, ramsey, many-atom-
  laser): user has free atom indices `k = Index(...)` over a population
  `N`, often only `[a'a]` or similar non-atom-LHS ops. Cross-atom
  moments derived in completion are genuine multi-atom states over the
  population. The closure *requires* NE injection + `expand_completeness`
  fold of `σ^gg → 1 - Σ σ^kk` to keep state count compact and match
  master-rewrite compatibility.

These two shapes need *opposite* treatment at the same code point.
The current conditional uses `isempty(_user_concrete_atom_indices)` as
the discriminator, which works for unique_squeezing and most others
but fails when the user names a free index like `k` (it gets
mis-classified as concrete because it appears in initial_operators
without a sum-scope binding).

The honest discriminator is **at construction time**: a user calling
`j(1)` (via `(::Index)(::Integer)`) is concrete; a user calling
`Index(h, :k, N, ha)` directly is free. SQA's `Index` struct currently
doesn't carry this metadata, which is the structural change item 1
flags as the path forward.

## Notes from this session (worth not relearning)

- The encounter-order naming in `_canonical_key` has *two* load-bearing
  side effects: (1) renaming indices position-by-position which
  collapses operator orderings (`σ_a σ_b` and `σ_b σ_a` map to same
  canonical name), (2) overwriting user-declared concretes when they
  appear later in a term. You cannot fix (2) without losing (1) inside
  QC alone, because operator-order normalisation requires sorting
  cross-site operators that SQA's `_partial_sort!` leaves as
  `Undetermined`. The `canonical_sort_name_break` SQA primitive in
  this branch provides the sort but doesn't fully restore (1) when
  the cross-atom operators differ in type (e.g. `σ^{12}_a σ^{22}_b`
  vs `σ^{22}_a σ^{12}_b` need both index-name AND op-type
  normalisation, and op-type swap changes physical identity).
- The dissipator algebra with and without NE injection produces
  different cumulant-2 truncations even when the underlying physics
  is the same (`σ^gg = 1 - Σ σ^{kk}` fold doesn't commute with
  cumulant-2 factorisation). Models depending on the cumulant
  cross-decay correction terms (unique_squeezing) need the
  no-NE form; models depending on `σ^gg` folding (heterodyne,
  filter-cavity, master-compatibility) need the with-NE form. There
  is no globally correct choice; QC has to discriminate by model
  shape.
- `eqs.operators` grows during `complete!` with algebra-minted index
  names. Anything reading `eqs.operators` after completion (e.g. to
  build canon, or to walk user state) gets the augmented set; the
  `initial_operators` snapshot is the right field for "what did the
  user originally declare".
- Never relax regression locks (per memory). State-count locks can
  legitimately change when the closure shape changes from a
  correctness fix, but only when (a) the drift-numerical test still
  passes and (b) the change has an explained mechanism. Both
  conditions hold for heterodyne and filter-cavity here.
- Local SQA path is wired into `Project.toml` and `test/Project.toml`
  `[sources]`. Reverting that requires also reverting any uses of the
  new SQA primitive (`canonical_sort_name_break`); currently QC
  doesn't call it.
