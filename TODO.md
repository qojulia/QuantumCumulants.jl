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

2. **`_collect_atom_indices_set!` discriminator**
   ([src/completion.jl](src/completion.jl)). The "is this index a
   user-pinned slot?" check now reads only QC state: any atom-space
   index that appears in `eqs.initial_operators` and isn't H/J-bound
   qualifies. This covers both free `j` (from `Index(h, :j, N, ha)`)
   and slot-minted `j(1)` shapes uniformly.

   (An earlier attempt added an `Index.concrete::Bool` flag to SQA to
   drive the discrimination at construction time. That flag turned out
   to be the wrong discriminator: the unique_squeezing example uses
   free `j`, which would have read `.concrete = false`, so the flag
   could not distinguish "user named this index" from "algebra minted
   this index" the way unique_squeezing needed. The flag was removed
   from SQA; the discriminator now reads `initial_operators` directly.)

The previous belief that "`Base.isequal` override breaks filter-cavity"
was a wrong conclusion; the actual breakage was the conjugate dedup
bug. With the dedup fix in place, no SQA-side change is needed at all.

## 2. `unique_squeezing.jl`: FIXED (free-`j` shape was broken)

The example file uses `meanfield([a, a'a, σ(2,2,j)], H, J; ...)` with
free `j` (declared via `Index(h, :j, N, ha)`); the TEST file uses
`σ(2,2,j(1))` (slot-minted via `(::Index)(::Integer)`). These two
shapes were not equivalent in the rewrite: an early attempt at a
construction-time tag (an SQA `Index.concrete` boolean flag) routed
the two shapes through different paths. The slot-mint path
correctly skipped NE injection in `_derive_for` (because NE blocks
the dissipator's cumulant cross-decay correction that bounds the
dynamics), while the free-`j` path triggered NE injection because
the construction-time tag was `false`, so `_user_concrete_atom_indices`
returned an empty set for it.

Fix: read QC state, not the construction-time tag. Any atom-space
index the user named in `initial_operators` that isn't H/J-bound is a
"user-pinned slot" regardless of whether the user wrote `j` or `j(1)`.
Both representations now yield the same closure shape and the same
numerical result (X=2.29, P=0.44 at N=100, matching master and the
Effective model).

Re-investigation (test parameters, j(1) shape):

- N=1: `⟨a'a⟩(t_end) = 30.149750050640776` with test parameters
  (`tend = 2π/ωd`, `reltol = 1e-8`). The locked value 30.149750050640534
  is matched to within `rtol = 1e-3`. The N=1 Full model legitimately
  differs from the Effective model (which is the high-N adiabatic limit
  at N_=69), so comparing the two at N=1 is not meaningful.
- N=100: X plateau `2.292`, P plateau `0.437`, settling smoothly over
  `[0, 4π/ωd]`. These are exactly the values that the previous session
  claimed only master and the Effective model produce. The rewrite
  reproduces them.
- Scale already emits the correct `(N-1)` prefactor on cross-atom
  moments (e.g. eq 6 of the post-scale system has
  `-0.5(-1+N)*g*im*⟨σ_{j_1_1}^{12} σ_{j_1_2}^{21}⟩`). The pre-scale
  user-concrete `j_1` is gracefully renamed to canonical slot `j_1_1`
  and the post-scale model is fully permutation-symmetric, as it
  should be (the Hamiltonian and jumps are symmetric across all N
  atoms; the user-concrete `j(1)` in `σ(2,2,j(1))` is just a readout
  label).

My initial "no fix needed" diagnosis was wrong because I had been
running the j(1) (test) variant via the MCP. The actual example file
uses free `j`, which the user-concretes detector then mis-classified
as "no user-pinned atoms" and incorrectly triggered NE injection.
The fix is one line in `_collect_atom_indices_set!`.

## 3. Det vs stoch NE asymmetry: RESOLVED (2026-06-02, structural discriminator)

**Resolution.** The path-asymmetric NE policy (det skips NE when the user
named an atom index; stoch always injects) is replaced by ONE policy selected
by a STRUCTURAL property of the channel set, applied identically to the
deterministic and noise paths:

> A system whose dissipators include a **dephasing channel** (a diagonal
> atomic jump `σ^{αα}`) is a *population* system and gets cross-atom NE
> injection (so SQA folds `σ^gg = 1 − Σ σ^kk` and the order-2 closure stays
> compact/bounded). A system with **no** dephasing channel (e.g.
> unique_squeezing's single `σ^{12}` decay) is a *concrete-site* system and
> skips NE (retaining the dissipator's cumulant cross-decay correction).

Implementation: `_has_dephasing_channel(jumps)` in
[completion.jl](src/completion.jl) inspects the (flattened) jump list once for
an indexed `σ^{αα}`; `meanfield` and both `_derive_for` methods pass it as the
`cross_ne` flag to `_build_op_drift`, which injects the NE per derived operator
(`_distinct_atom_indices([op])`, see Order-independence below). The
`user_concretes` heuristic (`_user_concrete_atom_indices` /
`_collect_atom_indices_set!`) is **deleted**.

Why this works where the heuristic was fragile: the old `isempty(user_concretes)`
discriminator mis-read a free population index like `k` (it appears in
`initial_operators` without a sum binding) as "concrete". The structural
signal does not depend on index naming. Under it, superradiant and
heterodyne-det FLIP from the heuristic's "concrete" to "population" and every
regression lock + numeric still holds, confirming the structural classifier is
the correct one.

Crucially this is NOT the per-channel "drop the diagonal" idea (which deleted
physical dephasing, see the reverted attempt below). The policy is the existing
post-hoc population NE, which drives the `σ^gg` fold and is a no-op on
already-emitted diagonals; `χ`/`ν` stay in every closure.

**Validated:** full suite `818/818`, JET clean (5/5). unique_squeezing N=1 =
30.1497; heterodyne SDE (seed 2, `T_end=0.1`) bounded (max ≈ 591); heterodyne
drift-only bounded; superradiant / filter-cavity / unique_squeezing closure
shapes unchanged (23-20 / 54-48 / 22-19). Tests:
`measurement_backaction_indices_comparison` is tightened to assert det and stoch
reach the SAME closure (`det_lhs == stoch_lhs`) AND every shared drift agrees;
new `test/dephasing_discriminator_test.jl` locks the discriminator classification
and the free-`j`-vs-slot-`j(1)` readout-shape invariance (the exact case the old
heuristic mis-classified).

**Order-independence (also fixed).** The NE is injected PER DERIVED OPERATOR:
`_build_op_drift` computes `_distinct_atom_indices([op])` from each op alone,
not from the whole `complete!` batch. Previously the batch-wide distinct set
let a sibling op's cross-atom indices leak NE into an unrelated equation (e.g.
`⟨aa⟩`'s coherent sum `Σ_j⟨aσ_j⟩` picked up a spurious `(j≠k, j≠k_2)` because
`k, k_2` came from a cross-atom moment derived in the same batch), and det/stoch
batch differently, so a handful of drifts differed in sum-split form. With
per-op injection det and stoch derive byte-identical drifts (the comparison
test now asserts ALL shared drifts agree, not a sample).

What did NOT work (kept for the record; see detail below): a single policy for
both regimes. The concrete-site and population regimes genuinely need opposite
diagonal treatment under order-2 truncation, and the only way to serve both
with one policy is master's lazy-moment / deferred-sum model, which is ruled
out for this rewrite. The structural discriminator instead selects the correct
one of the two existing (working) policies per system, robustly.

### Historical: path-asymmetric framing (superseded by the structural discriminator above)

The current rewrite used path-asymmetric NE injection (det path skips
NE when user named an atom index in initial_operators; stoch path
always injects NE). Both unique_squeezing and heterodyne SDE stay
bounded with this asymmetry. Master keeps both bounded with what
looks like a unified NE convention; the asymmetry is therefore
algebraic rather than physical and worth resolving.

### INVESTIGATION (2026-06-01): two layers — a σ^gg fold leak AND an irreducible truncation asymmetry

Investigated whether §3 is really an SQA problem. Verdict: the
"fix lives in SQA's `_accumulate_with_diag!` / `cumulant_expansion`
interaction" framing (the SUPERSEDED section below) is the wrong target
— `_accumulate_with_diag!` is algebraically correct. But the simpler
"it's purely a completeness-fold coverage bug" hypothesis is ALSO
wrong: it was tried in source and the heterodyne SDE still diverged.
The asymmetry has two distinct layers.

**Layer 1 — a real σ^gg fold leak (cosmetic, fold-fixable).** Direct
reproduction on an isolated cross-atom moment `⟨σ²²_{k1} σ²²_{k2}⟩`
under the heterodyne channel set `J=[σ¹²_j, σ²¹_j, σ²²_j]` (summed over
`j`):

1. The diagonal split is UNCONDITIONAL — `Σ_j` construction
   (`_accumulate_with_diag!`) emits diagonal self-terms regardless of
   NE. `assume_distinct_index` on an already-built sum is a no-op on
   the emitted diagonals; it only changes NE metadata. "NE toggles the
   diagonal split" (Attempt #2's mechanism story) is NOT what happens.
2. Single channel on a population moment → off-diagonal vanishes by
   commutation → NE vs no-NE cumulant(2) diff EXACTLY 0.
3. Full channel set → NE and no-NE differ, and the diff is entirely
   ground-state projectors (no-NE carries `⟨σ²²_{k1} σ^{11}_{k1}⟩`,
   `⟨σ^{11}⟩`, `−4⟨σσ⟩`; NE is pure `σ²²` with `−6⟨σσ⟩`).
4. `SQA.expand_completeness` applied to both forms makes the cumulant(2)
   of THIS term identical. Confirmed at the system level: with
   `get_adjoints=true`, the no-NE heterodyne closure leaks an explicit
   `⟨σ_k₁₁⟩` state (19 states vs NE's 18). It enters during cumulant
   factorisation: `_prod_ops([σ²¹_k, σ¹²_k]) = σ^{11}_k`, multiplied
   AFTER `_build_op_drift`'s one-shot `expand_completeness`, so never
   folded. Wrapping `_prod_ops` in `expand_completeness` removes the
   leak (19→18, no `σ^gg`) and is fast (0.025s vs 0.36s; the in-loop
   cost is negligible).

**Layer 2 — an irreducible cumulant-truncation asymmetry (NOT
fold-fixable).** After the `_prod_ops` fold, diffing the NE vs no-NE
heterodyne closures (`get_adjoints=true`, 18 states each, same state
set, no `σ^gg`) leaves EXACTLY ONE differing equation: the cross-atom
population moment `⟨σ²²_k σ²²_{k'}⟩`. The no-NE form keeps the diagonal
self-dephasing `(1−⟨σ²²⟩)` terms (the `j=k`, `j=k'` contributions of
the `σ²²_j` dephasing channel — physically complete); the NE form drops
them. These are two genuinely different order-2 closures of the same
operator, and `expand_completeness` does NOT reconcile them.

**Decisive experiment (reverted).** Applied BOTH the `_prod_ops` fold
AND symmetrised the noise path (`_derive_for(::NoiseMeanFieldEquations)`
to mirror the det rule, no-NE when `user_concretes` non-empty —
heterodyne's free `k` counts as user-concrete). Closure state sets
matched the NE config (14 eqs, no `σ^gg`). But the heterodyne SDE
DIVERGED: `ReturnCode.Unstable`, `⟨a'a⟩ → NaN` (suite went 785/2 →
782/5, the 3 SDE-bounded tests failing). So the single differing
`⟨σ²²σ²²⟩` equation is load-bearing: the NE form (dropping the dephasing
channel's diagonal self-terms) is the order-2 truncation that stays
BOUNDED; the physically-complete no-NE form is unstable under order-2
cumulant truncation.

**Conclusion.** It is NOT primarily an SQA problem, and NOT merely a
fold-coverage bug. It is an intrinsic cumulant-truncation STABILITY
tradeoff: for a dephasing-type channel, the order-2 closure that keeps
the diagonal self-action diverges, and dropping it (via NE) is what
bounds the SDE; for unique_squeezing's single decay channel the
opposite holds (NE drops the only damping → divergence). Neither
choice is "the algebra being wrong"; the path-asymmetric NE workaround
is effectively SELECTING the stable truncation per dissipator shape.
Master's "unified convention" corresponds to the NE-like choice applied
uniformly (it kept sums symbolic through `cumulant_expansion` and
expanded only at `evaluate`, so the truncation saw the full sum as one
unit). A genuine unification would need QC to pick the stable order-2
closure structurally (e.g. defer the diagonal split like master's
`evaluate`, or detect dephasing-type channels), not a fold.

### Redesign blueprint (2026-06-01): the defer-the-split path

Pinned the exact mechanism that makes master stable where the rewrite's
no-NE form diverges. It is an ORDER-OF-OPERATIONS difference on the
dissipator's diagonal (`j=k`) contribution:

- **Rewrite (eager):** SQA's `*`/`Σ` runs `_accumulate_with_diag!` at
  Σ-construction, so the diagonal substitution `j→k` happens at the
  OPERATOR level — `σ²¹_j σ²²_k σ¹²_j |_{j=k}` COLLAPSES via Transition
  algebra to a lower-order operator FIRST, and `cumulant_expansion`
  truncates that. → "collapse-then-truncate".
- **Master (deferred):** `cumulant_expansion(::IndexedAverageSum)`
  (master `src/cumulant_expansion.jl:195`) keeps `Σ_j` SYMBOLIC and
  truncates the generic-`j` inner moment (`j` treated distinct from the
  LHS atoms); `evaluate` (master `src/index_average.jl`) materialises
  `Σ_j` LAST, substituting `j→k` into the already-FACTORISED averages.
  → "truncate-then-substitute".

For a single-decay channel these agree; for a dephasing-type channel
they differ, and master's truncate-then-substitute diagonal is the one
that stays bounded under order-2 truncation. So master does NOT "drop
the diagonal" (that was an imprecise reading in the Conclusion above) —
it keeps a DIFFERENT, stable diagonal. The current NE workaround drops
the diagonal entirely; it happens to be stable too, but only because
dropping ≈ master's bounded diagonal for these channels, and the
path-asymmetry is needed because unique_squeezing's single channel
needs the diagonal KEPT.

### Bedrock (2026-06-02): the QC-level defer-the-split is BLOCKED by eager SQA evaluation

Pushed the defer-the-split all the way to the algebra and hit a hard wall.
The plan was: build the dissipator summand with a FREE jump index `j`,
`average` + `cumulant_expansion` the generic summand, then materialize the
diagonal by substituting `j → c` for each LHS atom `c` (master's
truncate-then-substitute). This CANNOT work on top of SQA's eager model:

- SQA evaluates a commutator of two differently-named free indices on the
  same subspace as ZERO at construction (it treats `j ≠ k` by default). So
  the generic decay summand `(γ/2)(σ²¹_j[σ²²_k, σ¹²_j] + [σ²¹_j, σ²²_k]σ¹²_j)`
  evaluates to `0` for the free-`j` form. Verified directly:
  `cumulant_expansion(average(generic decay summand)) == 0`, while the
  same-atom (`j→k` first) term is correctly `−γ σ²²_k`.
- Therefore substituting `j → k` INTO the truncated generic (which is `0`)
  yields `0`, not the physical `−γ⟨σ²²_k⟩`. The diagonal only exists if the
  `j = c` substitution happens BEFORE the commutator is evaluated (so the two
  operators are recognised as same-site). That is exactly eager
  collapse-then-truncate, the unstable order.
- Inconsistently, the dephasing summand's generic form is NONZERO (its
  `σ²²_j … σ²²_j` sandwich is a product, not a vanishing commutator), and it
  carries two-`j`-factor terms. So even where the generic survives, it needs a
  sum-over-products representation.

master sidesteps both by keeping the inner moment `⟨σ_j …⟩` UNEVALUATED and
symbolic through `cumulant_expansion` (its `IndexedAverageSum`), so the
same-site product is only formed at `evaluate` time AFTER the `j → c`
substitution. SQA v0.5 deliberately removed that lazy-moment machinery (the
rewrite evaluates products eagerly). So a correct unified fix requires
RE-INTRODUCING a lazy/symbolic indexed-moment representation into SQA (the
`IndexedAverageSum` model) and threading it through QC's `cumulant_expansion`
/ `complete!` / `scale` / `evaluate` / codegen. That is a major SQA + QC
re-architecture that walks back a central rewrite simplification, not a
localized patch.

**What IS banked toward that fix:** the sum-scope canonical key
(`_canonical_key_sum`, landed dormant + unit-tested) solves the dedup-
convergence blocker that previously made Σ-wrapped-state closures
non-terminating. When the SQA lazy-moment work is done, the dedup primitive
is ready to wire in.

**What a real fix requires (defer the split):**

1. Build the dissipator with the jump index FREE (no `Σ` wrap) so SQA
   does not eager-split: keep `D[σ_j](op)` as a free-`j` product. (QC
   already forms `Jdagger·comm(op,J)` with free `j` in `_lindblad_rhs`
   BEFORE `_sum_over_jump_indices` wraps it — the split fires only at
   that wrap.)
2. `average` + `cumulant_expansion` the free-`j` summand (generic `j`,
   standard factorisation; `j` distinct from LHS atoms). Verified: this
   produces a clean truncated expression (no eager diagonal collapse).
3. Re-attach `Σ_j` and perform the diagonal/off-diagonal split at the
   AVERAGE level: `Σ_j E(j) ↦ [Σ_{j∉LHS} E(j)] + Σ_c E(j→c)` for each
   LHS atom `c` on `j`'s subspace, where `E(j→c)` is `change_index` on
   the average LEAVES (collapse allowed there, post-truncation). This
   is master's `evaluate`. The off-diagonal sum carries range `N−|LHS|`.
4. Make `scale`/`evaluate` and `find_missing` treat a `Σ`-wrapped
   average leaf as first-class state identity.

**Known failure mode (must be designed around):** the ARCHIVED attempt
below hit non-convergence — deriving a `Σ`-wrapped LHS generates `Σ`–`Σ`
cross-moment RHS leaves whose dedup keys (two bound indices) fail to
coalesce in `find_missing`, growing ~21 states/iteration past iter 6.
So step 4's multi-bound `_canonical_key`/`_dedup_key_strip_free_ne`
must canonicalise two-bound-index sum scopes, or the closure never
terminates. This is the load-bearing hard part, not the algebra.

**SOLVED in prototype (2026-06-02): the sum-scope canonical key.** The
blocker is tractable. A prototype `_canonical_key` extension coalesces
alpha-equivalent `Σ`-wrapped leaves (including the two-bound permutation
symmetry that defeated prior attempts). Verified on real off-diagonal
`Σ`-wrapped leaves: 1-bound alpha-equivalence, 2-bound alpha-equivalence,
2-bound permutation / operator-swap invariance
(`Σ_{a≠b}σ^{12}_aσ^{21}_b` ≡ `Σ_{a≠b}σ^{21}_aσ^{12}_b`), and
1-bound-vs-2-bound distinctness all hold. Recipe:

1. `SQA.assume_distinct_index` over ALL atom-space (Transition) indices
   (free + bound) so SQA order-canonicalises commuting cross-atom factors
   (the same trick §6 used for conjugate dedup, extended to bound indices).
2. Rename free indices to canonical free-slots and bound indices to
   canonical bound-slots, both minted from a per-subspace ANCHOR (a
   deterministic index independent of the input's names, e.g. the
   lex-first declared index on that subspace), in DISJOINT namespaces
   (e.g. `anchor(900+p)` free, `anchor(500+rank)` bound).
3. Brute-force over all `k!` permutations of the bound indices (k ≤ 3 for
   order ≤ 3 cumulants, so ≤ 6 perms; cheap).
4. Signature = `(sorted bound-slot names, op-term string with the Σ-scope
   ORDER stripped)`. Stripping scope order is essential: `Σ_aΣ_b = Σ_bΣ_a`,
   so `.indices` must be treated as a SET, not a vector. Take the
   lexicographic min over permutations.

Why prior attempts failed: `_canonical_key` discarded `.indices`
entirely (`return QAdd(args, Index[])`), so it could neither keep the
scope nor canonicalise its permutation symmetry. With this key, find_missing
dedup of `Σ`-wrapped leaves terminates. (Proven for the canonical-key
primitive; end-to-end closure convergence still requires building steps
1-3, but the documented load-bearing blocker is retired.)

**Validation gate:** the payoff (stability) is NOT symbolically
checkable. It requires implementing steps 1–4, then SOLVING the
heterodyne SDE at `T_end=0.1` (seed 2, EM, `dt=T_end/2e5`) and
confirming `retcode==Success`, `max|⟨a'a⟩|<1e8`, AND unique_squeezing
N=1 `⟨a'a⟩≈30.15` / N=100 plateaus stay bounded — under ONE unified
(NE-free) path. Until both pass under the unified path, the
path-asymmetric workaround stays.

**Scope/architecture note:** this reintroduces a symbolic-sum-through-
cumulant path (master's `IndexedAverageSum` model) on top of SQA's
eager-split data model — in tension with the rewrite's "use SQA's data
model natively" principle (CLAUDE.md). Worth weighing whether the fix
belongs in SQA (a lazy/deferred `Σ` that skips `_accumulate_with_diag!`
until an explicit `evaluate`) rather than QC. A dedicated session on an
isolated branch is the right vehicle; do NOT attempt incrementally on a
green tree.

**What's worth keeping anyway:** the `_prod_ops` `expand_completeness`
fold is an independent correctness improvement (kills the `get_adjoints
=true` `⟨σ^gg⟩` leak, cheap) and is the right home for §4's defensive
concern — but on its own it does NOT change the suite result (still
785/2) and does NOT resolve §3. The path-asymmetric workaround stays.

(NOTE: the heterodyne `13 vs 12` regression failure in
`examples_regression_test` was a SEPARATE, pre-existing dedup bug, NOT
the NE asymmetry: with `get_adjoints=false` the scaled closure was 13
states with NO `σ^gg`, and the fold does not change it. RESOLVED
2026-06-01 — see §6.)

## 6. Heterodyne `13 vs 12` conjugate-dedup bug: RESOLVED (2026-06-01)

Root cause: under `get_adjoints=false`, a cross-atom moment and its
conjugate should dedup to one state. The conjugate's adjoint reverses
the operator factor order (`⟨σ¹²_k σ²²_{k'}⟩` → `⟨σ²²_{k'} σ²¹_k⟩`), and
`_canonical_key`'s encounter-order rename then assigned slot names by
position, so the two keyed to *different* operator arrangements
(`σ²²·σ²¹` vs `σ²¹·σ²²`) and both survived. They are equal only after
recognising the free atom slots are mutually distinct (the framework
convention already stated in `_dedup_key_strip_free_ne`): distinct
atoms commute, so SQA can reorder the factors to a canonical
arrangement.

Fix ([src/completion.jl](src/completion.jl)): `_canonical_key` now
calls `SQA.assume_distinct_index` over the free atom-space (Transition)
slots BEFORE the encounter-order rename, via the new
`_atomspace_distinct_pairs` helper. The injected NE is stripped
downstream by `_dedup_key_strip_free_ne`, so the key stays NE-free.
This makes a moment and its order-reversed conjugate land in the same
arrangement and dedup. Heterodyne pre-scale drops 14 → 13 (the genuine
minimal `get_adjoints=false` count; the redundant conjugate is gone)
and scaled 13 → 12. The pre-scale lock in
`examples_regression_test.jl` was updated 14 → 13 with the mechanism
documented inline.

Result: full suite 792/792 (was 785/2), JET clean. The change is
restricted to the dedup key, not the stored states, and only fires on
products carrying ≥2 free atom-space indices on the same subspace, so
single-atom and Fock-only states are untouched.

### Attempt #1 (reverted): cumulant-level `expand_completeness` fold

Hypothesis: same-site Transition algebra inside a cumulant block
(e.g. `σ^{12}_i · σ^{21}_i ↦ σ^{11}_i`) produces ground-state
projectors that slip past the one-shot `expand_completeness` in
`_build_op_drift`, so the closure carries an explicit `⟨σ^{gg}⟩`
state that's completeness-redundant with `⟨σ^{ee}⟩`.

Tried: wrap `_prod_ops` in `SQA.expand_completeness(reduce(*, block))`
so the fold fires after same-site collapses. Also dropped the
LHS-to-bound NE pair in `_assume_distinct_atom_indices` and made
`_derive_for` symmetric across det/stoch.

Result:
- ✓ Closure shape changed (heterodyne pre-scale grew by 1; the new
  `⟨σ^{gg}⟩` state was confirmed via state dump).
- ✓ My narrow SDE bounded test at `T_end=0.05` passed.
- ✗ The example's actual `T_end=0.1` SDE STILL diverges. Trajectory
  is fine until t≈0.038 (well past the measurement pulse ending at
  0.02), then ramps `⟨a'a⟩` through 1e9, 1e26, 1e54, 1e135, NaN over
  5 timesteps. So the cumulant-level fold ISN'T the (full) fix; it
  shifts the bug shape but doesn't eliminate it.

Reverted. Restored the path-asymmetric NE workaround. Updated the
`heterodyne_detection (single-trajectory SDE bounded)` test to use
`T_end=0.1` (the example's value) plus an explicit
`sol.retcode == ReturnCode.Success` check so this kind of latent
divergence can't pass through a too-short test window again.

### Attempt #2: trace which symbolic term differs

Set up heterodyne under BOTH configs (path-asymmetric vs symmetric
NE-skip) at the same closure shape (12 states, no `⟨σ^{gg}⟩` extra),
diff the equations: only ONE differs, the cross-atom moment
`⟨σ_{k_1}^{22} · σ_{k_2}^{22}⟩`. Path-asymmetric RHS is 325 chars;
symmetric-skip RHS is 544. The DIFF reveals:

- Path-asymmetric has linear-decay `-γ - 2χ - pulse·η` on
  `⟨σ^{22}·σ^{22}⟩`; symmetric-skip has `-2γ - 2pulse·η`. The `-2χ`
  dephasing term and the `pulse·η` instead of `2·pulse·η` differ.
- Path-asymmetric has TWO cumulant-correction blocks:
  `2·(-γ/2 + χ + pulse·η/2)·(2⟨σ²² σ²²⟩⟨σ²²⟩ + ⟨σ²²⟩² - 2⟨σ²²⟩³)`
  and `2·pulse·η·(⟨σ¹² σ²²⟩⟨σ²¹⟩ + ⟨σ²¹ σ²²⟩⟨σ¹²⟩ + (1 - ⟨σ²²⟩)⟨σ²²⟩
   - 2⟨σ²²⟩⟨σ²¹⟩⟨σ¹²⟩)`. Both are MISSING in the symmetric-skip form.

These correction blocks ARE the "cumulant cross-decay correction
terms" from TODO §1's discussion. They emerge in the path-asymmetric
form because SQA's `_accumulate_with_diag!` splits the dissipator's
`Σ_j` over the `σ^{22}_j`-dephasing channel into diagonal (`j ∈ {k_1,
k_2}`) plus off-diagonal (`j ∉ {k_1, k_2}`), and the diagonal pieces
fold via Transition algebra into the LHS atoms' own moments. That
split is what `assume_distinct_index(..., [(k_1, j_bound), (k_2,
j_bound)])` triggers via NE propagation.

Without NE (symmetric-skip), SQA leaves the `Σ_j` and the LHS-atom
products in "Undetermined" form: free atom indices on the same
subspace with no NE annotation don't fire the diagonal split. The
cumulant truncation of an Undetermined product is then mis-shaped
(missing the cumulant cross-decay correction), and the SDE diverges
because those correction terms are the damping that bounds
`⟨σ^{22} σ^{22}⟩` against the dephasing-driven growth.

### Attempt #3: symmetric NE with `user_concretes`-immunity rules

Hypothesis: pair-skipping inside `_assume_distinct_atom_indices` can be
made symmetric across det/stoch by treating user-pinned atom indices
(those in `eqs.initial_operators`) as NE-immune. Tried two variants:

1. Skip pairs where EITHER index is user-pinned (`d in user_concretes
   && continue` plus the existing `other in user_concretes && continue`).
   Result: unique_squeezing OK; heterodyne SDE diverges (the
   `(k, j_bound)` NE pair is load-bearing for heterodyne).

2. Symmetric NE injection (drop `_derive_for`'s user_concretes-conditional
   skip, both paths inject NE via `_distinct_atom_indices`, keep current
   `other in user_concretes` filter inside `_assume_distinct_atom_indices`).
   Result: heterodyne OK; unique_squeezing diverges to 1.4e7 at N=1.

So no symmetric rule on top of the existing `_assume_distinct_atom_indices`
+ `cumulant_expansion` primitives satisfies both examples simultaneously.
The deeper asymmetry is in the dissipator structure:

- unique_squeezing: ONE atomic decay channel `J = σ(1, 2, i)`. NE on
  `(LHS, i_bound)` drops the same-site contribution which is the only
  damping channel for the LHS atom moment.
- heterodyne: dephasing-included channel set `J = [..., σ¹²_j, σ²¹_j,
  σ²²_j]`. NE on `(LHS, j_bound)` is needed to extract the correct
  cumulant-correction structure from the dephasing channel.

Both setups have a user-pinned LHS atom index against a bound dissipator
sum on the same subspace. The "right" NE behavior depends on whether
the dissipator includes a dephasing-type channel. That's too case-specific
to encode as a one-liner in `_assume_distinct_atom_indices` without
inspecting the jump operator types.

### Where the algebraically-clean fix lives: SQA (SUPERSEDED — see CORRECTED DIAGNOSIS above)

`_accumulate_with_diag!` in `src/algebra/pipelines.jl` already fires
the diagonal split unconditionally for free indices on the same
subspace as a bound sum-index: with NE annotated, the diagonal pairs
get SKIPPED (line 110, `_ne_contains(ne, sum_idx, ext_idx) && continue`)
so only off-diagonal contributes; without NE, the diagonal pairs get
added and BOTH off-diagonal and the per-diagonal substitutions are
emitted (lines 123-152). So the trigger condition isn't the bug.

What IS the bug: the no-NE case emits 1 off-diagonal + N diagonal
contributions which, taken together, should be algebraically
equivalent to the NE-annotated case's single off-diagonal-with-NE
emission (after restricting both to the same sum range). After
`cumulant_expansion` truncates, they're not equivalent. The heterodyne
diff (path-asymmetric vs symmetric-skip on the same closure shape)
shows the no-NE case loses two cumulant-correction blocks that the
NE-annotated case preserves.

So the SQA-side fix is in how `cumulant_expansion` interacts with the
per-diagonal emissions (`_accumulate_with_diag!` lines 131-152). When
SQA substitutes `sum_idx -> ext_idx` for a diagonal contribution
(line 132), it produces a term where the bound index has been
collapsed onto a free index. The subsequent
`_canonicalize!`/`_emit_scaled_by_scope!` re-emits the term. But QC's
`cumulant_expansion` (in `src/cumulant.jl`) treats each diagonal
emission as a separate symbolic average leaf, not as a piece of a
diagonal-vs-off-diagonal decomposition. So the cumulant-2 truncation
of each piece independently doesn't reproduce the truncation of the
joint sum.

Possible fixes:

1. Make `_accumulate_with_diag!` not emit explicit diagonal
   contributions for the no-NE case. Instead, just emit the
   off-diagonal-with-augmented-NE and let the user/QC choose to add
   diagonals via explicit NE annotations. This would make no-NE
   behave like NE-annotated by default, eliminating the asymmetry.
   Risk: changes the semantics of existing call sites that may
   depend on the diagonal-included emission.

2. Add an SQA primitive that returns the canonical-form
   diagonal/off-diagonal decomposition WITHOUT cumulant truncation,
   and have QC's `cumulant_expansion` recognise that primitive and
   factor it correctly (so cumulant-of-sum is computed as a single
   joint cumulant, not as cumulant-of-each-emission).

3. Document the asymmetry and keep the path-asymmetric NE workaround.

Until SQA's `_accumulate_with_diag!` or `cumulant_expansion` is updated:

- The path-asymmetric NE handling stands.
- The `_assume_distinct_atom_indices` LHS-to-bound pairing must stay
  on the stoch path.
- The relaxed assertion in
  `measurement_backaction_indices_comparison_test` stays.

The honest test target: heterodyne SDE at `T_end=0.1`, seed=2,
EM solver, `dt = T_end/2e5`, must stay bounded (`retcode == Success`,
`max|⟨a'a⟩| < 1e8`). The new SDE bounded testset locks exactly this.

### Attempt (2026-06-02, reverted): per-channel diagonal split keyed on jump type

Hypothesis (genuinely new; every prior attempt keyed on *index*
classification, this keyed on the *jump operator's type*): a DEPHASING
channel (diagonal jump `σ^{αα}`) needs the dissipator-sum diagonal
(`j = LHS atom`) dropped, a DECAY/ladder channel (`σ^{αβ}, α≠β`) needs it
kept. Implement uniformly (no det/stoch branch, no `user_concretes`) by
injecting NE between the bound jump index and the LHS atom indices, inside
`_sum_over_jump_indices` *before* the `Σ` wrap, only for diagonal jumps.

What looked promising (live-session, monkeypatched): unique_squeezing N=1
= 30.1497 and N=100 X/P = 2.29/0.44 (both correct); heterodyne single-traj
SDE bounded (max 591, close to the example's ~600); superradiant and
filter-cavity closures byte-identical to baseline (23/20, 54/48); det vs
stoch fully symmetric (21=21, 0 drift mismatch). It looked like a clean §3
resolution, so it was implemented in source and the path-asymmetric
machinery deleted.

Why it is WRONG (caught by `make test`: 797 pass, **5 errors**). All five
fail with `AssertionError: Expected an Initial parameter to exist for
variable ν` (or `χ`) at `ODEProblem` build: the dephasing rate vanished
from the compiled system entirely. Direct check: under the patch the
superradiant closure has NO `ν` anywhere. Two compounding mistakes:

1. **Pre-Σ NE deletes physical terms; post-hoc NE does not.** Injecting NE
   *before* `Σ`-construction makes `_accumulate_with_diag!` skip the
   diagonal, so the `j = k` term is never emitted. For a dephasing channel
   on a single-atom coherence that term IS the physical dephasing:
   `(ν/2)(σ²²·[σ¹²,σ²²] + [σ²²,σ¹²]·σ²²) = −(ν/2)σ¹²_k`. Dropping it deletes
   the coherence's dephasing, so `ν` disappears from the system. The OLD
   workaround's `_assume_distinct_atom_indices` runs *after* `Σ` is built,
   where (per Layer 1 above) it is a no-op on already-emitted diagonals and
   only changes NE metadata to drive the `σ^gg` fold + dedup. The two are
   not interchangeable; the pre-Σ form is lossy.

2. **Diagonal-jump-type targets the wrong channel.** For heterodyne's
   unstable moment `⟨σ²²_{k1} σ²²_{k2}⟩` the dephasing channel `σ²²`
   contributes ZERO (populations are dephasing-invariant: both diagonal and
   off-diagonal terms vanish by commutation). The load-bearing
   `(1−⟨σ²²⟩)` diagonal of Layer 2 comes from the off-diagonal PUMP channel
   `σ²¹` (`α≠β`), which `_is_diag_atom_jump` returns false for. So the rule
   never touched the channel that actually matters; heterodyne only
   "bounded" because deleting the χ dephasing happened to keep `⟨a'a⟩`
   finite. Boundedness is not correctness; the SDE bound was a false
   positive, and the drift-only physics test correctly rejected it.

Lesson: validate a closure-shape change against the *physics* regression
tests (drift-only numerics, locked master values), not just SDE
boundedness. The path-asymmetric workaround stands; defer-the-split remains
the only correct unification.

### Current source (after the 2026-06-02 structural-discriminator resolution)

- `meanfield` / both `_derive_for` methods: pass one `cross_ne =
  _has_dephasing_channel(jumps)` flag to `_build_op_drift` (identical for det
  and stoch).
- `_build_op_drift(...; cross_ne)`: when `cross_ne`, injects NE PER OP via
  `_assume_distinct_atom_indices(rhs, _distinct_atom_indices([op]))` (the op
  alone, never the `complete!` batch) → order-independent.
- `_has_dephasing_channel(jumps)` / `_is_diag_atom_jump`: the structural
  discriminator (indexed `σ^{αα}` jump present?).
- `_distinct_atom_indices(new_ops)`: unchanged logic (2+ atom-space indices on
  a subspace); `user_concretes` param removed.
- `_assume_distinct_atom_indices(q, distinct)`: `user_concretes` param removed.
- `_user_concrete_atom_indices` / `_collect_atom_indices_set!`: DELETED.
- Tests: `measurement_backaction_indices_comparison` TIGHTENED to
  `det_lhs == stoch_lhs` + every shared drift agrees; new
  `test/dephasing_discriminator_test.jl` (classification + readout-shape
  invariance).

## ARCHIVED: previous session's diagnosis (kept for audit)

**Root cause (this session's investigation)**: scale collapses
user-pinned atoms together with population slots, treating the mixed
model (one specific atom + N-1 symmetric atoms) as fully-symmetric.

Concretely, for unique_squeezing the RHS of `d⟨a⟩/dt` has
`Σ(i=1:N) ⟨σ_i₁₂⟩` (correct from meanfield: collective sum over all
N atoms). After `scale`, this becomes `N * ⟨σ_j_1_1₁₂⟩` (multiply by
range N, then canonical-rename `i` to `j_1_1`). The leaf
`⟨σ_j_1_1₁₂⟩` corresponds to the user's pinned atom 1 (via
`σ(2, 2, j(1))`). So the scaled expression says "N copies of atom 1's
dynamics" rather than the physically correct
`⟨σ_j_1⟩ + (N-1) ⟨σ_pop⟩` (atom 1 plus N-1 symmetric others).

At N=1 the bug is masked numerically (N=1·⟨σ_j_1⟩ is correct since
atom 1 is the only atom), but the closure still carries phantom
cross-atom states (`⟨σ_j_1 σ_j_1_2⟩`) that have spurious dynamics
contributing through ⟨a⟩, ⟨a'a⟩. At N=100, the all-atoms-equal-atom-1
approximation breaks completely, hence the oscillation.

**Fix surface (proposed, not yet implemented)**:

1. In SQA, the `concrete` flag distinguishes user-minted slot labels
   (`j(1)`) from algebra-internal phantoms, BUT QC's algebra-internal
   phantom mint sites all go through `(::Index)(::Integer)`, so they
   also set `concrete=true` and are indistinguishable from user
   atoms. Either:
   (a) Switch QC's phantom mints to the 4-arg `SQA.Index(name, range,
       space_index, sym)` shim (defaults `concrete=false`), turning
       `.concrete` into the discriminator.
   (b) Track user-concretes via `eqs.initial_operators` and pass them
       explicitly through scale/canonical-key.
2. `_min_slot_assignment` in [src/completion.jl](src/completion.jl)
   must NOT permute through user-concretes. Currently it does
   `slot_reps[sp][k] = first_idx(k)` for all `k`, renaming everything
   to phantoms; user-concretes (like `j_1`) should occupy fixed slots
   and only the algebra-phantom companions should permute.
3. `_scale_avg` in [src/scaling.jl](src/scaling.jl) must split sums
   against user-concretes BEFORE computing the `(range - |NE|)`
   prefactor:
     `Σ_b expr_b` becomes `expr_{b↦c} + assume_distinct_index(Σ_b
     expr_b, [(b, c)])` for each user-concrete `c` on `b`'s subspace.
   Then the prefactor on the constrained sum becomes `N - 1` and the
   `expr_{b↦c}` term emits the user-concrete contribution explicitly.

**Scope warning**: this change affects every example with a mixed
(concrete + population) shape, at least `heterodyne_detection` (uses
`σ(2,2,k(1))` concrete + `Σ_j` in H, same shape as unique_squeezing).
Heterodyne's current "12 states" lock would shift; the new closure
shape is *physically more correct* but the existing numerical-
end-state tests may need their locked values updated, with the
update justified by physics (compare master).

**Caveat for any reattempt**: the test's N=1 lock at 30.15 is a
broken-but-stable value. The example's *Effective model* (cavity-only
adiabatic-elimination Hamiltonian, exact at this order) gives
`⟨a'a⟩(t_end) = 28.038` at N=1 with the test's parameters; that is
the fix target.

**Status from this session's attempt (reverted)**:

Tried implementing steps 1-3 in src. Got the closure to split sums
correctly in `_scale_avg` (scaled `d⟨a⟩/dt` RHS shows
`⟨σ_j_1⟩ + (N-1)·⟨σ_j_1_1⟩` instead of `N·⟨σ_canon⟩`), but `mtkcompile`
then errors because the post-split RHS references `⟨σ_j_1_1⟩` (the
population representative) as a state, and that state isn't in the
pre-scale closure. To add it, `_canonical_key` must distinguish a
sum-bound leaf (`Σ_i ⟨σ_i⟩`) from an unbound leaf (`⟨σ_j_1⟩`) and
emit the sum-bound case as a separate missing state (Σ-wrapped LHS).
Tried this too: the closure no longer converges in 200 iterations
because deriving a Σ-wrapped LHS generates more Σ-wrapped RHS leaves
(at cumulant order 2: Σ-Σ cross moments). Those should still close
deterministically but the dedup keys for two-bound-index leaves
fail to coalesce, growing ~21 new states per iteration past iter 6.

The architectural problem: QC's current closure mechanism is built on
"every leaf is a single-atom moment (possibly cross-atom with
multiple specific indices)". Supporting "Σ-wrapped LHS state for the
population" requires the closure machinery to track sum-scope as a
first-class part of state identity, which spans:
- `_canonical_key` and `_dedup_key_strip_free_ne` (state identity).
- `_undo_for_derivation` (derive Σ-wrapped LHS).
- `_meanfield_deterministic` derivation of Σ-wrapped ops.
- Cumulant expansion of Σ-Σ cross moments.
- Scale handling of Σ-wrapped equations.
- `find_missing` dedup of multi-bound canonical keys.

This is a multi-file refactor; an isolated session focused on a
clean redesign is the right call, not incremental patches.

**Side-finding**: after reverting the unique_squeezing attempt to the
session's commit state (23e5700), `make test` shows heterodyne at
13/12 (one extra state) instead of the 12/12 the commit's
`make test` run logged. The standalone heterodyne reproducer
(`/tmp/het.jl`) gives 12. Likely a test-context state effect
(precompile cache, Manifest resolution, or test-ordering through
the find_missing's conjugate dedup). The earlier "776/776 passing"
state with commit 23e5700 may have depended on dep versions in a
Manifest snapshot since drifted.

## 4. Architectural: should closures auto-enforce algebraic constraints? (option 3 done; deeper options deferred)

The §3 bug class came from a completeness-redundant state
(`⟨σ^{gg}⟩` alongside `⟨σ^{ee}⟩`) leaking into the closure. The fix
calls `SQA.expand_completeness` inside `_prod_ops` so the fold fires
right after same-site collapse produces `σ^{gg}` inside a cumulant
block. This works because `expand_completeness` is the SQA primitive
that knows the N-level identity `σ^{gg} = 1 − Σ σ^{kk}`.

But the fix is point-local. It relies on every operator-product site
that can produce a completeness-redundant operator calling
`expand_completeness` afterwards. Sites today that DO call it correctly:

- `_build_op_drift` at construction.
- `_noise_drift_one` (via `expand_completeness` wrap inside
  `_build_noise_equations_forward`).

(The `_prod_ops` fold was tried — see §3 Layer 1 — but reverted: no
current config leaks `σ^{gg}` into a closure, and folding inside every
cumulant block product adds hot-loop cost for no observable benefit.
Per the recommendation below, it stays out until a leak actually
appears.)

If a future feature adds a new operator-multiplication site (e.g. a
new correlation function builder, a new `evaluate` codepath, a new
order>2 cumulant variant) and forgets the fold, the same closure-
overcounting bug reappears in that codepath. The current discipline
is "every caller remembers to fold".

### The architectural question

Should completeness-redundant states be eliminated AT THE TYPE LEVEL,
not via callers calling `expand_completeness`? Options:

1. **Make `expand_completeness` an internal invariant of operator
   products.** Lift the fold into SQA's `*(QSym, QSym)` /
   `_canonicalize!` so any same-site `σ^{gg}` that the algebra
   produces gets folded immediately. Callers stop needing to call
   `expand_completeness` explicitly; the operator algebra never
   carries an `σ^{gg}` to begin with. Cost: every product reduction
   pays the fold check; need to confirm no current code depends on
   `σ^{gg}` staying explicit (e.g. for displaying equations).

2. **Make the closure machinery enforce constraints across states.**
   Track algebraic identities (`⟨σ^{gg}⟩ = 1 − Σ ⟨σ^{kk}⟩`,
   `[a, a†] = 1`, etc.) as MTK algebraic equations alongside the
   ODE/SDE. Let the solver dispatch on DAEProblem rather than
   ODEProblem when redundant states are present. Cost: significantly
   more invasive; MTK's DAE path has different perf characteristics
   and may not work cleanly with SDE/Stochastic flows.

3. **Status quo: callers fold explicitly.** Add a quality test that
   asserts no `σ^{gg}` (and analogous completeness-redundant ops)
   appears in any post-cumulant closure across the example suite.
   Catches future leaks without changing architecture.

Recommended path forward: option 3 first (cheap defensive test).
Revisit option 1 if a similar leak shows up in a different codepath
despite the test, or if perf benchmarks show `expand_completeness` in
`_prod_ops` is a hot spot worth pushing deeper into SQA.

**DONE (2026-06-01): option 3 implemented.**
[test/completeness_invariant_test.jl](test/completeness_invariant_test.jl)
asserts no ground-state projector `σ^{gg}` (a Transition with
`i == j == ground_state`) survives into any state or drift/noise RHS
leaf, across representative scalar + indexed (2-level, Dicke, 3-level)
models, both `get_adjoints` settings, pre- and post-scale. Detection
guards on `_is_leaf_average` so a *product* of averages (whose
`undo_average` would reassemble a spurious same-site `σ^{gg}`) is not a
false positive. Currently 9/9 green on all configs, so option 1 (the
`_prod_ops` fold / lifting the fold into SQA) is deferred per the rule
above. §4 is parked here unless the test fires.

Other algebraic-redundancy classes that might exhibit similar bugs:

- Bosonic `[a, a†] = 1`. SQA's `_canonicalize!` handles this via
  normal-ordering, but if a cumulant block produces an un-normal-
  ordered product, the average could carry redundant operator
  orderings.
- Spin algebra (if added later): `S²_x + S²_y + S²_z = S(S+1)`.
- Multi-mode bosonic with cross-mode constraints.

None of these are known bugs today; they're plausible future
analogues if the architecture stays "callers fold explicitly".

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
| Per-channel pre-`Σ` NE for diagonal jumps `σ^{αα}` (unified path) | US N=1/N=100 + heterodyne SDE bounded, superradiant/filter-cavity closures identical, det/stoch symmetric; but **5 `make test` errors**: `ν`/`χ` vanish from compiled systems | Pre-`Σ` NE removes the `j=k` diagonal, which for a dephasing channel IS the physical single-atom dephasing `−(ν/2)σ¹²`; and the load-bearing diagonal for `⟨σ²²σ²²⟩` is the off-diagonal PUMP `σ²¹`, not the diagonal `σ²²`. Boundedness ≠ correctness. See the 2026-06-02 writeup above. |

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
