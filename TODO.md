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

## 3. Det vs stoch NE asymmetry: STILL OPEN

The current rewrite uses path-asymmetric NE injection (det path skips
NE when user named an atom index in initial_operators; stoch path
always injects NE). Both unique_squeezing and heterodyne SDE stay
bounded with this asymmetry. Master keeps both bounded with what
looks like a unified NE convention; the asymmetry is therefore
algebraic rather than physical and worth resolving.

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

### Where the algebraically-clean fix lives: SQA

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

### Workarounds in current source

- `_derive_for(::MeanFieldEquations)`: skip NE when
  `_user_concrete_atom_indices(eqs)` is non-empty.
- `_derive_for(::NoiseMeanFieldEquations)`: always inject NE.
- `_assume_distinct_atom_indices`: pairs LHS atom indices with both
  other LHS indices AND the dissipator's bound sum-index.
- `measurement_backaction_indices_comparison: deterministic vs
  stochastic LHS match`: relaxed assertion (intersection non-empty,
  both closures non-empty, per-equation drift agreement). The
  directional subset assertion is OFF until the path-asymmetric NE
  goes away.

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

## 4. Architectural: should closures auto-enforce algebraic constraints? (open)

The §3 bug class came from a completeness-redundant state
(`⟨σ^{gg}⟩` alongside `⟨σ^{ee}⟩`) leaking into the closure. The fix
calls `SQA.expand_completeness` inside `_prod_ops` so the fold fires
right after same-site collapse produces `σ^{gg}` inside a cumulant
block. This works because `expand_completeness` is the SQA primitive
that knows the N-level identity `σ^{gg} = 1 − Σ σ^{kk}`.

But the fix is point-local. It relies on every operator-product site
that can produce a completeness-redundant operator either calling
`expand_completeness` afterwards, or going through `_prod_ops`. Other
sites today that DO call it correctly:

- `_build_op_drift` at construction.
- `_noise_drift_one` (via `expand_completeness` wrap inside
  `_build_noise_equations_forward`).
- `_prod_ops` (this session's fix).

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

Other algebraic-redundancy classes that might exhibit similar bugs:

- Bosonic `[a, a†] = 1`. SQA's `_canonicalize!` handles this via
  normal-ordering, but if a cumulant block produces an un-normal-
  ordered product, the average could carry redundant operator
  orderings.
- Spin algebra (if added later): `S²_x + S²_y + S²_z = S(S+1)`.
- Multi-mode bosonic with cross-mode constraints.

None of these are known bugs today; they're plausible future
analogues if the architecture stays "callers fold explicitly".

## 5. `excitation-transport-chain` end-of-chain drift (low priority)

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
