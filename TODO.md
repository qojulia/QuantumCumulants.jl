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
- `heterodyne_detection.jl` (deterministic pulse only; SDE part regressed,
  see below)
- `cavity_antiresonance_indexed.jl` (antiresonance dip restored; the
  literal-key fix in `_canonicalise_avg_leaves` resolved the apparent
  dipole-dipole bug, which was actually an alpha-rename collapse of
  distinct per-atom states)
- `filter-cavity_indexed.jl` (filter populations stable at ~0.01,
  intensity spectrum is a clean Lorentzian; after canon-slot lookup
  update in the example file)

### Open regressions

- `heterodyne_detection.jl`: deterministic pulse (plot 1) matches master,
  every SDE trajectory (plots 2, 3, 4) blows up to ~1e220 within ~0.02 ms.

  Root cause: state-space size mismatch from differing
  Hermitian-conjugate-pair handling in `find_missing`.

  - Master on its native stack (Julia 1.11 + QC v0.4.3 + ModelingToolkit
    v10.32.1 + StochasticDiffEq v6.87.0): `scaled_eqs` has **12
    unknowns**, all 200 trajectories complete, `max n_avg = 442.4`
    photons. Confirmed by `cd /var/tmp/qc-master && julia +1.11
    --project=examples examples/heterodyne_detection.jl`.
  - Rewrite (Julia 1.12 + QC v0.5.0 + MTKBase 1.36.2 + StochasticDiffEq
    v7.0.0): same operator list compiles to **18 unknowns**. The 6
    extras are conjugates of states that are already tracked, plus
    `⟨σ_11⟩`:

  | rewrite extra | partner already tracked |
  | ------------- | ----------------------- |
  | `⟨a*σ_12⟩`    | `⟨a†*σ_21⟩` (conj)      |
  | `⟨a†*σ_21⟩`   | `⟨a*σ_12⟩`              |
  | `⟨σ_12*σ_12⟩` | `⟨σ_21*σ_21⟩` (conj)    |
  | `⟨σ_21*σ_21⟩` | `⟨σ_12*σ_12⟩`           |
  | `⟨σ_22*σ_22⟩` | (Hermitian, redundant)  |
  | `⟨σ_11⟩`      | `1 - ⟨σ_22⟩`            |

  Master's `find_missing` (`src/utils.jl`) filters conjugates out of
  the missing-state pool when `get_adjoints=true`. Rewrite's
  `find_missing` (`src/completion.jl:56-57`) does the opposite: with
  the default `get_adjoints=true` it explicitly pushes BOTH the state
  and its conjugate into `missing_states`, doubling the
  conjugate-active dimension.

  Why this is an SDE blowup and not a drift bug. Drift preserves the
  conjugate-pair invariant `u_i(t) = conj(u_j(t))` exactly (drifts for
  `i` and `j` are complex conjugates by construction; plot 1 confirms
  the ODE is fine). The diffusion column for `i` and `j` is also a
  symbolic conjugate pair, but at runtime the compiled `g(u, p, t)`
  evaluates them on two independent floating-point paths. Both terms
  multiply the same `dW`, so any roundoff between `g[i]` and
  `conj(g[j])` drifts `u[i] - conj(u[j])` away from zero each step.
  Atom-cavity coupling `Ng = 2.3e5` amplifies any `1e-12` slip within
  about 100 steps; once the invariant is broken, drift no longer
  damps and `⟨a*a⟩` / `⟨σ_12*σ_22⟩` grow exponentially in opposite
  directions. Master never exhibits the drift because there is no
  `u[j]` to drift against, only `conj(u[i])`. The instability is
  laser-regime-specific (stable for `N <= 5000`).

  Hypotheses confirmed-not-the-cause along the way
  (`/tmp/diag_het*.jl`):
  noise-formula reordering, dropping `expand_completeness`, adding
  `Symbolics.simplify` before/after cumulant truncation,
  `brownians`-vs-`noise_eqs` construction path, RNG seed, `dt`
  refinement, solver choice (`EM`, `EulerHeun`, `SROCK1`, `SRA1`,
  `SRA2`, `SOSRI`, `SOSRI2`, `SRIW1`, `ImplicitEM`,
  `ImplicitEulerHeun`, `SKenCarp`), MTKBase version (1.36.2 vs 1.41.0
  identical), `exp(0)` / `sin² + cos²` symbolic leftovers,
  time-dependent jump (`J = a*exp(iωlt)` vs `J = a`), strict Itô vs
  Stratonovich. Minimal complex SDE with conj() compiles and runs
  cleanly under the same stack, so the MTKBase pipeline itself works
  in isolation.

  Concrete fix path:
  1. Flip the default in `src/completion.jl::find_missing` to
     `get_adjoints=false` so completion produces only one of each
     conjugate pair (matches master).
  2. Tighten `_conj_substitution_dict` in `src/mtk.jl` so it keys on
     the same canonical form the noise-rhs leaves use (walk via
     `_canonical_key` instead of raw `_avg_conj_for_codegen`).
     Without (2), `complete(eqs; get_adjoints=false)` produces 13
     unknowns but `mtkcompile` then refuses with `Brownian _qc_dW
     appears non-linearly` because literal `⟨σ_21 σ_22⟩` averages on
     the noise rhs are not rewritten to
     `conj(var_for_⟨σ_12 σ_22⟩)`.
  3. Audit other examples and tests for assumptions that both members
     of every conjugate pair are tracked as separate states (e.g.
     `get_solution` lookups, initial-condition mapping). Anywhere that
     reads `eqs.states[i]` expecting the conjugate to also be there
     needs to go through `conj(...)` instead.

  Documenting unfixed in this round; plot 1 stays correct, plots 2-4
  remain a known regression until the above lands.

### Minor numerical divergence

- `excitation-transport-chain.jl`: dashed "end of chain" trace settles at
  ~0.10 on rewrite vs ~0.11 on master. Smaller than pre-fix gap and well
  within the noise of a JumpProblem ensemble; could be a remaining
  algebraic difference. Low priority.
- `retrodiction_homodyne.jl`: docs build entry in `docs/make.jl` is still
  commented out to match master (PR #266 history; suspected 5 x 10^4 step
  SDE solve being too heavy for docs).

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
