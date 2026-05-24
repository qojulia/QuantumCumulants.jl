# TODO

Open work items for the v1 rewrite. Each entry names the failure mode and where it
shows up so the fix can be verified end to end.

## Numerics / MTK bridge

### Cluster correlations lost: scale collapses distinct atoms to one slot

In the superradiant laser, the expected RHS of `d/dt ⟨a' σ_1^{12}⟩` includes
a two-atom cluster term `+ i g_1 (N-1) ⟨σ_1^{12} σ_2^{21}⟩`. After
`scale`, this term currently shows up as `(N-1) ⟨σ_i^{22}⟩` because:

1. `scale` aliases every cluster-space free `Index` to canonical-first
   ([src/scaling.jl](src/scaling.jl) `_scale_qadd`). For 2-atom correlations
   `⟨σ_i^{12} σ_j^{21}⟩` with `i ≠ j`, this renames `j → i`.
2. The `(i, j)` NE entry on the renamed term becomes `(i, i)`. SQA's
   [`_substitute_ne`](../SecondQuantizedAlgebra.jl/src/expressions/qterm.jl)
   silently drops contradictory pairs instead of zeroing the term.
3. SQA's same-site reduction then fires: `σ_i^{12} σ_i^{21} → σ_i^{22}`.

Net result: a 2-op cluster average is mis-folded into a single-atom
average. The system still closes and runs but the dynamics is wrong for
multi-atom correlations.

A separate but related problem: when `complete` derives an equation for an
observable like `⟨a' σ_i^{12}⟩` whose free index `i` collides with `H`/`J`'s
sum index, SQA's `commutator(im*H, op)` loses the `Σ_{k≠i}` off-diagonal
contribution (the cluster term never reaches `scale` in the first place).

**Sketch of the fix (scope is significant, ~half-day):**

- *Scale side:* replace the canonical-first collapse with per-Hilbert-space
  graph coloring: free indices that any term marks NE-distinct get
  separate canonical "slot" reps (`i`, `i_2`, `i_3`, ...) minted from the
  user's first-declared index. Two-atom cluster averages survive as
  `⟨σ_i^{12} σ_{i_2}^{21}⟩`.
- *Completion side:* alpha-rename `H`/`J`'s bound sum indices to fresh
  names before calling `meanfield` in `_derive_for`, so the commutator
  diagonal split keeps the `Σ_{k≠i}` part.
- *Dedup side:* `_canonical_key` (used by `find_missing`) and
  `_scale_state_key` need a normalization step that uses the
  "different-slot-implies-distinct" convention to fold operator orderings
  and strip irrelevant NE bookkeeping. Without this, the new cluster
  states multiply via alpha-equivalence variants.
- *Test updates:* `indexed_filter_cavity: per-Hilbert-space
  evaluate/scale commute` expects the *old* collapsed equation count;
  update to the new count once the design is in.

Investigation log on the `rewrite` branch (see git log for the σᵢᵢ leak fix)
confirmed the diagnosis; a from-scratch implementation pass is what's
missing. Until then, mean-field results that rely on two-atom cluster
correlations under scale are quantitatively off.

Surfaces in: `examples/superradiant_laser_indexed.jl` (cluster term is
silently folded into a single-atom average; trajectory shape changes vs.
master).

### Solver-accuracy follow-up: `conj(cosh(ξ))` in symbolic RHS

`examples/unique_squeezing.jl`: full vs. effective-model squeezing should
agree at large `N`, but the RHS carries unresolved `cosh(ξ)·conj(cosh(ξ))`
/ `cosh(ξ)*` (conjugate) factors even though `ξ` is real-typed. Check
whether MTK simplifies `conj(cosh(ξ)) → cosh(ξ)` for `ξ::Real` declared
via `@variables ξ`, and if not, either pre-simplify in QC's `System`
builder or document the workaround.

Surfaces in: `examples/unique_squeezing.jl`.

## Examples status

Examples still blocked on a numerical / plotting step:

- `examples/superradiant_laser_indexed.jl`: trajectory off-master because
  `scale` folds the two-atom cluster correlation `⟨σ_1^{12} σ_2^{21}⟩` into
  a single-atom average (see "Cluster correlations lost" above). The
  σᵢᵢ·σⱼⱼ idempotency leak that previously affected this example is
  fixed (per-jump dissipator now correctly summed over the jump index).
- `examples/waveguide.jl`: MTK v10 rejects `Ω_i_i` (the diagonal the
  Hamiltonian never references). Wrap the param dict in
  `parameter_map(sys, …)` to drop the diagonal entries.
- `examples/retrodiction_homodyne.jl`: still commented out of the docs build
  (`docs/make.jl`); depends on the `Backward()` SDE path and `modify_equations`,
  which have not been ported.
