# QuantumCumulants.jl v1, architecture design

> Status: design proposal, awaiting implementation. Both SQA and QC are
> rewriting for 1.0 and can take breaking API changes.

## Why we are here

The current `rewrite` branch has 2,893 lines under `src/`, up from master's
~2,200, and the diff vs the last clean commit (`01ced1c`) is +1638/-369 across
9 files. The growth came from porting master's algorithms line by line on top
of SQA v0.5's new data model, then patching every place the new model
surprised the literal port. The result is workable but ugly:

- `scaling.jl` (551 lines) hand-builds `QAdd(QTermDict(QTerm(ops, ne) => c), Index[])`
  to bypass SQA's diagonal splitting, then has 12 helpers to undo the damage:
  `_renormalize_qadd`, `_stable_op_sort`, `_strip_spurious_ne`,
  `_canonicalize_op`, `_apply_index_rename`, `_distinct_indices`,
  `_index_permutations`, `_permutations_of`, `_safe_substitute`, `_qfield_key`,
  `_op_shape_key`, `_op_sort_key`.
- `evaluate.jl` (320 lines) re-implements per-operator dispatch tables
  (`_swap_op_index(::Destroy, ...)`, `(::Create, ...)`, `(::Transition, ...)`,
  etc.) that already exist inside SQA's `change_index`.
- `meanfield.jl` (295 lines) carries three near-duplicate derivation paths
  (`meanfield`, `meanfield`-noise, `meanfield`-backward). Same
  Heisenberg + Lindblad + cumulant pipeline, three times.
- `correlation.jl` (520 lines) re-derives a parallel `to_system` plus
  parameter bookkeeping rather than reusing the main bridge.

The deeper cause is that the v1 code was written scalar-first, then patched
for indices. Master had separate types
(`MeanfieldEquations`/`IndexedMeanfieldEquations`,
`scale`/`scaleME`, `evaluate`/`evalME`,
`CorrelationFunction`/`IndexedCorrelationFunction`) and separate code paths
for each pair. v1 collapsed the types but kept the parallel code paths.

This document throws away that structure and starts from "what does SQA give
us, and how do we want a user to drive a calculation."

## What SQA v0.5 actually provides

Exhaustive inventory in `docs/src/devdocs.md` of the SQA repo and the API
page. The relevant primitives, with the architectural role each plays:

| Primitive | What it does | Role for QC |
|---|---|---|
| `Index(h, name, range, space)` | symbolic summation index, user-owned | the *only* way to refer to a position in an indexed family |
| `IndexedOperator(op, i)`, `Destroy(name, space, i)` etc. | every `QSym` carries an `.index` field; `NO_INDEX` for "no index" | indexing is **always on**, never a separate type |
| `Σ(expr, i, [ne])`, `∑` | constructs sum metadata on `QAdd.indices`; eager scalarisation when expr is independent of `i` | input syntax for sums in `H`, `J`, observables |
| `change_index(expr, from, to)` | substitutes an `Index`, routing through `_canonicalize!` (or `_accumulate_with_diag!` if the QAdd still has sum scope). Handles ne propagation, diag splitting, projector squashing | **the single substitution primitive** for scale, evaluate, completion |
| `*`, `+`, `commutator`, `adjoint` | operator arithmetic that always returns a canonicalised `QAdd` | building blocks for the Heisenberg derivation |
| `simplify`, `normal_order`, `expand`, `expand_completeness`, `assume_distinct_index` | term-level rewrites | used selectively in completion / scaling |
| `average(qfield)`, `undo_average(avg)` | round-trip between operator and symbolic-scalar layers; preserves sum metadata | the leaf operation in cumulant expansion and equation derivation |
| `has_sum_metadata`, `get_sum_indices`, `get_sum_non_equal` | introspect sum scope on an average leaf | sum-collapse prefactor computation |
| `to_numeric`, `numeric_average`, `expect` | symbolic to matrix / expectation value | initial-value plumbing |
| `qadjoint`, `inner_adjoint` | distributed adjoint over BasicSymbolic trees | conjugate-pair handling in completion |
| `find_operators` (re-exported into QC) | enumerate basis of operators of a given order | used by user code, not internal pipeline |

The five operator hooks (`_site_compare`, `_can_commute`, `_commute_pair`,
`_reduce_pair`, `_ground_state_expand`) are internal extension points.
QC must not touch them.

### What is missing from SQA's public API

These are the patches we keep doing by hand in v1, each of which is a
candidate for promotion into SQA (we own SQA too):

1. **`strip_sum_scope(q::QAdd) -> QAdd`**: drop `.indices`, keep arguments.
   Currently spelled `QAdd(q.arguments, Index[])` in QC, in 5+ places.
   Public because evaluate and scale both need "treat this as a non-sum now."

2. **`set_index(op::QSym, new_idx::Index) -> QSym`**: typed replacement of
   one operator's `.index` field, all subtypes dispatched in one helper.
   QC currently re-implements this 7 times under `_swap_op_index` in
   `evaluate.jl` and inside `change_index` body in SQA. Promote SQA's
   internal dispatch and export.

3. **`fresh_index(seed::Index, label) -> Index`**: produce a new `Index`
   sharing `space_index` and `range` but with a derived name. The naming
   policy says SQA never mints these on the user's behalf; this helper is
   a "user-authored helper" that goes in QC, but takes a single canonical
   shape so QC doesn't reinvent it per file.

4. **`canonicalise_undetermined(q::QAdd; tiebreaker)`**: when two operators
   on the same Hilbert subspace have free symbolic indices with no `ne`
   constraint, SQA leaves them `Undetermined` and `_partial_sort!` does not
   reorder. For dedup of evaluate/scale results we want a purely
   lexicographic tiebreaker on (space, index-name, op-shape) so
   `σ_(b=1) σ_(a=1)` and `σ_(a=1) σ_(b=1)` collapse. This is a NEW SQA
   public function: it is NOT part of the algebra (the order matters in
   general), but it is part of the canonical-key machinery for
   `find_missing` / `evaluate` dedup. Live in SQA so QC can call it once.

5. **`enumerate_sum(q::QAdd, bound::Index, n::Integer; name_fn) -> QAdd`**:
   unroll one bound index: strip it from `.indices`, then for `k=1..n`,
   `change_index(stripped, bound, name_fn(bound, k))`, sum. Goes in SQA
   because it is purely an SQA-level transformation; QC supplies the
   `name_fn` to honour naming policy (typically "user's first-declared
   index for this space, suffixed").

6. **`assume_distinct_index` extension for "all of these are pairwise
   distinct"**: the API exists but takes an explicit list of pairs.
   A `pairwise_distinct(q, indices::Vector{Index})` convenience is useful
   in completion after `find_missing` resolves a tuple of fresh canonical
   indices.

After these additions, QC's role becomes: walk the symbolic equation tree,
call SQA primitives, build MTK systems. No hand-built `QAdd`s. No
per-operator dispatch. No re-inventing canonicalisation.

## What QC master actually does

Master's public surface, organised by what it produces:

- **Equation derivation**: `meanfield`, `indexed_meanfield`, `meanfield_backward`,
  `commutator`, `cumulant_expansion`, `cumulant`, `get_order`, `find_operators`,
  `find_missing`, `unique_ops`. Different signatures for indexed/scalar and
  forward/backward, but the core operation is `average(commutator(im*H, op) + L(op, J, rates))`
  followed by `cumulant_expansion` and `complete`.

- **Equation manipulation**: `complete!`, `indexed_complete!`, `scale`,
  `scaleME`, `evaluate`, `evalME`, `modify_equations`,
  `translate_W_to_Y`, `subst_reds`, `split_sums`. Each is a transformation
  `equations -> equations`.

- **MTK / numerics bridge**: `MTK.System(eqs)`, `initial_values`,
  `get_solution`, `get_scale_solution`.

- **Correlations**: `CorrelationFunction`, `IndexedCorrelationFunction`,
  `Spectrum`, `correlation_u0`, `correlation_p0`. Builds an extended
  Hilbert space + ancilla operators, computes a derived equation set,
  exposes the Laplace-transform spectrum.

- **Measurement / retrodiction**: `MeanfieldNoiseEquations`,
  `BackwardMeanfieldNoiseEquations`, `meanfield_backward`,
  `translate_W_to_Y`, `_master_noise_dY`. Adds a stochastic-noise field
  to the equations.

- **Indexing helpers** (mostly internal): `insert_index`, `insert_indices`,
  `change_index`, `order_by_index`, `reorder`, `inorder!`, `_inconj`,
  `inadjoint`, `ismergeable`, `create_index_arrays`, `create_value_map`,
  `eval_term`, `evalEq`, `count_eq_number`, `subst_reds_eval`. These
  live in master because master had a parallel indexed type system; v1
  obsoletes most of them by collapsing the type system.

### What we keep, drop, merge

- **Keep, in unified form**: `meanfield`, `complete!`, `scale`, `evaluate`,
  `cumulant_expansion`, `find_missing`, `find_operators`,
  `CorrelationFunction`, `Spectrum`, `to_system`, `initial_values`,
  `get_solution`, `modify_equations`, `translate_W_to_Y`,
  `meanfield_backward`. Each one function, dispatching uniformly on
  `MeanFieldEquations` regardless of indexing.

- **Drop entirely**: `indexed_meanfield`, `indexed_complete!`, `scaleME`,
  `evalME`, `IndexedCorrelationFunction`, `IndexedMeanfieldEquations`,
  `EvaledMeanfieldEquations`, `ScaledMeanfieldEquations`,
  `IndexedMeanfieldNoiseEquations`, `BackwardMeanfieldNoiseEquations`,
  `insert_index`, `insert_indices`, `eval_term`, `evalEq`, `evalME`,
  `subst_reds_eval`, `count_eq_number`, `check_arr`, `order_by_index`,
  `reorder`, `inorder!`, `_inconj`, `inadjoint`, `ismergeable`, `getAvrgs`,
  `containsIndexedOps`, `writeNeqs`, `ClusterSpace`, `cluster_expand`,
  `plotME`. Most of these are duplicates of an existing capability or
  artefacts of the dual indexed/scalar type system.

- **Merge**: `meanfield` and `meanfield_backward` into one function taking
  `direction::Forward()` or `Backward()`. `complete!` and
  `indexed_complete!` already merged in v1; same for `scale`.

## The new architecture

### One equation type

```julia
struct MeanFieldEquations{O, Op, OpEq, J, Jd, Iv, S, D <: EvolutionDirection}
    equations::Vector{Symbolics.Equation}      # d⟨op⟩/dt ~ rhs
    operator_equations::Vector{Symbolics.Equation}  # dop/dt ~ rhs (for derivation reuse)
    states::Vector{S}                          # the ⟨op⟩ LHS, in order
    operators::Vector{Op}                      # the op corresponding to each state
    hamiltonian::Op
    jumps::Vector{J}
    jumps_dagger::Vector{Jd}
    rates::Vector
    iv::Iv                                     # MTK independent variable (t)
    order::O                                   # cumulant truncation, Int or Vector{Int}
    direction::D                               # Forward() or Backward()

    # Optional noise channel
    noise_equations::Vector{Symbolics.Equation}    # empty Vector if no noise
    operator_noise_equations::Vector{Symbolics.Equation}
    efficiencies::Vector                           # empty if no noise

    # Canonicalisation alphabet (computed once, immutable)
    canon::Dict{Int, Vector{SQA.Index}}        # space_index to user's declared indices

    # Provenance: which transforms have been applied (for dispatch and
    # for catching incompatible composition like `scale(scale(eqs))`).
    history::Vector{Symbol}                    # :derived, :completed, :scaled, :evaluated
end
```

Rationale:
- **Indexing is always present.** An "unindexed" equation is just one whose
  operators all carry `NO_INDEX`. No separate type.
- **Noise is a presence/absence flag.** `isempty(efficiencies)` is the
  predicate. No separate `NoiseMeanFieldEquations` type. (Reconsidered:
  keep two types if the type parameterisation pays performance, or collapse
  if not. Decide via prototype.)
- **Direction is a type parameter** so `Forward` and `Backward` dispatch
  on `_lindblad_term`, and the user-facing constructors stay simple.
- **`canon` is computed once** at construction time from
  `get_indices(hamiltonian) ∪ jumps ∪ operators` and is invariant through
  every subsequent transform. No more recomputing `_build_canonical_indices`
  inside `find_missing`, `scale`, `evaluate`.
- **`history`** lets `scale` refuse to operate on `:evaluated` equations
  (concrete positions, no symmetry to collapse), and lets `to_system`
  validate that `evaluate` baked in all symbolic ranges.

### The five transformations

Each is one function. Each operates on `MeanFieldEquations`. Each calls SQA
primitives only.

```julia
# 1. Derive equations from physical input.
derive(ops, H, J=[]; rates=ones(length(J)), efficiencies=Float64[],
       order=nothing, direction=Forward(), simplify=true,
       mix_choice=maximum, iv=_default_iv()) -> MeanFieldEquations

# 2. Close the system: derive equations for any RHS-only averages.
complete!(eqs; max_iter=200, simplify=true, filter_func=nothing,
          mix_choice=maximum, get_adjoints=true) -> MeanFieldEquations

# 3. Symmetry collapse (permutation symmetry across same-Hilbert-subspace atoms).
scale!(eqs; h=nothing) -> MeanFieldEquations           # h: per-subspace selector

# 4. Materialise symbolic ranges to integers.
evaluate(eqs; limits, h=nothing) -> MeanFieldEquations

# 5. Emit an MTK System ready for ODEProblem/SDEProblem.
to_system(eqs; name) -> MTK.System
```

Plus three glue functions:

```julia
initial_values(eqs) -> Dict{Num, ComplexF64}        # defaults to 0
initial_values(eqs, state) -> Vector{ComplexF64}    # state::Ket via SQA.numeric_average
get_solution(sol, op, eqs) -> (τ -> Number)         # reverse lookup
```

And five user-driven utilities (each one function, dispatching on
`MeanFieldEquations`):

```julia
find_missing(eqs; get_adjoints=true) -> Vector{<:BasicSymbolic}
modify_equations!(eqs, f) -> MeanFieldEquations
translate_W_to_Y(eqs) -> MeanFieldEquations          # noise-only
cumulant_expansion(expr, order; ...)               # tree-walker helper
get_order(expr) -> Int
```

That is the entire public surface. Roughly 13 functions vs master's 30+.

### How indexing flows through every stage

Indexing is treated identically at every stage because SQA does it for us.
The flow:

1. `derive(ops, H, J)`: user constructs `H` using `Σ(...)`. SQA puts the sum
   scope into `QAdd.indices`. Heisenberg derivation calls `commutator(H, op)`,
   which is SQA's `*` and `+`. These propagate scope correctly.
   `cumulant_expansion` walks the result and factorises averages. The result
   is `MeanFieldEquations` whose equations may carry sum-metadata on average
   leaves.

2. `complete!(eqs)`: walks RHS, finds leaf averages not present in
   `eqs.states`. Dedup is a single `state_key(avg, canon)` function: undo
   the average, canonicalise free indices to the canon-first-declared index
   for each space, then return the resulting QAdd as the key. Conjugate
   pair handling: same key after `qadjoint`.

3. `scale!(eqs)`: walks every equation, alpha-renames every free `Index`
   in every operator to the canon-first index for that space (a single
   `change_index` call per index). Sum scope `.indices` collapses to a
   prefactor `(range - count_ne)` *only when* the bound index appears in
   the operator. Otherwise the scope is spurious; see "diagonal splitting"
   in SQA's devdocs. Dedup by `state_key`.

4. `evaluate(eqs; limits)`: walks every equation, finds every `Index`
   anywhere whose `range` is a key of `limits`. For LHS-free indices,
   enumerates the cartesian product `1..limits[range]` and emits one
   equation per assignment. For RHS-only indices (implicitly summed,
   inherited from the original `Σ` that cumulant expansion may have
   broken up), wraps the leaf in `enumerate_sum`. Dedup by `state_key`.

5. `to_system(eqs)`: builds a `Dict{avg, MTK.Num}` map. For each equation,
   `Differential(t)(map[lhs]) ~ rewrite_rhs(rhs, map)`. Parameter discovery
   walks the rhs and collects every BasicSymbolic that is not in `map` and
   not the independent variable. Returns `MTK.System`. Noise channel adds
   `MTK.@brownians _qc_dW` and noise-RHS expressions.

The dedup function `state_key(avg, canon)` is the single source of truth
for "are these two states the same physical observable." It is defined once
in `completion.jl` and used by `find_missing`, `scale!`, `evaluate`, and
`to_system`'s reverse lookup.

### How CorrelationFunction stops being bespoke

Master's `CorrelationFunction` builds an extended Hilbert space, constructs
new operators on it, calls `meanfield` to derive the τ-equation, then plumbs
the ambient steady-state via `correlation_u0` / `correlation_p0`. That logic
gets re-implemented inside `to_system(CorrelationFunction)` in v1, which is
why correlation.jl is 520 lines.

The cleaner design:

```julia
function CorrelationFunction(op1, op2, eqs::MeanFieldEquations;
                             steady_state=false)
    # extended space and ancilla operator construction: lives here
    h_ext = SQA.extend_hilbert(eqs.hamiltonian, op2)
    op2_anc = SQA.create_ancilla(op2, h_ext)

    # Re-derive on extended space: this is just `derive` again
    eqs_τ = derive([op1 * op2_anc], lift(eqs.hamiltonian, h_ext), ...)
    complete!(eqs_τ)

    # Ambient averages become parameters at to_system time.
    # We just record which states in eqs_τ correspond to ambient
    # (non-τ-evolving) averages from eqs.
    ambient = identify_ambient_states(eqs_τ, eqs)

    return CorrelationFunction(op1, op2_anc, eqs_τ, eqs, ambient, steady_state)
end
```

Then `to_system(::CorrelationFunction)` is `to_system(c.eqs_τ)` plus a
parameter-substitution layer that maps `ambient` averages to the user's
chosen steady-state values. Both `correlation_u0` and `correlation_p0`
disappear into a single `setup_correlation_problem(c, ss_sol)` that returns
`(u0_dict, p_dict)` ready for `ODEProblem`.

`Spectrum` is unchanged in spirit but operates on the unified `eqs_τ`. It
is Laplace-transform on a generic `MeanFieldEquations`, not a special method
for `CorrelationFunction`.

### MTK boundary

Exactly one file, `mtk.jl`, owns the MTK interface. Three functions:

```julia
to_system(eqs::MeanFieldEquations; name) -> MTK.System
initial_values(eqs, args...) -> Dict or Vector
get_solution(sol, op, eqs) -> (τ -> Number)
```

Parameter discovery is **explicit**, not tree-walk. When `derive` builds the
equation set, it also records the set of free symbolic parameters
(`Set{BasicSymbolic}`) appearing in `H`, `J`, `rates`, `efficiencies`. Every
subsequent transform (`scale`, `evaluate`, `modify_equations`) updates this
set. `to_system` reads it directly. This eliminates the entire class of
"Expected an Initial parameter to exist for variable N3" bugs: post-
evaluate, the `N3` was removed from the parameter set, so MTK never sees it.

The variable map `Dict{Average, MTK.Num}` is built once in `to_system` and
returned alongside the System for `get_solution`'s reverse lookup, instead
of re-derived from the equation list each call.

### Provenance and composition rules

`history::Vector{Symbol}` constrains composition:

| Operation | Allowed history | Resulting history |
|---|---|---|
| `derive(...)` | (any) | `[:derived]` |
| `complete!` | ends in `:derived`, `:scaled`, `:evaluated` | append `:completed` |
| `scale!` | `:completed` last, `:evaluated` absent | append `:scaled` |
| `evaluate` | `:completed` last | append `:evaluated` |
| `to_system` | `:completed` or `:evaluated` last | (MTK.System produced) |
| `modify_equations!` | any | append `:modified` |
| `translate_W_to_Y` | noise-mode, `:completed` last | append `:translated` |

Violations raise an early, informative error. The current spaghetti happened
in part because no rule said "you cannot scale an evaluated system", so the
code accumulated patches to handle that case anyway.

## Concrete file layout

```
src/
  QuantumCumulants.jl    (45 lines)   exports + includes
  equations.jl           (~150)        MeanFieldEquations + history + canon
  cumulant.jl            (~200)        cumulant_expansion, cumulant, get_order
  derive.jl              (~150)        derive() (was meanfield.jl)
  complete.jl            (~150)        complete!, find_missing, state_key
  scale.jl               (~100)        scale!
  evaluate.jl            (~80)         evaluate
  correlation.jl         (~250)        CorrelationFunction, Spectrum
  mtk.jl                 (~200)        to_system, initial_values, get_solution
  measurement.jl         (~80)         translate_W_to_Y, modify_equations
  noise.jl               (~100)        noise drift terms used by derive.jl
  latexify.jl            (~30)         display recipes
```

Total target: ~1,500 lines (vs current 2,893). Roughly half the code,
with the same feature set.

Key file responsibilities:

- **equations.jl**: the struct, constructors, `_copy`, `_append!`, history
  validation. No algorithms.
- **cumulant.jl**: the Wick factorisation. Tree walker, no SQA-specific
  branching.
- **derive.jl**: ONE `derive(ops, H, J; ...)` function. Dispatches on
  `direction` and `isempty(efficiencies)` for the four variants.
- **complete.jl**: `state_key` is here. `find_missing` uses it. `complete!`
  uses `find_missing` + `derive` on the missing ops.
- **scale.jl**: alpha-rename free indices to canon; sum-scope prefactor
  collapse. One pass per equation.
- **evaluate.jl**: enumerate cartesian product of LHS-free index ranges;
  unroll RHS-implicit summed indices via `enumerate_sum`. One pass.
- **correlation.jl**: build extended-space equations via `derive`; expose
  `Spectrum` as a Laplace operation on the resulting `MeanFieldEquations`.
- **mtk.jl**: bridges everything. Knows nothing about indexing. By the time
  it sees an equation set, all indices are concrete (post-evaluate) or
  sym-canonical (post-scale).
- **measurement.jl**: post-derivation transforms specific to the
  noise/retrodiction workflow.

## SQA changes we'll request

Promote-to-public in SQA, with breaking-change tolerance:

1. `strip_sum_scope(q::QAdd) -> QAdd`: new export.
2. `set_index(op::QSym, idx::Index) -> QSym`: new export, single dispatch
   table.
3. `canonicalise_undetermined(q::QAdd; tiebreaker=default_lex) -> QAdd`: new
   export. Used for state-key dedup in QC.
4. `enumerate_sum(q::QAdd, bound::Index, n::Integer; name_fn) -> QAdd`:
   new export. The naming policy is preserved by user-supplied `name_fn`.
5. Make `_substitute_ne`, `_canonicalize!` officially private and not
   re-exported even by accident.

Plus a few smaller items:
- A version of `change_index` that accepts a `Vector{Pair{Index,Index}}`
  for batch substitution. Avoids the N-call loop in QC.
- `pairwise_distinct(q::QAdd, idxs::Vector{Index})` convenience.
- Public access to the `_partial_sort!` tiebreaker so QC can compose with
  its own dedup key without duplicating SQA's order rules.

## Migration plan from current `rewrite` branch

This is not a rewrite-from-scratch. We keep all the *correct* code and
delete only what's redundant or fighting SQA. Ordered by blast radius:

### Phase 1: Refactor scaling.jl (3 commits, ~1 day)
1. Delete `_renormalize_qadd`, `_stable_op_sort`, `_strip_spurious_ne`,
   `_canonicalize_op`, `_apply_index_rename`, `_distinct_indices`,
   `_index_permutations`, `_permutations_of`, `_qfield_key`,
   `_op_shape_key`, `_op_sort_key`, `_safe_substitute`.
2. Replace with a single pass: for each equation, for each free index `i`
   in the LHS operator, `change_index(eq, i, canon[i.space_index][1])`.
   For each average leaf in the RHS with non-empty `.indices`, prefactor
   `(i.range - count_ne)` only when `i` appears in `op`. Dedup by
   `state_key`.
3. Run `make test`. Expect scaling regressions; fix the dedup-vs-canon
   interaction.

Success criterion: `scaling.jl` drops to ≤120 lines and all of
`test/scaling_test.jl` plus `test/indexed_meanfield_test.jl` pass.

### Phase 2: Rewrite evaluate.jl on the new primitives (1 day)
1. Same shape as Phase 1's scale but enumerating over `limits[range]`
   integer values.
2. Distinguish LHS-free (enumerate, produce N equations) vs RHS-only-free
   (implicitly summed, wrap in `enumerate_sum`).
3. Restore the new `test/indexed_evaluate_test.jl` as a passing test.

Success criterion: `evaluate.jl` ≤100 lines; the three pending tests
(`indexed_scale_test`, `indexed_mixed_order_test`,
`indexed_filter_cavity_test`) port over and pass with equation-count
assertions matching master.

### Phase 3: Unify meanfield.jl into derive.jl (1 day)
1. Collapse three near-duplicate paths into one parameterised function.
2. Move the noise-drift helpers to `noise.jl` as composable terms.
3. `direction::Forward()/Backward()` controls sign; `efficiencies=[]`
   vs non-empty controls noise.

Success criterion: `derive.jl` (new) ≤200 lines; all of
`test/meanfield_test.jl` and `test/noise_test.jl` pass.

### Phase 4: Refactor correlation.jl (1 to 2 days)
1. Express `CorrelationFunction` construction via `derive` on extended
   Hilbert space.
2. Move ambient-parameter plumbing into a single `setup_correlation_problem`
   helper.
3. `to_system(::CorrelationFunction)` delegates to `to_system(c.eqs_τ)`.

Success criterion: `correlation.jl` drops from 520 to ≤250 lines; all of
`test/systems/spectrum_test.jl` passes.

### Phase 5: SQA additions (parallel, 0.5 day)
1. Add `strip_sum_scope`, `set_index`, `canonicalise_undetermined`,
   `enumerate_sum`, batch `change_index`, `pairwise_distinct` to SQA.
2. Update QC to use the new public API in place of the local helpers.

Success criterion: every `_apply_*` / `_swap_*` / manual `QAdd(...)` call
in QC `src/` is gone.

### Phase 6: Port master tests (1 day)
1. SQA-level tests in `test/pending/` (operator algebra, fock, nlevel,
   spin, average_sums, double_sums, index_basic, parameters) move to
   SQA's test suite.
2. QC-surface tests in `test/pending/` (numeric_conversion,
   measurement_*, higher_order, indexed_*, scaling, spin) port to
   v1 surface and unskip.

Success criterion: `test/pending/` is empty.

## Risks and open questions

1. **Numerical agreement with master.** The new evaluate / scale may produce
   equations with different surface form (operator ordering, NE constraints)
   than master's. We assert agreement *only* at the ODE solution level
   (steady states, transients up to ε), not at the symbolic level.
   Master's exact equation counts are NOT a contract.

2. **`Forward`/`Backward` parameterisation.** Currently a type parameter on
   `NoiseMeanFieldEquations`. Keeping that as a type parameter on the unified
   `MeanFieldEquations` adds compilation cost for users who don't need
   backward. Alternative: a runtime field. Decide via prototype.

3. **`MeanFieldEquations` vs `NoiseMeanFieldEquations` collapse.** Putting
   `noise_equations` always-present (empty when absent) means every
   `MeanFieldEquations` carries the extra fields. Probably fine; benchmark
   `_copy` and serialisation overhead.

4. **Parameter discovery: explicit vs tree-walk.** Explicit is robust but
   requires every transform to maintain the set. Tree-walk is what we have
   and is buggy in edge cases. Recommend explicit + a `validate_parameters`
   call in `to_system` that tree-walks and asserts equality, catching
   regressions.

5. **What if SQA's `change_index` is the bottleneck?** Current `change_index`
   re-canonicalises every call. Batch substitution + `set_index` may give
   real speedups on indexed completion of 1000+ equation systems. Profile
   before optimising.

## What I want feedback on

- **Direction parameterisation**: type parameter vs runtime field.
- **`NoiseMeanFieldEquations` collapse**: yes/no.
- **SQA additions**: any of the six I proposed feel wrong-shaped?
- **Provenance enforcement**: too rigid? Should `evaluate(scale(eqs))`
  be allowed?
- **MTK parameter discovery**: explicit vs tree-walk?
- **Test porting scope**: which master tests are user-facing contracts
  vs which were master-internal unit tests we don't owe?

Once these are pinned down, the migration plan executes one phase at a
time with green tests between each phase.
