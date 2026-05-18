# v1 regressions discovered during test porting

These were identified while porting the master test suite. Each item is a
behavior the master branch handled but the v1 rewrite either drops, weakens,
or gets wrong. Listed in rough priority order.

## `_filter_rhs!` originally tried to reapply `complex(...)`

**Symptom:** when `complete!(...; filter_func=phase_invariant)` recursed
into a `complex(re, im) * ⟨...⟩` subtree, it called `complex(::SymReal,
::SymReal)` to rebuild the node — no method.

**Fixed:** `_filter_expr` now short-circuits via `_has_average(x)`. Calls
that contain no average leaves are returned unchanged.

## `find_missing` flagged products of averages as missing states

**Symptom:** master's `_collect_missing!` walked the average leaves only;
v1's used `is_average(x)` which is true for any AvgSym-typed subtree (so
`⟨a⟩·⟨σ⟩` was treated as a single average to derive an equation for).

**Fixed:** `_is_leaf_average` predicate checks for `operation === sym_average`
so products of averages are recursed into rather than pushed.

## `get_order` returned the operator-product order, not the leaf-moment order

**Symptom:** for cumulant-truncated forms like `⟨a⟩·⟨a*a⟩` (a 1st-order
moment times a 2nd-order moment), `get_order` returned 3 (the order of
`undo_average(⟨a⟩·⟨a*a⟩) = a*a*a`). This made `cumulant_expansion(_, 2)`
flagged-as-untruncated even when the result was already a product of
sub-order-2 leaf moments.

**Fixed:** `get_order` only follows the `undo_average` branch when the head
is `sym_average` (leaf average). Products are now `maximum(child orders)`.

## `cumulant(op, n)` was wrong: only summed partitions of size exactly `n`

**Symptom:** `cumulant(a*b)` returned `-⟨a⟩⟨b⟩` instead of
`⟨a*b⟩ - ⟨a⟩⟨b⟩`. The `_term_cumulant` loop iterated only
`partitions(ops, n)` (exactly `n` blocks) instead of `k = 1:n`.

**Fixed:** loop is now `for k in 1:n`.

## Forward measurement-noise drift double-subtracted the cumulant mean

**Symptom:** the noise expression had `-2·⟨J† + J⟩·⟨op⟩` where it should
have had `-⟨J† + J⟩·⟨op⟩`. Manifested as a factor-of-2 mismatch in any
analytic comparison test (single-atom damped cavity etc.).

**Fixed:** `src/noise.jl::_noise_drift_one` no longer multiplies `avg_term`
by 2. Backward path was already correct.

## `complete!` did not accept `order=` kwarg

**Symptom:** master allowed `complete(eqs; order=2)` to override the
cumulant-truncation order at closure time. v1 dropped this — the order had
to be fixed at `meanfield(...; order=...)` construction.

**Fixed:** `complete(eqs; order=N)` now runs `cumulant_expansion(_, N)`
first, then completes.

## SQA v0.5 strict index conflict in collective indexed Hamiltonians

**Symptom:** master test_indexed_meanfield's collective-decay JC setup
(`H` builds `Σ(g(i)*(a'*σ(1,2,i) + ...), i)` and ops use `σ(_,_,k)`) errors
with `ArgumentError: Summation index i appears in both factors` during the
inner `commutator(im*H, op)` step in v1, even though master accepted it.

**Status:** appears to be a SQA-level change (v0.5 stricter index hygiene)
rather than QC code. The same construction worked under SQA in master.

**Workaround:** use scalar (non-indexed) ops in collective tests, or
explicitly declare index disjointness via `assume_distinct_index` if that
mechanism exists in v0.5.

## `Σ(..., i, j; non_equal=true)` kwarg removed

**Symptom:** master constructed `DoubleSum(_, i, j; non_equal=true)` to
build dipole-dipole sums excluding the diagonal. SQA v0.5 removed the
`non_equal` keyword.

**Status:** SQA-level API change. QC consumers must build the `i != j`
exclusion manually (e.g. via `Σ` with explicit `[non_equal_indices]`
argument).

## `scale(::NoiseMeanFieldEquations)` not defined

**Symptom:** master scaled indexed noise systems via `scale(stoch_eqs)`.
v1 only has `scale(::MeanFieldEquations)`.

**Status:** missing. Would need to scale both `eqs.equations` and
`eqs.noise_equations` consistently.

## No public `conj` / `_conj` of an `Average`

**Symptom:** master had `QuantumCumulants._conj(average(2im*a'*σ))` returning
`(-2im)*average(a*σ')` (i.e. it pushed the conjugation inside the average
and conjugated the operator). v1's `conj(average(...))` returns a literal
`conj(...)` wrapper; `SQA.qconj(average(...))` errors with `MethodError:
no method matching conj(::QAdd)` because `qconj` on a SymReal-typed average
tries to `qconj` each Mul-arg, which recurses into `qadjoint` on a SymReal
and falls into Symbolics' generic `conj`.

**Status:** the master test_average assertion is `@test_skip`-marked in
test/cumulant_test.jl pending a v1 public surface for this operation.

## Deferred from CHANGELOG (not regressions per se, just not yet ported)

- `evaluate(eqs; limits=(N=>n,))` to unroll indexed sums to fixed-N
- `initial_values(eqs, ψ::Ket)` for spin-coherent / Fock state initialisation
- Rate-matrix collective decay (`J = [[J1, J2]]; rates = [[Γ11 Γ12; Γ21 Γ22]]`)
- `meanfield_backward(...)` (use `meanfield(...; direction=Backward())`)
- `modify_equations(eqs, f)` / `translate_W_to_Y(eqs)` (retrodiction helpers)
- `find_operators(h, order)` returning a deduplicated, canonically-ordered
  operator set so closure tests can assert exact counts
