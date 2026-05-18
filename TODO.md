# v1 regressions discovered during test porting

These were identified while porting the master test suite. Each item is a
behavior the master branch handled but the v1 rewrite either drops, weakens,
or gets wrong. Listed in rough priority order.

## RHS rebuild leaves `complex(re, im)` literal calls in equations

**Symptom:** equation RHSes are stored as e.g. `⟨a⟩ * complex(-κ/2, -Δ)`
because `complex(BasicSymbolic{SymReal}, BasicSymbolic{SymReal})` materialises
a literal `complex(...)` call instead of unfolding to `re + im*1im`.

**Consequence:** `simplify(rhs - analytic_target; expand=true)` will not
reduce to 0 even when the difference is identically zero, because
`complex(...)` is opaque to Symbolics' polynomial rules.

**Workaround applied:** `src/meanfield.jl` now runs
`_rewrite_complex_literals` on every RHS before storing it
(`_meanfield_forward`, `_meanfield_noise` both directions, noise drift +
noise eq). Tests use a local `_iz` helper that also runs the rewrite before
`simplify`.

**Why it's still on this list:** the rewrite is reactive. The right fix is
that v1's `Symbolics`/SQA glue should never produce the literal `complex(...)`
call in the first place (or `simplify` should expand it). Once a user
constructs an analytic comparison target like `-1im*Δ*⟨a⟩ - 0.5κ*⟨a⟩`,
they hit the same problem and need to know about the rewrite.

## `complete!()` does not preserve conjugate-pair states needed by codegen

**Symptom:** for two-level systems where initial ops are
`[σ(:e,:g), σ(:e,:e)]`, `complete!()` exits with only those two states because
`find_missing` correctly recognises `⟨σ_12⟩ = conj(⟨σ_21⟩)` and skips it. The
RHS of the `⟨σ_22⟩` equation still references `⟨σ_12⟩`, but it's not a state.

**Consequence:** `to_system(eqs)` then `solve(prob, Tsit5())` raises
`MethodError: objects of type SecondQuantizedAlgebra.AvgFunc are not callable`
at runtime — the codegen emits a literal `avg(σ_12)` call because there is
no state variable to substitute.

**Root cause:** v1's `to_system` codegen substitutes RHS averages via
`_avg_to_var_dict(eqs)`, which only knows about averages that appear in
`eqs.states`. It has no rule to map `avg(adj(op))` to `conj(state_of_op)`.

**Master's behavior:** master also tracked only canonical states (see
`vs_adj` in `master:src/utils.jl`), but its System codegen knew how to
express conjugate observables via the primary states. This is a real gap.

**Workaround:** user must call `complete!(eqs; get_adjoints=true)` (the
current default in v1 to keep ODE flows working), which inflates closure
counts (e.g. V-level: 32 vs master's 16).

**Proper fix:** teach `_avg_to_var_dict` / `to_system` codegen to emit
`avg(adj(op)) → conj(state_var_of(avg(op)))` substitution rules.

## `complete!()` over-counts the canonical closure

**Symptom:** master's V-level test closes at 16 equations. v1's `complete!`
with `get_adjoints=false` produces 18 — two extra second-order moments where
master's ordering canonicalisation picks a single representative of a
conjugate pair but v1 doesn't.

**Consequence:** master-equivalent count assertions need a tolerance band.

**Root cause:** SQA v0.5 normal-ordering / `find_operators` heuristics
differ subtly from the old in-tree QC implementation. When both `⟨a*σ⟩`
and `⟨a'*σ'⟩` appear in a single RHS scan, master would pick one based on
its hashing/`<`-ordering rules; v1's `_collect_missing!` picks based on
which it encounters first, which depends on the SQA term order.

**Proper fix:** match master's `lt_reference_order` choice when picking
which conjugate of a pair becomes the canonical state.

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
