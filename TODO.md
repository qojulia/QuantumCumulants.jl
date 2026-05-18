# v1 regressions discovered during test porting

These were identified while porting the master test suite. Each item is a
behavior the master branch handled but the v1 rewrite either drops, weakens,
or gets wrong. Listed in rough priority order.

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

## Deferred from CHANGELOG (not regressions per se, just not yet ported)

- `evaluate(eqs; limits=(N=>n,))` to unroll indexed sums to fixed-N
- `initial_values(eqs, ψ::Ket)` for spin-coherent / Fock state initialisation
- Rate-matrix collective decay (`J = [[J1, J2]]; rates = [[Γ11 Γ12; Γ21 Γ22]]`)
- `meanfield_backward(...)` (use `meanfield(...; direction=Backward())`)
- `modify_equations(eqs, f)` / `translate_W_to_Y(eqs)` (retrodiction helpers)
- `find_operators(h, order)` returning a deduplicated, canonically-ordered
  operator set so closure tests can assert exact counts
