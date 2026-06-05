# TODO

## cavity_antiresonance: collective indexed dissipation (the one remaining suite failure)

`QCNew/test/examples_regression_test.jl` `@testset "cavity_antiresonance_indexed"` (line ~604)
is the single failing test (suite is otherwise 852 pass / 0 fail / 0 broken). It uses a
`DoubleIndexedVariable` rate matrix `Γ(i,j)` for the indexed jump `σ12_i` (collective decay)
plus a coherent dipole-dipole term `Σ_i Σ_{j≠i} Ω(i,j) σ21_i σ12_j` with `Ω` declared
`identical=false`.

### Symptom

`evaluate(eqs; limits=(N=>2))` followed by `System` throws

```
Could not evaluate value of parameter i_2_2. Missing values for variables in expression i_2_2.
```

i.e. a bound summation index (`i`/`j`) leaks into the generated ODE as an MTK parameter,
because the matrix-rate / dipole coefficient is left outside its sum with a dangling index.

### Root cause

QCNew does not implement **collective indexed dissipation**. In
`QCNew/src/operator_drift.jl`, `_lindblad_rhs` only special-cases a *concrete* `AbstractMatrix`
rate. A symbolic `DoubleIndexedVariable` rate `Γ(i,j)` on a singly-indexed jump `σ12_i` falls
through to the scalar path, which multiplies `Γ(i,j)` as a scalar and sums only over `i`,
leaving `j` dangling.

The correct form (see master QuantumCumulants v0.4.3 `indexed_master_lindblad`,
`src/index_meanfield.jl:134`) is the cross-jump dissipator

```
Σ_{i,j} Γ(i,j) D[σ_i, σ_j]
  = Σ_{i,j} (Γ(i,j)/2) ( σ_i† [O, σ_j] + [σ_i†, O] σ_j )
```

with `J†` placed at the rate's first index and `J` at its second.

### What is needed (two layers)

1. **Operator layer (`_lindblad_rhs`, operator_drift.jl).** Detect a 2-arg
   `DoubleIndexedVariable` rate on an indexed jump and emit the collective double sum above
   (build the partner `Index` from the rate's second arg via
   `SQA.Index(name, range, space_index, Num(sym))`, then `SQA.change_index` the jump onto it,
   and `SQA.Σ` over both indices). A drafted version of this exists in the session history; it
   removed the leaks and moved the numeric value toward master's, so the operator-layer piece
   is largely understood.

2. **Moment layer (`average_and_truncate`, moments.jl).** The index-dependent coefficient must
   stay INSIDE its sum so it follows the sum's diagonal split. Two concrete requirements:
   - Coherent `Ω` term: `Σ_i Ω(i,k) σ_k σ_i σ_k` collapses on the diagonal `i=k` to
     `Ω(k,k) σ_k`, which must be **zero** (Ω is `identical=false`). A coefficient applied
     outside the sum keeps a free `i` and never vanishes. Routing `coeff·op` through
     `SQA.Σ(...)` before `average` produces the correct `0`.
   - Collective `Γ` term: the off-diagonal sum constraint (`i_2 ≠ j` on `Γ(i_2,j) σ_j`) lives
     on the free LHS index `i_2`, which is in the coefficient, not the operator, so NE pairs
     must be carried using coefficient indices too (else the diagonal `j=i_2` double-counts).

   A drafted coefficient-inside-sum rewrite of `average_and_truncate` (split the coefficient
   into additive parts, sum each over op-used ∪ coeff-used indices, carry NE accounting for
   coefficient indices) was attempted but **reverted** because it (a) regressed
   `quality/quality_test` (ExplicitImports), and (b) changed cavity's own `complete` count
   from 4 to 5 (a spurious state). A correct version must be guarded so it does not alter the
   already-passing indexed tests (e.g. restrict to `DoubleIndexedVariable` coefficients, and
   verify no spurious states appear).

3. **SQA double-sum diagonal (likely a SecondQuantizedAlgebra fix).** Even with 1+2, the
   collective Γ **diagonal** came out with coefficient `-1` instead of master's `-1/2`: SQA's
   `Σ` over the double sum appears to double-count the `i=j=i_2` point. Master's SQA (0.4.5)
   gets `-1/2` from the same `∑(c, ind1, ind2)`, so this is a regression in the pinned
   `SecondQuantizedAlgebra` `redesign-v2` branch's double-sum diagonal split. Reproduce with:

   ```julia
   # in QCNew/test env
   using QCNew; const SQA = QCNew.SecondQuantizedAlgebra
   # build Σ_{i,j} Γ(i,j) D[σ_i,σ_j] acting on σ12_k and compare the diagonal
   # coefficient of ⟨σ12_k⟩ against master (should be -1/2 · Γ_kk).
   ```

### Ground-truth check

Master QC v0.4.3 for this exact system: `complete` = 3 eqs, `evaluate(N=>2)` = 5 eqs,
`abs2(⟨a⟩)(t=30)` ≈ `0.001984`. The current lock `0.0019889934766074693` is itself ~3e-3 off
master, so when this is fixed, re-derive the value lock against master rather than trusting the
existing number. Do NOT weaken the test (no `@test_broken`/skip) without re-deriving the
correct count and value first.

### Definition of done

`make test` (or the QCNew suite) green with `cavity_antiresonance_indexed` passing on a
re-derived count and value, and no regression to the other 852 tests (especially
`quality/quality_test`, `indexed_scope`, `indexed_meanfield`, `dicke`, `multilevel`).

Resolved investigation write-ups are archived in
[INVESTIGATION_NOTES.md](INVESTIGATION_NOTES.md).
