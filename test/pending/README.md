# Pending master test ports

Tests in this directory are master-branch tests that have NOT been ported to
the v1 surface. Each is one of:

## SQA-level (should be ported into SecondQuantizedAlgebra.jl, NOT here)

Per [../../CLAUDE.md](../../CLAUDE.md): "Operator-algebra/index/sum unit tests
belong in **SecondQuantizedAlgebra.jl**, not here." These files exercise
operator algebra, indexing, sums, or low-level construction. The right home
is SQA's test suite. They are kept here as a porting source-of-truth.

- `average_sums_test.jl`: `IndexedAverageSum`, `IndexedAverageDoubleSum`,
  `SpecialIndexedAverage`, `NumberedOperator`. All SQA-level.
- `cluster_test.jl`: `ClusterSpace` is **removed** in v1. Reference only.
- `code_quality_test.jl`: Superseded by `test/quality/quality_test.jl` and
  JET in `test/quality/`.
- `double_sums_test.jl`: `DoubleSum`, `SingleSum`, multi-index `Σ` ordering.
  SQA-level.
- `fock_test.jl`: `Destroy`/`Create` commutators, Fock algebra. SQA-level.
- `index_basic_test.jl`: `Index`, `IndexedOperator`, `IndexedVariable`,
  `SingleSum`/`DoubleSum`, `change_index`, `order_by_index`, `reorder`.
  SQA-level.
- `nlevel_test.jl`: N-level operator algebra (`Transition` composition,
  `ha.ground_state`, level identities). SQA-level. Some QC-level workflow
  bits should migrate to `test/systems/multilevel_test.jl` if not already
  covered.
- `spin_test.jl`: Mixes SQA-level Pauli/Spin algebra with QC-level Dicke
  meanfield + ODE. The QC parts are already covered in
  `test/meanfield_test.jl` (Pauli closed two-spin). Algebra parts go to SQA.

## Blocked on v1 features still to land

These test QC surface features but cannot run yet because the corresponding
v1 implementation gap is open. See [../../TODO.md](../../TODO.md) for
status of each blocker.

- `higher_order_test.jl`: 4th/6th-order cumulant + `Spectrum`. Blocked on
  Spectrum numerical stability in v1.
- `indexed_filter_cavity_test.jl`: Blocked on `scale(eqs; h = [hilb])` and
  `evaluate(eqs; limits, h = [...])` per-Hilbert-space variants.
- `indexed_mixed_order_test.jl`: Blocked on the v1 `complete!`
  mixed-order behavior diverging from master (master derives 8 base
  equations, v1 derives 23 for the same input; both are valid closures
  under cumulant=[1,2] but the test-asserted count is master-specific).
- `indexed_scale_test.jl`: Blocked on `scale(eqs; h = [k])` per-Hilbert
  scaling.
- `measurement_backaction_indices_test.jl`: Indexed noise + measurement
  backaction. Blocked on `scale(::NoiseMeanFieldEquations)` indexed path
  parity.
- `measurement_backaction_indices_comparison_test.jl`: Comparison
  workflow; same blocker as above.

## Awaiting port (QC-surface; not blocked, just unported)

- `indexed_correlation_test.jl`: Should port directly. Indexed
  `CorrelationFunction` is unified into `CorrelationFunction` in v1.
- `measurement_retrodiction_test.jl`: `meanfield_backward` is implemented
  in v1; this should port after a `_master_lindblad_backward` audit.
- `numeric_conversion_test.jl`: `to_numeric`, `numeric_average`,
  `initial_values(eqs, ψ::Ket)` are exercised. Most of the surface is in
  SQA now; the QC bits (initial_values via Ket) are already implemented.
  A subset of this test belongs in v1 once `numeric_average` extensions
  for `OpticalState` are settled.
- `parameters_test.jl`: Symbolic parameter handling in `meanfield`. Mostly
  works in v1 but the master test uses `eqs.operator_equations[i].rhs`
  comparisons that v1 stores differently (SymReal symbolic vs QAdd).
- `scaling_test.jl`: Master's `scale` test. Partially superseded by v1
  `test/scaling_test.jl`. Remaining bits exercise indexed scaling beyond
  what v1 has yet (see indexed_filter_cavity_test).
