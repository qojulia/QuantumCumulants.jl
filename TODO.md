# v1 rewrite, remaining work

Outstanding work on the v1 rewrite. See [DESIGN.md](DESIGN.md) for the
target architecture and [CHANGELOG.md](CHANGELOG.md) for what has landed.

Current state: 658 pass / 658 total (`make test 2>&1 | tee /tmp/maketest.log`).
All 14 examples run end-to-end. All non-SQA master tests are ported,
test-strengthening pass is done (noise, parameters, higher_order,
meanfield, indexed_scale, indexed_correlation, indexed_scope all
assert physics, not just pipeline shape). Remaining shape-only files
(completion, cumulant, find_operators, numeric_conversion, scaling)
are correctly scoped to algebraic/structural properties; downstream
physics is covered by the strong files. Remaining work is feature
gaps and a small set of SQA primitive promotions.

## Open v1 feature gaps that block deeper test ports

- **`numeric_average(op, state; level_map=…)` kwarg passthrough**:
  master's `initial_values(eqs, ψ::Ket; level_map=…)` translates symbolic
  level names (`:g`, `:e`) to integer indices in `NLevelBasis`. v1
  delegates to `SQA.numeric_average` which dropped the `level_map`
  kwarg. Either re-add `level_map` to SQA or document that symbolic
  levels must be replaced with integer levels in v1 user code.

## SQA primitive promotions

- **SQA additions**: `strip_sum_scope`, `set_index`,
  `canonicalise_undetermined`, `enumerate_sum`, `pairwise_distinct`.
  Promote to SQA's public API and have QC call them instead of carrying
  local equivalents.
- **SQA op-scalar product overload**:
  `*(::BasicSymbolic{<:SymReal}, ::QAdd)` and the reverse are missing,
  so multiplying an averaged quantity into a `QAdd` requires wrapping
  in `Symbolics.Num(...)` (e.g. `measurement_retrodiction_test.jl`'s
  `f_measure` callback for Kalman drives). Adding the overload would
  let `modify_equations` callbacks be written naturally.

## Definition of "done" for this branch

- [ ] Ported tests assert master-equivalent (or stronger) physics, not
      just pipeline shape.