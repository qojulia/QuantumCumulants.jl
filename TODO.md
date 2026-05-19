# v1 rewrite, remaining work

Outstanding work on the v1 rewrite. See [DESIGN.md](DESIGN.md) for the
target architecture and [CHANGELOG.md](CHANGELOG.md) for what has landed.

Current state: 625 pass + 1 broken / 626 total (`make test 2>&1 | tee /tmp/maketest.log`).
All 14 examples run end-to-end. All non-SQA master tests are ported,
test-strengthening pass is done (noise, parameters, higher_order,
meanfield, indexed_scale, indexed_correlation, indexed_scope all
assert physics, not just pipeline shape). Remaining shape-only files
(completion, cumulant, find_operators, numeric_conversion, scaling)
are correctly scoped to algebraic/structural properties; downstream
physics is covered by the strong files. Remaining work is
architectural cleanup and SQA-side feature gaps.

### Known representation-diff blocker

| File | What's not reachable | Why |
|---|---|---|
| `test/parameters_test.jl::detuned two-level commutator sign` (currently `@test_broken`) | `iszero(simplify(rhs + 2im * Δ * ⟨s₁₂⟩))` | SymbolicUtils stores `2im * Δ` as `complex(0, 2Δ)` literal vs our factored form. Both code-gen to the same numeric values but the symbolic difference doesn't simplify. Leave `@test_broken` until SymbolicUtils representation settles. |

## Open v1 feature gaps that block deeper test ports

- **`evaluate(eqs; h = [hilb])` per-Hilbert-space filter**: master
  accepts `h::Vector` to evaluate only specific Hilbert subspaces. v1
  emits a warning and ignores `h`.
- **`scale(eqs; h = [k])` per-Hilbert-space scaling**: master accepts
  `h::Vector` to scale a subset of Hilbert subspaces, leaving others as
  symbolic indices.
- **`scale` per-atom rate coefficients pick up an extra `N`**: scaling
  the rate-equation for ⟨σ_k₂₂⟩ in the indexed JC laser at order=1
  produces `R + N*(-R-Γ)*⟨σ_k₂₂⟩` instead of the per-atom form
  `R + (-R-Γ)*⟨σ_k₂₂⟩` (the LHS is a *specific* atom's population, not
  the sum over atoms, so the `N` factor on the decay term is wrong). At
  `N=1` the formula reduces correctly. Until fixed, physics steady-state
  assertions must go through `evaluate(eqs_c; limits = (N => 1))`
  rather than `scale`; see `test/indexed_correlation_test.jl::order=1
  laser steady state via evaluate(N=>1)`. **Root cause:** SQA's
  `meanfield` derivation puts a `Σ_k` sum-scope on the RHS σ22 term
  even though the LHS uses the same free `k` for a specific atom: at
  the operator level the inner average is `Σ(k=1:N) σ_k₂₂` (sum scope
  `Index[k]`) while the LHS is just `σ_k₂₂` (no scope). `scale` then
  correctly attaches an `(N - 0) = N` prefactor via
  `_sum_scope_prefactor`. The fix lives in SQA: when a jump
  `σ(2,1,k)` and an observable `σ(2,2,k)` share the *same* free index,
  the Lindblad contraction should produce a per-atom term (no sum
  scope) rather than a bound-`k` sum. Workaround in tests: use
  `evaluate(N=>1)` or use different free indices for ops vs jumps.
- **`complete!` mixed-order parity**: for `order = [1, 2]` on the
  indexed JC laser, v1 derives more equations than master (23 vs 8).
  Both are valid closures; the assertion-count divergence blocks
  straight equation-count ports.
- **`Spectrum` numerical stability** for higher-order cumulants. Affects
  spectral-equality assertions; ODE-steady-state convergence still works.
- **`numeric_average(op, state; level_map=…)` kwarg passthrough**:
  master's `initial_values(eqs, ψ::Ket; level_map=…)` translates symbolic
  level names (`:g`, `:e`) to integer indices in `NLevelBasis`. v1
  delegates to `SQA.numeric_average` which dropped the `level_map`
  kwarg. Either re-add `level_map` to SQA or document that symbolic
  levels must be replaced with integer levels in v1 user code.

## Architectural follow-ups from [DESIGN.md](DESIGN.md)

File-size and maintainability wins, NOT user-facing feature changes.

- ~~**Step 4**: consolidate `_canonical_key` / `_build_canonical_indices`
  into a single source of truth shared by `find_missing`, `scale!`,
  `evaluate`, `to_system`.~~ *Done, and the framing oversold the
  duplication.* `_build_canonical_indices` is already centralised in
  `completion.jl` and called from `scaling.jl` and `evaluate.jl`.
  `_canonical_key` (alpha-equivalence, position-indexed) and
  `_scale_qadd` + `_scale_state_key` (permutation collapse to
  canonical-first) implement *different policies* and shouldn't be
  unified. The only genuine duplication was the index-collection loop
  inside `_canonical_key` vs `_free_op_indices`; `_free_op_indices` now
  lives in `completion.jl` and is shared. `to_system` doesn't use this
  machinery.
- **Step 5**: audit `correlation.jl` (520 lines) and `mtk.jl` (320
  lines) for the same scalar-first-then-patched pattern. Target:
  correlation drops to ≤250 lines by reusing `to_system`.
- **SQA additions**: `strip_sum_scope`, `set_index`,
  `canonicalise_undetermined`, `enumerate_sum`, batched `change_index`,
  `pairwise_distinct`. Promote to SQA's public API and have QC call
  them instead of carrying local equivalents.
- **SQA op-scalar product overload**:
  `*(::BasicSymbolic{<:SymReal}, ::QAdd)` and the reverse are missing,
  so multiplying an averaged quantity into a `QAdd` requires wrapping
  in `Symbolics.Num(...)` (e.g. `measurement_retrodiction_test.jl`'s
  `f_measure` callback for Kalman drives). Adding the overload would
  let `modify_equations` callbacks be written naturally.

## Definition of "done" for this branch

- [ ] Ported tests assert master-equivalent (or stronger) physics, not
      just pipeline shape.
