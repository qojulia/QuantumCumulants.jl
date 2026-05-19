# v1 rewrite, remaining work

Outstanding work on the v1 rewrite. See [DESIGN.md](DESIGN.md) for the
target architecture and [CHANGELOG.md](CHANGELOG.md) for what has landed.

Current state: 625 pass + 1 broken / 626 total (`make test 2>&1 | tee /tmp/maketest.log`).
All 14 examples run end-to-end. All non-SQA master tests are ported and
the bound-index coefficient orphaning bug is resolved; the remaining
work is test-strengthening and architectural cleanup.

## Test strengthening (next pass)

Most ported tests cover the same scenarios as master *in pipeline shape*
but do NOT verify the same physics agreement that master asserted. The
next pass strengthens them so v1 demonstrably reproduces master where v1
has the features, and explicitly marks the gaps where it doesn't.

### Helpers already in tests, reuse them

- `test/parameters_test.jl::_is_zero(x)`: symbolic-vs-Bool zero detection:
  ```julia
  _is_zero(x::SQA.QAdd) = iszero(x)
  _is_zero(x::Number) = iszero(x)
  function _is_zero(x::SymbolicUtils.BasicSymbolic)
      s = simplify(x; expand = true)
      s isa Number && return iszero(s)
      SymbolicUtils.isconst(s) && return iszero(s.val)
      return isequal(s, 0)
  end
  ```
  Lift to `test/common.jl` or duplicate per file.
- `test/higher_order_test.jl::ϕ` and `phase_invariant`: the
  phase-invariance filter pattern (counts adjoint excess per term, used
  with `complete(...; filter_func=phase_invariant)`). Same pattern in
  `test/indexed_correlation_test.jl` (commented-out master form would
  need it for order=2).
- `parameter_map(eqs, [g(i) => 0.1, κ => 1.0, …])` in [src/mtk.jl](src/mtk.jl):
  user-facing pmap helper for `IndexedVariable` parameters. Scalar values
  broadcast to N-element vectors; per-atom vectors pass through.

### Concrete strengthening patterns (use these, don't reinvent)

**ODE numerical equality between two derivation paths:**
```julia
SQA = QuantumCumulants.SecondQuantizedAlgebra
sol_a = solve(prob_a, Tsit5(); abstol = 1e-10, reltol = 1e-10)
sol_b = solve(prob_b, Tsit5(); abstol = 1e-10, reltol = 1e-10)
# For scale-vs-evaluate: scale collapses to 1 atom-state, evaluate
# produces N atom-states by symmetry; use get_solution on each.
obs_op_a = SQA.undo_average(eqs_a.states[k_a])
obs_op_b = SQA.undo_average(eqs_b.states[k_b])
val_a = get_solution(sol_a, obs_op_a, eqs_a)(sol_a.t[end])
val_b = get_solution(sol_b, obs_op_b, eqs_b)(sol_b.t[end])
@test isapprox(val_a, val_b; atol = 1e-6)
```

**Reference-RHS equality for noise equations (master pattern):**
```julia
@variables Δ::Real κ::Real η::Real
H = Δ * a' * a
me = meanfield(a, H, [a]; rates = [κ], efficiencies = [η])
expected = average(-1im * Δ * a - 0.5κ * a)
@test _is_zero(simplify(me.equations[1].rhs - expected; expand = true))
expected_noise = √(η * κ) * (
    average(a' * a) + average(a * a) -
    average(a)^2 - average(a) * average(a')
)
@test _is_zero(simplify(me.noise_equations[1].rhs - expected_noise; expand = true))
```

**Steady-state cumulant-order convergence:**
```julia
he4 = complete(meanfield(a'*a, H, J; rates=rates); order=4, filter_func=phase_invariant)
he6 = complete(meanfield(a'*a, H, J; rates=rates); order=6, filter_func=phase_invariant)
# build sys, prob, solve to t=long
n_ss_4 = real(get_solution(sol4, a'*a, he4)(sol4.t[end]))
n_ss_6 = real(get_solution(sol6, a'*a, he6)(sol6.t[end]))
@test abs(n_ss_4 - n_ss_6) / abs(n_ss_6) < 0.05
```

**Independent reference for indexed_mixed_order (brute-force N=3):**
Construct the 3-atom Hilbert space explicitly (`tensor(ha for _ in 1:3)`
plus the cavity), build the master operator, evolve the density matrix,
compare. Use `QuantumOpticsBase.timeevolution.master` or similar. If too
involved, ODE-solve `evaluate(eqs; limits=(N=>3))` and compare against
ODE-solve of scale of the same system; both should converge to the same
steady-state observable for permutation-symmetric initial conditions.

**SDE pipeline for measurement_retrodiction:**
Use `StochasticDiffEq.SDEProblem` + `EM()`. v1 has SDESystem code-gen
wired (see `examples/heterodyne_detection.jl` and
`examples/retrodiction_homodyne.jl`). The Kalman forward+backward
agreement is the master assertion; reproducing it requires
`meanfield_backward` to produce the correct backward dynamics (verify
via `direction = Backward()`).

### Per-file audit (still to do)

| File | Master's strongest assertion | What we have now | Action |
|---|---|---|---|
| `test/indexed_correlation_test.jl` | order=2 indexed JC with phase-invariant filter, `evaluate(corr, ...)`, `split_sums` | order=1 only, no filter | Try order=2 (may be slow). If the v1 `complete!` mixed-order issue blocks order=2, add a ParallelTestRunner timeout and assert "completes in N seconds" as a regression marker. |

### Tests where master's assertion is NOT reachable (representation diff)

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
Gated on the test-strengthening above landing first.

- ~~**Step 3**: unify `meanfield.jl` three derivation paths (forward,
  noise-forward, backward) into one parameterised `derive()`. ~120
  lines saved.~~ Done: `_meanfield_deterministic` and the (Forward,
  Backward) `_meanfield_noise` methods now share `_build_op_drift`,
  `_build_avg_drift_eqs`, and `_finalize_noise_eqs`; the two
  `_meanfield_noise` methods collapse to one. Saved ~32 lines and
  removed three near-identical drift loops.
- **Step 4**: consolidate `_canonical_key` / `_build_canonical_indices`
  into a single source of truth shared by `find_missing`, `scale!`,
  `evaluate`, `to_system`. Currently `scaling.jl` and `completion.jl`
  carry near-duplicate machinery.
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
