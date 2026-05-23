# TODO

Open work items for the v1 rewrite. Each entry names the failure mode and where it
shows up so the fix can be verified end to end.

## Numerics / MTK bridge

### MTK v10 rejects unused-parameter dict entries

Master's MTK (v8/v9) silently ignored entries in the `Dict(ps .=> p0)` that
referred to parameters not actually in the compiled system. MTK v10 asserts:

```
AssertionError: Expected an `Initial` parameter to exist for variable `Ω_2_2`,
but did not find one.
```

Two cases that surface this:

- `examples/waveguide.jl`: builds `ps` from `[Ωp(i, j) for i in 1:M for j in 1:M]`,
  but the Hamiltonian excludes `i == j`, so `Ω_i_i` is never a system parameter.
- `examples/heterodyne_detection.jl`: deterministic comparison uses
  `System(scaled_eqs; noise=false)`. The flag drops the Brownian column, which
  also drops `ξ` (it appeared only in the efficiency factor `sqrt(ξ κ)` of the
  noise term). The original `p .=> p0` still passes `ξ`, which MTK rejects.

Fix sketch: provide a helper that filters a user-built parameter dict to the
intersection with `parameters(sys)`, or accept a relaxed-mode keyword on
`ODEProblem`. Without it, every example author has to know which parameters
are actually live in the compiled system.

### `AvgFunc not callable` at solve time after `evaluate` of `DoubleIndexedVariable`

With `parameter_map`/`evaluate` now flattening `Γ(i, j)` into a 2-D
Symbolics-array parameter `Γ[i, j]` (fixed in this branch), `ODEProblem`
builds cleanly but `solve` crashes during RHS codegen with:

```
MethodError: objects of type SecondQuantizedAlgebra.AvgFunc are not callable
```

The `getindex(Γ, i, j)` term inside an averaged operator product is being
wrapped through `AvgFunc` and then invoked as a function. Likely a
`_safe_substitute` / averaging-pass issue when the rewrite from
`Γ(i_1, i_2)` to `Γ[1, 1]` happens inside `⟨…⟩`.

Surfaces in: `examples/cavity_antiresonance_indexed.jl` (transmission sweep).

Fix sketch: trace the codegen of an `⟨op · Γ[1,1]⟩` term through
`MTK.System(::MeanFieldEquations)`'s `_safe_substitute` / `_collect_params!`
passes and ensure scalar `getindex(arr, i, j)` is hoisted out of the
averaged-product wrapper.

### `σᵢᵢ · σⱼⱼ = 0` simplification missing before cumulant expansion

Orthogonal-projector idempotency (`σ_22·σ_11 = 0` always) is not applied at
the operator level before the second-order cumulant expansion. Result: the
RHS of `⟨σ_22⟩` in a scaled, phase-invariant superradiant laser comes out as

```
d/dt ⟨σ22⟩ = -R·⟨σ22⟩ + R·⟨σ22⟩·⟨σ11⟩ + … + (terms with ⟨σ22⟩² and ⟨σ22⟩³)
```

instead of the expected `R·⟨σ11⟩ - (Γ + ν)·⟨σ22⟩ + …`. Cumulant truncation
then approximates `⟨σ_22·σ_11⟩ ≈ ⟨σ_22⟩·⟨σ_11⟩` which is numerically
non-zero even though the exact value is identically 0. The σ22² and σ22³
terms are the same idempotency leak. `⟨σ_11⟩` also appears as a redundant
state on completion (master eliminated it via `σ_11 + σ_22 = 1`).

Surfaces in: `examples/superradiant_laser_indexed.jl` (trajectories run but
σ22 stays at 0 from the default `u0 = 0` initial state).

Fix sketch: per CLAUDE.md's "ground-state projectors stay atomic" guidance,
this should not auto-expand into `1 - Σ σ_kk`. Instead, ensure SQA's
`_canonicalize!` collapses `σ_ii · σ_jj → δ_ij σ_ii` for transitions on the
same site (`acts_on` match), or add the simplification at the meanfield /
scale stage before the cumulant truncation runs. Note that
`assume_distinct_index` is the right primitive when the two factors are on
different atoms; that already implies `σ_ii^{(p)} · σ_jj^{(q≠p)} ≠ 0`, but
same-site `σ_ii · σ_jj` must still reduce.

### Solver-accuracy follow-up: `conj(cosh(ξ))` in symbolic RHS

`examples/unique_squeezing.jl`: full vs. effective-model squeezing should
agree at large `N`, but the RHS carries unresolved `cosh(ξ)·conj(cosh(ξ))`
/ `cosh(ξ)*` (conjugate) factors even though `ξ` is real-typed. Check
whether MTK simplifies `conj(cosh(ξ)) → cosh(ξ)` for `ξ::Real` declared
via `@variables ξ`, and if not, either pre-simplify in QC's `System`
builder or document the workaround.

Surfaces in: `examples/unique_squeezing.jl`.

## Refactors deferred from the /simplify pass

### Replace `noise::Bool` kwarg on `MTK.System(::NoiseMeanFieldEquations)` with `drop_noise`

The `noise=false` path on `MTK.System(::NoiseMeanFieldEquations; …)` branches
inside `_to_system_sde` (the Brownian column, the per-equation maketerm, and
the final `MTK.System(...)` constructor signature all differ). Both reviewers
flagged this as the noise/drift split being wedged into one function instead
of two.

Cleaner shape: add a `drop_noise(::NoiseMeanFieldEquations) → MeanFieldEquations`
accessor that rewraps `.equations` / `.states` / `.operators` / `.hamiltonian`
/ `.jumps` / `.jumps_dagger` / `.rates` / `.iv` / `.order` / `.direction` into
a `MeanFieldEquations`. Users then write
`System(drop_noise(eqs); name=:sys)` instead of `System(eqs; name=:sys, noise=false)`,
and `_to_system_sde` collapses back to a single straight-line function.

Surfaces in: `examples/heterodyne_detection.jl` (deterministic comparison
section).

Fix sketch: add the accessor to [src/equations.jl](src/equations.jl) next to
the `MeanFieldEquations` / `NoiseMeanFieldEquations` definitions, then drop
the `noise::Bool` kwarg in [src/mtk.jl](src/mtk.jl) and revert
`_to_system_sde` to its single-branch form.

### Hoist `T_uw` / `signed_rhs` out of the noise=true / noise=false branches in `_to_system_sde`

Both branches compute the same `T = typeof(...)` and the same `signed_rhs`
(sign flip on `rhs` when `direction == Backward`). Only the `noise_term`
construction and the `_collect_params!(ps_set, noise_rhs, …)` call are
noise-specific. Hoisting cuts the per-equation body in half and removes one
of the two `_collect_params!(ps_set, rhs, …)` calls.

Not done because the simpler accessor refactor above subsumes it (the
function would no longer have two branches at all).

## Examples status

All twelve docs examples run to the derivation step. Numerical / plotting
step still blocked on the items above:

- `examples/superradiant_laser_indexed.jl`: runs end to end; trajectory still
  anomalous due to the σᵢᵢ·σⱼⱼ idempotency issue above.
- `examples/filter-cavity_indexed.jl`: runs end to end.
- `examples/waveguide.jl`: blocked on MTK v10 rejecting `Ω_i_i` (the diagonal
  the Hamiltonian never references).
- `examples/heterodyne_detection.jl`: blocked on MTK v10 rejecting `ξ` from
  the deterministic (`noise=false`) system's parameter dict.
- `examples/cavity_antiresonance_indexed.jl`: blocked on `AvgFunc not callable`
  at solve time (downstream of the `DoubleIndexedVariable` flattening fixed
  in this branch).
- `examples/retrodiction_homodyne.jl`: still commented out of the docs build
  (`docs/make.jl`); depends on the `Backward()` SDE path and `modify_equations`,
  which have not been ported.
