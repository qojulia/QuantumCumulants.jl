# Changelog

All notable changes to QuantumCumulants.jl will be documented in this file.

## [0.5.0]

Breaking release. Clean rewrite atop SecondQuantizedAlgebra v0.5 and ModelingToolkitBase. Migration from 0.4 is hand-driven; there are no deprecation shims.

### Removed

- `simplify` keyword from every QC API that exposed it: `meanfield`, `complete!`, `complete`, `cumulant`, `cumulant_expansion`, `evaluate`, `translate_W_to_Y`, `CorrelationFunction`, `Spectrum`. QC never runs `SymbolicUtils.simplify` on its derived RHS expressions; canonical-form post-processing is the user's responsibility (`Symbolics.simplify(eq.rhs; expand=true)`). The default `simplify=true` was a 99% time sink at higher cumulant orders (e.g. order=6 `CorrelationFunction` went from ~100s to ~10ms with this change) for no downstream numerical benefit, and the `Spectrum` kwarg was unused dead code.
- `ClusterSpace` removed in SQA v0.5. Use `Index` and `Σ` directly.
- `heisenberg` deprecation alias dropped. Call `meanfield` directly.
- `plotME` removed. Its target `EvaledMeanfieldEquations` no longer exists; latexify on the unified `MeanFieldEquations` covers the same use case.
- `get_scale_solution` removed. The unified `get_solution(sol, op, eqs)` handles both scalar and scaled-indexed systems.
- `subst_reds` removed. Redundant-conjugate substitution is folded into the `System(eqs)` codegen path automatically.
- `find_missing_sums` removed. Unified `find_missing` covers both scalar and indexed cases.
- Indexed-helper exports `value_map`, `split_sums`, `scale_term`, `eval_term`, `insert_index`, `order_by_index` removed. The equivalent transformations are expressed through SQA's `change_index` and `Σ` directly.
- Parameter-machinery exports `Parameter`, `CNumber`, `RNumber`, `RealParameter`, `Average` (type), and the `@cnumbers`/`@rnumbers`/`cnumbers`/`cnumber`/`rnumbers`/`rnumber` macros removed. Use `Symbolics.@variables x::Real` (or `::Complex`) for symbolic parameters.
- Indexed-sum types `IndexedAverageSum`, `IndexedAverageDoubleSum`, `SingleSum`, `DoubleSum`, `AvgSums`, `NumberedOperator`, `SpecialIndexedTerm` removed. They were SQA-internal in master; SQA v0.5 represents indexed sums uniformly as `QAdd` with `.indices` metadata, so the type zoo collapsed.

### Renamed / Unified

- `indexed_meanfield`, `indexed_complete`, `indexed_complete!`, `IndexedCorrelationFunction`, `scaleME`, `evalME`, `IndexedMeanfieldEquations`, `EvaledMeanfieldEquations`, `ScaledMeanfieldEquations` collapsed into the unified `meanfield`, `complete`/`complete!`, `CorrelationFunction`, `scale`/`scale!`, `MeanFieldEquations`.
- `MeanfieldNoiseEquations`, `IndexedMeanfieldNoiseEquations`, `BackwardMeanfieldNoiseEquations` collapsed into `NoiseMeanFieldEquations{...,Direction}` with `Direction <: Union{Forward, Backward}`.
- `meanfield_backward(...)` replaced by `meanfield(...; direction=Backward())`.
- `MeanfieldEquations` renamed `MeanFieldEquations` (capital F).

### Changed

- `meanfield(ops, H, J; ...)` is the single derivation entry. Noise is opt-in via `efficiencies=...`; retrodiction via `direction=Backward()`.
- `scale(eqs; h=Int[])` and `evaluate(eqs; limits, h=Int[])` accept a `h::Vector{Int}` of `space_index` values selecting which Hilbert subspaces are collapsed / unrolled. The empty default targets every subspace (previous behaviour). Hybrid systems can unroll some subspaces with `evaluate` and collapse others with `scale` in any order.
- `scale` drops additive constants from the LHS when SQA's commutator pipeline produces a `c + ⟨...⟩` shape (e.g. `⟨a*a'⟩` renamed to a same-index pair collapses to `1 + ⟨a'a⟩`). The constant-offset equation deduplicates against its normal-ordered sibling.
- Equation LHS in the user-facing struct is the raw `Average` `BasicSymbolic`. The `Differential(t)(u(t))` form is built only inside `ModelingToolkitBase.System(eqs; name)`, which QC extends to dispatch on `MeanFieldEquations`, `NoiseMeanFieldEquations`, and `CorrelationFunction`.
- The `System(eqs)` conjugate substitution builds `conj(avg_op(t))` via `SymbolicUtils.term(conj, var; type=Number)` rather than `Base.conj(var)`, so the symbolic `conj` node survives `mtkcompile` and `build_function`. Without this, every `⟨X⟩ - ⟨X†⟩` driving term silently zeroed on the compiled RHS because Symbolics simplifies `conj(::SymReal)` to identity. Related ModelingToolkit gap tracked at https://github.com/SciML/ModelingToolkit.jl/issues/4548.
- `cumulant_expansion` redistributes the sum-scope metadata (`SumIndices`/`SumNonEqual`) onto each factorised product term so `scale` can recover the per-index range prefactor after Wick factorisation. Bound indices that no factored leaf references keep SQA's "spurious bound idx, factor 1" convention.
- Quality gates (Aqua + ExplicitImports + CheckConcreteStructs + JET) are part of CI.
- Updated to ModelingToolkitBase 1.36, SymbolicUtils 4, Symbolics 7, SecondQuantizedAlgebra 0.5.

### Migration

- Replace `using ModelingToolkit` with `using ModelingToolkitBase`.
- Replace `indexed_meanfield(...)` / `indexed_complete!(...)` with the unified `meanfield(...)` / `complete!(...)`.
- Replace `IndexedCorrelationFunction(...)` with `CorrelationFunction(...)`.
- Replace `meanfield_backward(...)` with `meanfield(...; direction=Backward())`.
- Replace `@cnumbers x y` / `@rnumbers x y` with `@variables x::Real y::Real` (or `::Complex` as appropriate).
- Build numeric systems via `System(eqs; name=:sys)` then `mtkcompile(sys)`. QC extends `ModelingToolkitBase.System` to dispatch on `MeanFieldEquations`, `NoiseMeanFieldEquations`, and `CorrelationFunction`; the standalone `to_system` entry point is gone.
- Access trajectories via `get_solution(sol, op, eqs)(t)` instead of `sol[op]` or `get_scale_solution(...)`.
- For numeric initial values from symbolic level names, pre-translate the `:g`/`:e` labels to integer levels before calling `initial_values(eqs, ψ)`. The master `level_map` kwarg was dropped.

### Examples

All 14 examples from the v0.4 series are ported to v0.5 and pass an automated
smoke run:

- `damped_cavity.jl` (new minimal demo with 12-digit analytical agreement)
- `ramsey_spectroscopy.jl`, `optomechanical-cooling.jl`, `many-atom-laser.jl`,
  `excitation-transport-chain.jl` exercise the full meanfield to MTK to solve workflow.
- `mollow.jl`, `single-atom-laser-spectrum.jl`, `superradiant_laser_indexed.jl`,
  `filter-cavity_indexed.jl` exercise the symbolic + numeric pipeline.
  `Spectrum(c, ps)(ω, u_end, p0)` solves the Laplace-domain system
  `(iω·I - A) X̃ = u_τ - x_∞` where `A` is the Jacobian of the τ-equations
  (extracted once via symbolic finite differences), `u_τ` is the τ=0 state
  built from `u_end`, and `x_∞ = -A⁻¹ b_const` centres on the constant
  source. The vector overload `S(ω_vec, u_end, p0)` reuses one matrix
  extraction across all ω; the 4th/6th-cumulant spectrum agreement test
  in `test/higher_order_agreement_test.jl` matches master's 0.2 tolerance.
- `cavity_antiresonance_indexed.jl`, `waveguide.jl`, `unique_squeezing.jl`:
  the symbolic equation setup is intact; the heavier downstream features
  (`evaluate(; limits=...)` to unroll N to fixed sizes; `initial_values(eqs, ψ)`
  for spin coherent states; collective-decay rate-matrix solves) are slated
  for a follow-up minor release.
- `heterodyne_detection.jl`, `retrodiction_homodyne.jl` construct the
  forward/backward `NoiseMeanFieldEquations` via the unified `meanfield(...;
  efficiencies=..., direction=Forward()/Backward())` API. The downstream
  `SDESystem` integration via `System(eqs; name)` is built; the ensemble-averaging
  examples remain a follow-up.


<!-- Links generated by Changelog.jl -->

