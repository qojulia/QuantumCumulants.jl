# Changelog

All notable changes to QuantumCumulants.jl will be documented in this file.

## [1.0.0]

Breaking release: clean rewrite atop SecondQuantizedAlgebra v0.5 and ModelingToolkitBase. Migration is hand-driven; there are no deprecation shims.

### Fixed

- `scale` undercounted an `N` factor on cumulant-factorised products of averages
  derived from multi-index sums (`Σ_{i_1,…,i_k}`). After cumulant expansion the
  outer Σ scope (stored as `SumIndices`/`SumNonEqual` metadata on the original
  averaged term) was lost when `⟨A_{i_1}·B_{i_2}⟩` factorised to `⟨A_{i_1}⟩·⟨B_{i_2}⟩`,
  so `scale` saw a bare product with no sum information to translate into a
  range prefactor. `cumulant_expansion` now redistributes the metadata onto the
  resulting product by stamping each bound index onto the first averaged leaf
  whose op references it; bound indices that no leaf uses keep SQA's existing
  "spurious bound idx → factor 1" convention. The Ising-XX nonlinear-interaction
  agreement test in `test/indexed_scale_test.jl` (three assertions previously
  `@test_broken`) now passes.

- `to_system`'s conjugate substitution `⟨op†⟩ → conj(avg_op(t))` previously folded to `⟨op†⟩ → avg_op(t)` because the state variable carries `SymReal` symtype and Symbolics simplifies `conj(::Real)` to identity. This silently zeroed every `⟨X⟩ - ⟨X†⟩` driving term on the RHS of the compiled ODE, breaking driven cavities and any dissipative system with phase-sensitive coherent drive. Fix: build the conjugate via `SymbolicUtils.term(conj, var; type=Number)` instead of `Base.conj(var)`, bypassing the simplifier so the conj node survives through `mtkcompile` and `build_function`. ModelingToolkit's lack of a clean complex-state interface is tracked at https://github.com/SciML/ModelingToolkit.jl/issues/4548.

### Removed

- `simplify` keyword removed from every QC API that exposed it: `meanfield`, `complete!`, `complete`, `cumulant`, `cumulant_expansion`, `evaluate`, `translate_W_to_Y`, `CorrelationFunction`, `Spectrum`. QC never runs `SymbolicUtils.simplify` on its derived RHS expressions; canonical-form post-processing is the user's responsibility (`Symbolics.simplify(eq.rhs; expand=true)`). The default `simplify=true` was a 99% time sink at higher cumulant orders (e.g. order=6 `CorrelationFunction` went from ~100s to ~10ms with this change) for no downstream numerical benefit, and the `Spectrum` kwarg was unused dead code.
- `indexed_meanfield`, `indexed_complete!`, `IndexedCorrelationFunction`, `scaleME`, `evalME`, `IndexedMeanfieldEquations`, `EvaledMeanfieldEquations`, `ScaledMeanfieldEquations` collapsed into the unified `meanfield`, `complete!`, `CorrelationFunction`, `scale!`, `MeanFieldEquations`.
- `ClusterSpace` removed in SQA v0.5; use `Index` and `Σ` directly.
- `MeanfieldNoiseEquations`, `IndexedMeanfieldNoiseEquations`, `BackwardMeanfieldNoiseEquations` collapsed into `NoiseMeanFieldEquations{...,Direction}`.

### Changed

- `meanfield(ops, H, J; ...)` is the single derivation entry. Noise is opt-in via `efficiencies=...`; retrodiction via `direction=Backward()`.
- `scale(eqs; h=Int[])` and `evaluate(eqs; limits, h=Int[])` accept a `h::Vector{Int}` of `space_index` values selecting which Hilbert subspaces are collapsed / unrolled. The empty default targets every subspace (previous behaviour). Hybrid systems can unroll some subspaces with `evaluate` and collapse others with `scale` in any order.
- `scale` now drops additive constants from the LHS when SQA's commutator pipeline produces a `c + ⟨...⟩` shape (e.g. `⟨a*a'⟩` renamed to a same-index pair collapses to `1 + ⟨a'a⟩`). The constant-offset equation deduplicates against its normal-ordered sibling.
- Equation LHS in the user-facing struct is the raw `Average` BasicSymbolic. The `Differential(t)(u(t))` form is built only inside `to_system(...)`.
- Quality gates (Aqua + ExplicitImports + CheckConcreteStructs + JET) are part of CI.
- Updated to ModelingToolkitBase 1.36, SymbolicUtils 4, Symbolics 7, SecondQuantizedAlgebra 0.5.

### Migration

- Replace `using ModelingToolkit` with `using ModelingToolkitBase`.
- Replace `indexed_meanfield(...)` and `indexed_complete!(...)` with the unified `meanfield(...)` / `complete!(...)`.
- Replace `IndexedCorrelationFunction` with `CorrelationFunction`.
- Build numeric systems via `to_system(eqs; name=:sys)` then `mtkcompile(sys)`.
- Access trajectories via `get_solution(sol, op, eqs)(t)` instead of `sol[op]`.

### Examples

All 14 examples from the v0.4 series are ported to v1.0 and pass an automated
smoke run:

- `damped_cavity.jl` (new minimal demo with 12-digit analytical agreement)
- `ramsey_spectroscopy.jl`, `optomechanical-cooling.jl`, `many-atom-laser.jl`,
  `excitation-transport-chain.jl` — full meanfield → MTK → solve workflow
- `mollow.jl`, `single-atom-laser-spectrum.jl`, `superradiant_laser_indexed.jl`,
  `filter-cavity_indexed.jl` exercise the symbolic + numeric pipeline.
  `Spectrum(c, ps)(ω, u_end, p0)` solves the Laplace-domain system
  `(iω·I - A) X̃ = u_τ - x_∞` where `A` is the Jacobian of the τ-equations
  (extracted once via symbolic finite differences), `u_τ` is the τ=0 state
  built from `u_end`, and `x_∞ = -A⁻¹ b_const` centres on the constant
  source. The vector overload `S(ω_vec, u_end, p0)` reuses one matrix
  extraction across all ω; the 4th/6th-cumulant spectrum agreement test
  in `test/higher_order_agreement_test.jl` matches master's 0.2 tolerance.
- `cavity_antiresonance_indexed.jl`, `waveguide.jl`, `unique_squeezing.jl` —
  the symbolic equation setup is intact; the heavier downstream features
  (`evaluate(; limits=...)` to unroll N to fixed sizes; `initial_values(eqs, ψ)`
  for spin coherent states; collective-decay rate-matrix solves) are slated
  for a follow-up minor release.
- `heterodyne_detection.jl`, `retrodiction_homodyne.jl` — construct the
  forward/backward `NoiseMeanFieldEquations` via the unified `meanfield(...;
  efficiencies=..., direction=Forward()/Backward())` API. The downstream
  `SDESystem` integration via `to_system` is built; the ensemble-averaging
  examples remain a follow-up.
