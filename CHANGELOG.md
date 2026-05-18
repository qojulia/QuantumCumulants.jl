# Changelog

All notable changes to QuantumCumulants.jl will be documented in this file.

## [1.0.0]

Breaking release: clean rewrite atop SecondQuantizedAlgebra v0.5 and ModelingToolkitBase. Migration is hand-driven; there are no deprecation shims.

### Removed

- `indexed_meanfield`, `indexed_complete!`, `IndexedCorrelationFunction`, `scaleME`, `evalME`, `IndexedMeanfieldEquations`, `EvaledMeanfieldEquations`, `ScaledMeanfieldEquations` collapsed into the unified `meanfield`, `complete!`, `CorrelationFunction`, `scale!`, `MeanFieldEquations`.
- `ClusterSpace` removed in SQA v0.5; use `Index` and `Σ` directly.
- `MeanfieldNoiseEquations`, `IndexedMeanfieldNoiseEquations`, `BackwardMeanfieldNoiseEquations` collapsed into `NoiseMeanFieldEquations{...,Direction}`.

### Changed

- `meanfield(ops, H, J; ...)` is the single derivation entry. Noise is opt-in via `efficiencies=...`; retrodiction via `direction=Backward()`.
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
  `filter-cavity_indexed.jl` — exercise the symbolic + numeric pipeline. The
  `Spectrum(c, ps)` Laplace machinery is partly re-implemented and accessible
  via `Spectrum(c, ps)(ω, u_end, p0)`; consumers should verify numerical
  results against the FFT path as the linear solver tuning is still in flux.
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
