# CLAUDE.md

## Common commands

All workflows go through the `Makefile`. Drive the whole test suite via `make test`, not bare `julia -e 'include(...)'` on individual files (per repo convention).

- `make test` — full suite. Uses `ParallelTestRunner` (autodiscovery from [test/](test/)).
- `make jet` — JET static analysis ([test/quality/JET.jl](test/quality/JET.jl)).
- `make format` — Runic formatter over `src/ test/ benchmark/ examples/ docs/`.
- `make docs` / `make servedocs` — build / live-serve documentation.
- `make bench` / `make benchlocal` — benchmarks under [benchmark/](benchmark/).
- `make setup` — installs `Changelog`, `LiveServer`, and the `Runic` app once.

Running a single test file: pass a name filter to `ParallelTestRunner` via `Pkg.test` args, e.g.
`julia --project -e 'using Pkg; Pkg.test(test_args=["completion"])'`.
The runner skips anything under [test/pending/](test/pending/) (ports blocked on v1 regressions, see [TODO.md](TODO.md)).

## Project layout and dependencies

- The package is a thin symbolic layer on top of [SecondQuantizedAlgebra.jl](https://github.com/qojulia/SecondQuantizedAlgebra.jl) (SQA), which provides the operator algebra (`FockSpace`, `NLevelSpace`, `Destroy`, `Transition`, `Index`, `Σ`, …). `QuantumCumulants` `@reexport`s SQA, so users `using QuantumCumulants` get the full algebra surface.
- `Project.toml` pins SQA via a local `[sources]` path (`/home/oameye/Documents/SecondQuantizedAlgebra.jl`). When debugging issues that look algebraic (commutators, indexed sums, operator ordering), check whether the bug actually lives in SQA, not here.
- MTK conversion uses `ModelingToolkitBase` (aliased `const MTK`), not full `ModelingToolkit`. Examples and docs say `using ModelingToolkitBase`.

## Architecture (big picture)

Entry point: [src/QuantumCumulants.jl](src/QuantumCumulants.jl). The pipeline is symbolic derivation → cumulant expansion → completion → (optional) scaling/noise → MTK system → numeric solve.

- [equations.jl](src/equations.jl) — `AbstractMeanFieldEquations`, `MeanFieldEquations`, `NoiseMeanFieldEquations{...,Direction}`, `EvolutionDirection` (`Forward`/`Backward`). LHS in the user-facing struct is the raw `Average` `BasicSymbolic`; the `Differential(t)(u(t))` form is built only inside `to_system`.
- [meanfield.jl](src/meanfield.jl) — unified `meanfield(ops, H, J; rates, order, efficiencies, direction)`. There is no longer a separate `indexed_meanfield`; indexed and scalar systems share one entry point. Noise is opt-in via `efficiencies=...`; retrodiction via `direction=Backward()`.
- [cumulant.jl](src/cumulant.jl) — `cumulant_expansion`, `cumulant`, `get_order` (Wick-style factorization to a given order).
- [completion.jl](src/completion.jl) — `complete`/`complete!`/`find_missing` close the system by deriving equations for any averages appearing on the RHS but missing from the LHS.
- [scaling.jl](src/scaling.jl) — `scale`/`scale!` for permutation-symmetric (indexed) systems. Currently only defined for `MeanFieldEquations`; scaling `NoiseMeanFieldEquations` is a known gap (see [TODO.md](TODO.md)).
- [mtk.jl](src/mtk.jl) — `to_system(eqs; name)` builds the MTK `System`; `initial_values(eqs, u0)`, `get_solution(sol, op, eqs)(t)` are the user-facing bridges to numerical solutions. `sol[op]` indexing is gone; always go through `get_solution`.
- [correlation.jl](src/correlation.jl) — `CorrelationFunction`, `Spectrum`, `correlation_u0`, `correlation_p0` (unified, no separate `IndexedCorrelationFunction`).
- [noise.jl](src/noise.jl) — diffusion/noise machinery for `NoiseMeanFieldEquations`.
- [latexify.jl](src/latexify.jl) — `@latexrecipe`s for pretty-printing equations.

## v1 rewrite context

This is the `rewrite` branch building the breaking 1.0 release on top of SQA v0.5 + ModelingToolkitBase. There are no deprecation shims; migration from `master` is hand-driven. Key collapses (see [CHANGELOG.md](CHANGELOG.md)):

- `indexed_meanfield`, `indexed_complete!`, `IndexedCorrelationFunction`, `scaleME`/`evalME`, `IndexedMeanfieldEquations`, `EvaledMeanfieldEquations`, `ScaledMeanfieldEquations` → unified `meanfield` / `complete!` / `CorrelationFunction` / `scale!` / `MeanFieldEquations`.
- `MeanfieldNoiseEquations`, `IndexedMeanfieldNoiseEquations`, `BackwardMeanfieldNoiseEquations` → `NoiseMeanFieldEquations{...,Direction}`.
- `ClusterSpace` removed (SQA v0.5); use `Index` + `Σ` directly.

When in doubt about whether something is "missing" vs. "moved," check [CHANGELOG.md](CHANGELOG.md) Migration section and [TODO.md](TODO.md) (lists known v1 regressions still to port: `non_equal=true` kwarg, `scale(::NoiseMeanFieldEquations)`, `evaluate(eqs; limits=...)`, `meanfield_backward`, rate-matrix collective decay, etc.).

## Test layout

- [test/runtests.jl](test/runtests.jl) — autodiscovers via `ParallelTestRunner.find_tests` and excludes [test/pending/](test/pending/).
- [test/systems/](test/systems/) — end-to-end physical models (cavity, Dicke, two-level laser, spectrum, multilevel) used as regression scenarios.
- [test/quality/](test/quality/) — `quality_test.jl` (Aqua + ExplicitImports + CheckConcreteStructs) and `JET.jl` (run separately via `make jet`; skipped from the default suite when filtered).
- [test/pending/](test/pending/) — kept in-tree as a porting source-of-truth; do not enable until the corresponding v1 feature lands.
- Operator-algebra/index/sum unit tests belong in **SecondQuantizedAlgebra.jl**, not here. Tests in this repo should exercise the QC surface (mean-field derivation, cumulant expansion, completion, scaling, MTK bridge, correlation/spectrum), not low-level SQA behavior.
