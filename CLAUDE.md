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

**Always log `make test` output to a file** (e.g. `make test 2>&1 | tee /tmp/maketest.log`) so follow-up filters and triage queries grep the log instead of re-running the suite. The full run takes 90+ seconds; rerunning to extract one stack trace is wasteful. Same applies for `make jet` / `make format` if you need to inspect their output more than once.

## Project layout and dependencies

- The package is a thin symbolic layer on top of [SecondQuantizedAlgebra.jl](https://github.com/qojulia/SecondQuantizedAlgebra.jl) (SQA), which provides the operator algebra (`FockSpace`, `NLevelSpace`, `Destroy`, `Transition`, `Index`, `Σ`, …). `QuantumCumulants` `@reexport`s SQA, so users `using QuantumCumulants` get the full algebra surface.
- `Project.toml` pins SQA via `[sources]` to the `redesign-v2` branch on `https://github.com/qojulia/SecondQuantizedAlgebra.jl`. When debugging issues that look algebraic (commutators, indexed sums, operator ordering), check whether the bug actually lives in SQA, not here.
- MTK conversion uses `ModelingToolkitBase` (aliased `const MTK`), not full `ModelingToolkit`. Examples and docs say `using ModelingToolkitBase`.

## Architecture (big picture)

Entry point: [src/QuantumCumulants.jl](src/QuantumCumulants.jl). The pipeline is symbolic derivation → cumulant expansion → completion → (optional) scaling/noise → MTK system → numeric solve.

- [equations.jl](src/equations.jl) — `AbstractMeanfieldEquations`, `MeanfieldEquations`, `NoiseMeanfieldEquations{...,Direction}`, `EvolutionDirection` (`Forward`/`Backward`). LHS in the user-facing struct is the raw `Average` `BasicSymbolic`; the `Differential(t)(u(t))` form is built only inside `MTK.System(eqs; name)`.
- [meanfield.jl](src/meanfield.jl) — unified `meanfield(ops, H, J; rates, order, efficiencies, direction)`. There is no longer a separate `indexed_meanfield`; indexed and scalar systems share one entry point. Noise is opt-in via `efficiencies=...`; retrodiction via `direction=Backward()`.
- [cumulant.jl](src/cumulant.jl) — `cumulant_expansion`, `cumulant`, `get_order` (Wick-style factorization to a given order).
- [completion.jl](src/completion.jl) — `complete`/`complete!`/`find_missing` close the system by deriving equations for any averages appearing on the RHS but missing from the LHS.
- [scaling.jl](src/scaling.jl) — `scale`/`scale!` for permutation-symmetric (indexed) systems. Defined for both `MeanfieldEquations` and `NoiseMeanfieldEquations`.
- [mtk.jl](src/mtk.jl) — `MTK.System(eqs; name)` is QC's extension of `ModelingToolkitBase.System`, dispatched on `MeanfieldEquations`/`NoiseMeanfieldEquations`/`CorrelationFunction`; `initial_values(eqs, u0)`, `get_solution(sol, op, eqs)(t)` are the user-facing bridges to numerical solutions. `sol[op]` indexing is gone; always go through `get_solution`.
- [correlation.jl](src/correlation.jl) — `CorrelationFunction`, `Spectrum`, `correlation_u0`, `correlation_p0` (unified, no separate `IndexedCorrelationFunction`).
- [noise.jl](src/noise.jl) — diffusion/noise machinery for `NoiseMeanfieldEquations`.
- [latexify.jl](src/latexify.jl) — `@latexrecipe`s for pretty-printing equations.

## v1 rewrite context

This is the `rewrite` branch building the breaking 1.0 release on top of SQA v0.5 + ModelingToolkitBase. There are no deprecation shims; migration from `master` is hand-driven. Key collapses (see [CHANGELOG.md](CHANGELOG.md)):

- `indexed_meanfield`, `indexed_complete!`, `IndexedCorrelationFunction`, `scaleME`/`evalME`, `IndexedMeanfieldEquations`, `EvaledMeanfieldEquations`, `ScaledMeanfieldEquations` → unified `meanfield` / `complete!` / `CorrelationFunction` / `scale!` / `MeanfieldEquations`.
- `MeanfieldNoiseEquations`, `IndexedMeanfieldNoiseEquations`, `BackwardMeanfieldNoiseEquations` → `NoiseMeanfieldEquations{...,Direction}`.
- `ClusterSpace` removed (SQA v0.5); use `Index` + `Σ` directly.

When in doubt about whether something is "missing" vs. "moved," check [CHANGELOG.md](CHANGELOG.md) Migration section.

## What the v1 rewrite is *for*

The point of the rewrite is **not** to reproduce master line-by-line on top of new dependency versions. The point is a **maintainable, scalable, performant** implementation that uses SQA's new data model natively. The starting reading list is SQA's [docs/src/devdocs.md](https://github.com/qojulia/SecondQuantizedAlgebra.jl/docs/src/devdocs.md) (architecture, naming policy, canonicalisation pipeline, diagonal splitting) and [docs/src/symbolic_sums.md](https://github.com/qojulia/SecondQuantizedAlgebra.jl/docs/src/symbolic_sums.md). Read these *first* before porting anything from master; the answers to "what should this function do" come from there, not from master's source.

### Key SQA invariants QC code must respect

- **Naming policy (operational).** `Index` objects are *user-owned*. SQA never mints `Index` on the user's behalf because an algebra-invented index breaks downstream pattern-matching, initial-condition substitution, and follow-on `evaluate(...; limits=...)` calls (devdocs *Naming policy*). When QC needs new indices (e.g. inside `scale!`, `evaluate`, `complete!`'s alpha canonicalisation), it must derive names from the user's existing vocabulary (typically the first-declared index per Hilbert subspace) so the result stays traceable and usable.
- **`*(QAdd, QAdd)` is strict.** Two `QAdd`s with overlapping `.indices` raise `ArgumentError` instead of silently alpha-renaming (devdocs *Disjoint bound indices in products*). Multiply only after a `SQA.change_index` rename on one side.
- **Undetermined free-index pairs stay put.** Two ops on the same subspace with different free symbolic indices and no `ne` are `Undetermined`; `_partial_sort!` deliberately leaves them in physical order, same-site collapse does not fire (devdocs *Operator sorting*, *Free indices and `assume_distinct_index`*). Use `assume_distinct_index(q, [(i, j)])` when QC code knows the user means "distinct atoms"; do not invent a private NE convention.
- **Sum metadata round-trip.** `average(QAdd_with_indices)` stamps each averaged term with `SumIndices`/`SumNonEqual` metadata; `undo_average` restores it. Comparing or canonicalising averages must not silently strip the metadata, or constant-operator sums (e.g. `⟨1⟩` inside `Σ_i`) collapse incorrectly (devdocs *Averaging*).
- **`change_index` is the only knob.** Substituting an `Index` is the single primitive that carries algebraic normalisation: `_canonicalize!` for non-sum results (projector squashing, commuting reorder, ne-violated drop), `_accumulate_with_diag!` for sum results (off-diagonal + diagonal split under existing `.indices`), `_substitute_ne` for NE propagation. If a QC pass touches indices, it should be doing it through `change_index`, not by hand.
- **One compound type.** Everything is `QSym` (atoms) or `QAdd` (sums of products with prefactors). No `QMul`. `QTerm` is a dict key, not a value type. Don't write code that branches on `QMul` or assumes a 3-level expression tree (devdocs *Type hierarchy*).
- **Ground-state projectors stay atomic.** `σᵍᵍ` is a legitimate canonical-form atom, not auto-expanded into `1 - Σ σᵏᵏ`. Use `expand_completeness` only when the basis change is the explicit goal. QC's mean-field code already relies on this — don't accidentally trigger the expansion (devdocs *Simplification vs. normal ordering*).

### Concrete coding-style rules

- **Read SQA's source first**, not master's, when porting a feature. The right port is "what do SQA primitives imply this should look like," not "what would master's lines translate to."
- **Reach for `SQA.change_index`, `SQA.Σ`, `*`, `+` first.** Manual `QAdd(QTermDict(QTerm(ops, ne) => c), idx)` construction is a last resort and must be justified (typically only inside SQA itself).
- **Mint Indices, when required, from the user's vocabulary.** Derive names from declared indices (`Symbol(i.name, "_", k)`); do not invent fresh prefixes like `:_eval_1`. Per-Hilbert-space canonicalisation lives in [completion.jl](src/completion.jl)'s `_build_canonical_indices`.
- **`_canonical_key` in completion.jl is for bound-index alpha-equivalence**, not for distinguishing concrete positions. After `evaluate`, the indices are *concrete* (the user constructed them, or we minted them from the user's vocabulary), and the dedup key should be the operator itself, optionally normalised for commuting-op ordering — not `_canonical_key`.
- **Don't keep parallel copies of master helpers.** If SQA is missing a primitive QC genuinely needs, add it to SQA (or write a one-liner here that composes existing ones). Do not duplicate.
- **Failure mode to watch for**: building a `QAdd` by hand and then being surprised that `σ_ee² ≠ σ_ee`, or that `σ_2 σ_1 ≠ σ_1 σ_2`. The fix is to route through SQA's product (`*`) or `change_index`, which run the canonicalisation pipeline; never to hand-roll a sort.
- **Use `Symbolics.IM`, not Julia's `im`, in symbolic expressions.** The meanfield codepath emits the symbolic imaginary unit; Julia's `im` becomes a `complex(0, …)` literal that SymbolicUtils does not unify with the factored form, so `simplify(rhs - expected)` will not cancel even when the two are numerically identical. Tests asserting symbolic equality of RHS expressions must write `2 * Symbolics.IM * Δ * ⟨…⟩`, not `2im * Δ * ⟨…⟩`.

## Test layout

- [test/runtests.jl](test/runtests.jl) — autodiscovers via `ParallelTestRunner.find_tests`.
- [test/systems/](test/systems/) — end-to-end physical models (cavity, Dicke, two-level laser, spectrum, multilevel) used as regression scenarios.
- [test/quality/](test/quality/) — `quality_test.jl` (Aqua + ExplicitImports + CheckConcreteStructs) and `JET.jl` (run separately via `make jet`; skipped from the default suite when filtered).
- Operator-algebra/index/sum unit tests belong in **SecondQuantizedAlgebra.jl**, not here. Tests in this repo should exercise the QC surface (mean-field derivation, cumulant expansion, completion, scaling, MTK bridge, correlation/spectrum), not low-level SQA behavior.
