# Changelog

All notable changes to QuantumCumulants.jl will be documented in this file.

## [0.5.2]

### Fixed

`complete(...; get_adjoints = false)` now keeps the conjugation-canonical member of each conjugate pair, instead of whichever side the right-hand-side traversal happened to reach first. That traversal order is `objectid`-seeded and so varies between Julia processes, which made the surviving representative (and therefore the equation count after `scale`) nondeterministic across sessions: the heterodyne-detection system, for instance, scaled to 12 or 13 equations depending on the process. The closed set is now process-independent. See issue [#295](https://github.com/qojulia/QuantumCumulants.jl/issues/295).

## [0.5.0]

Breaking release: a clean rewrite of QuantumCumulants.jl atop **SecondQuantizedAlgebra.jl v0.6** (which now owns the operator algebra) and **ModelingToolkitBase.jl** (in place of the full ModelingToolkit). Migration from the v0.4 series is hand-driven; there are no deprecation shims. The substantive changes fall into three groups: direct renames, constructs replaced by more general machinery, and behavioural changes in shared API names. All 14 examples from the v0.4 series are ported and pass an automated smoke run. This entry is the migration reference for users with code written against **QuantumCumulants.jl v0.4**. Algebra-surface migration (Hilbert spaces, operators, indices, sums) lives in the [SecondQuantizedAlgebra.jl changelog](https://qojulia.github.io/SecondQuantizedAlgebra.jl/stable/changelog/).

### Added

**Graph-native core (`MomentGraph`).** The cumulant hierarchy is now carried by a single `MomentGraph`: nodes are canonicalised moments, edges are the drift couplings between them, and per-subspace `treatments` record whether a subspace is left free, unrolled, or scaled. `meanfield` seeds and closes it, and `complete`/`complete!`, `scale`/`scale!`, and `evaluate` transform it; `System` reads from it to emit code. Every `*Equations` struct now stores `graph` as the source of truth, with `equations` a derived view regenerated from the graph rather than a hand-maintained parallel list. Correlation functions and spectra are derived on the same graph.

**Introspection helpers.** New exported accessors over the graph: `states`, `moments`, `moment_variable_map` (moment to `BasicSymbolic` variable), `closure_report` (which moments were truncated and how), `noise_channels` (per-jump noise drift on a `NoiseMeanfieldEquations`), and `parameter_map`.

**Direction types.** `EvolutionDirection` with `Forward` and `Backward` are exported, replacing the separate `meanfield_backward` entry point (`meanfield(...; direction=Backward())`).

**In-place variants.** `scale!` and `simplify!` join the existing `complete!` / `modify_equations!`.

**Correlation functions of time-dependent Hamiltonians.** `CorrelationFunction` now supports systems whose Hamiltonian, jumps, or rates depend on the time variable (e.g. a drive `f(t)`); previously these failed with `UndefVarError: t not defined`. See issues [#171](https://github.com/qojulia/QuantumCumulants.jl/issues/171) and [#93](https://github.com/qojulia/QuantumCumulants.jl/issues/93). Per the quantum regression theorem the τ-evolution is governed by the generator at `t₀+τ`, so the parent time variable is substituted by `iv0 + τ`. Pass the parent evolution's stop time `t₀` via the new `iv0` keyword: `@variables t0::Real; CorrelationFunction(op1, op2, eqs; iv0 = t0)`, and supply its value in `correlation_p0`. Time-independent systems are unaffected (the keyword is ignored when no time dependence is present).

### Renamed

One-to-one renames. Replace the left column with the right column anywhere it appears in your code.

| v0.4 | v0.5 |
|---|---|
| `indexed_meanfield(...)` | `meanfield(...)` |
| `indexed_complete(...)` / `indexed_complete!(...)` | `complete(...)` / `complete!(...)` |
| `find_missing_sums(...)` | `find_missing(...)` |
| `scaleME(...)` | `scale(...)` / `scale!(...)` |
| `evalME(...)` | `evaluate(...)` |
| `IndexedCorrelationFunction(...)` | `CorrelationFunction(...)` |
| `heisenberg(...)` | `meanfield(...)` |
| `meanfield_backward(...)` | `meanfield(...; direction=Backward())` |
| `IndexedMeanfieldEquations` / `EvaledMeanfieldEquations` / `ScaledMeanfieldEquations` | `MeanfieldEquations` |
| `MeanfieldNoiseEquations` / `IndexedMeanfieldNoiseEquations` / `BackwardMeanfieldNoiseEquations` | `NoiseMeanfieldEquations{…,Direction}` |
| `to_system(eqs)` | `System(eqs; name)` |
| `sol[op]` / `get_scale_solution(...)` | `get_solution(sol, op, eqs)(t)` |
| `@cnumbers a b c` | `@variables a b c` *(reexported from Symbolics)* |
| `@rnumbers a b c` | `@variables a::Real b::Real c::Real` |
| `Parameter(:x)` | `@variables x`|

### Removed

A handful of constructs that lived at the type level in v0.4 have been replaced by mechanisms that scale better and avoid bespoke types.

**`ClusterSpace`.** Indexed families of identical subsystems are now expressed through the SQA `Index` and `Σ` machinery. Where v0.4 wrote

```julia
hc = ClusterSpace(NLevelSpace(:atom, 2), N, k)
σ(i, j) = Transition(hc, :σ, i, j)
H = Δ * sum(σ(2, 2))
```

v0.5 writes

```julia
ha = NLevelSpace(:atom, 2)
i  = Index(ha, :i, N, ha)
σ(α, β, k) = IndexedOperator(Transition(ha, :σ, α, β, 2), k)
H = Δ * Σ(σ(2, 2, i), i)
```

`N` stays symbolic through equation derivation. The cumulant truncation order, which used to be the `k` parameter on `ClusterSpace`, is now passed to the derivation: `meanfield(…; order=k)`.

**The indexed-sum type zoo.** `IndexedAverageSum`, `IndexedAverageDoubleSum`, `SingleSum`, `DoubleSum`, `AvgSums`, `NumberedOperator`, and `SpecialIndexedTerm` are gone. They were SQA-internal in v0.4; SQA v0.5 represents indexed sums uniformly as a `QAdd` with `.indices` metadata, so the type zoo collapsed. Code that pattern-matched on these should iterate the resulting `QAdd` and inspect `term.ne` directly.

**Indexed-helper exports.** `value_map`, `split_sums`, `scale_term`, `eval_term`, `insert_index`, and `order_by_index` are removed. The equivalent transformations are expressed through SQA's `change_index` and `Σ` directly.

**Parameter machinery.** `Parameter`, `CNumber`, `RNumber`, `RealParameter`, the `Average` *type*, and the `@cnumbers`/`@rnumbers`/`cnumbers`/`cnumber`/`rnumbers`/`rnumber` macros are removed. Use `@variables x::Real` (or `::Complex`) for symbolic parameters.

**The `simplify` keyword.** Removed from every QC entry point that exposed it (`meanfield`, `complete`/`complete!`, `cumulant`, `cumulant_expansion`, `evaluate`, `translate_W_to_Y`, `CorrelationFunction`, `Spectrum`). QC no longer runs `SymbolicUtils.simplify` on its derived right-hand sides; canonical-form post-processing is the caller's responsibility (`Symbolics.simplify(eq.rhs; expand=true)`). The old default `simplify=true` was a 99% time sink at higher cumulant orders (e.g. an order-6 `CorrelationFunction` dropped from ~100s to ~10ms) for no downstream numerical benefit.

**`subst_reds`.** Redundant-conjugate substitution is folded into the `System(eqs)` codegen path automatically.

### Changed

A few API names survive but their behaviour has shifted. The differences are intentional and pay off in correctness or performance, so be aware of them when porting tests.

**`meanfield` is the single derivation entry.** `meanfield(ops, H, J; …)` replaces the scalar/indexed/noise/backward split. Noise is opt-in via `efficiencies=…` (yielding a `NoiseMeanfieldEquations`); retrodiction via `direction=Backward()`.

**`scale` and `evaluate` take a subspace selector.** `scale(eqs; h=Int[])` and `evaluate(eqs; limits, h=Int[])` accept a `h::Vector{Int}` of `acts_on` indices selecting which Hilbert subspaces are collapsed / unrolled. The empty default targets every subspace (the previous behaviour). Hybrid systems can unroll some subspaces with `evaluate` and collapse others with `scale`, in any order.

**`scale` drops additive constants from the LHS.** When SQA's commutator pipeline produces a `c + ⟨…⟩` shape (e.g. `⟨a a'⟩` collapsing to `1 + ⟨a'a⟩`), the constant-offset equation deduplicates against its normal-ordered sibling.

**Equation LHS is the raw `Average`.** In the user-facing struct, the LHS is the bare `Average` `BasicSymbolic`. The `Differential(t)(u(t))` form is built only inside `ModelingToolkitBase.System(eqs; name)`, which QC extends to dispatch on `MeanfieldEquations`, `NoiseMeanfieldEquations`, and `CorrelationFunction`.

**`System` conjugate substitution survives compilation.** `System(eqs)` builds `conj(avg_op(t))` via `SymbolicUtils.term(conj, var; type=Number)` rather than `Base.conj(var)`, so the symbolic `conj` node survives `mtkcompile` and `build_function`. Without this, every `⟨X⟩ - ⟨X†⟩` driving term silently zeroed on the compiled RHS because Symbolics simplifies `conj(::SymReal)` to identity. Related ModelingToolkit gap tracked at <https://github.com/SciML/ModelingToolkit.jl/issues/4548>.

**`cumulant_expansion` redistributes sum-scope metadata.** The `SumIndices`/`SumNonEqual` metadata is pushed onto each factorised product term so `scale` can recover the per-index range prefactor after Wick factorisation. Bound indices that no factored leaf references keep SQA's "spurious bound idx, factor 1" convention.

**Dependencies and CI.** Updated to ModelingToolkitBase 1.36, SymbolicUtils 4, Symbolics 7, SecondQuantizedAlgebra 0.6.3. Quality gates (Aqua + ExplicitImports + CheckConcreteStructs + JET) are now part of CI.

### Migration

Side-by-side, the single-atom-laser workflow from the tutorial.

```julia
# v0.4
using QuantumCumulants, ModelingToolkit
@cnumbers Δ g γ κ ν
hf = FockSpace(:cavity); ha = NLevelSpace(:atom, (:g, :e)); h = hf ⊗ ha
@qnumbers a::Destroy(h)
σ(i, j) = Transition(h, :σ, i, j)
H = Δ*a'*a + g*(a'*σ(:g, :e) + a*σ(:e, :g))
J = [a, σ(:g, :e), σ(:e, :g)]
eqs  = meanfield([a'*a, σ(:e, :e), a'*σ(:g, :e)], H, J; rates=[κ, γ, ν], order=2)
eqs  = complete(eqs)
sys  = to_system(eqs)
prob = ODEProblem(sys, u0, (0.0, 10.0), p0)
sol  = solve(prob, Tsit5())
n    = sol[a'*a]

# v0.5
using QuantumCumulants, ModelingToolkitBase
@variables Δ::Real g::Real γ::Real κ::Real ν::Real
hf = FockSpace(:cavity); ha = NLevelSpace(:atom, (:g, :e)); h = hf ⊗ ha
@qnumbers a::Destroy(h)
σ(i, j) = Transition(h, :σ, i, j)
H = Δ*a'*a + g*(a'*σ(:g, :e) + a*σ(:e, :g))
J = [a, σ(:g, :e), σ(:e, :g)]
eqs  = meanfield([a'*a, σ(:e, :e), a'*σ(:g, :e)], H, J; rates=[κ, γ, ν], order=2)
eqs  = complete(eqs)
sys  = mtkcompile(System(eqs; name=:laser))
prob = ODEProblem(sys, merge(initial_values(eqs, u0), Dict(p0)), (0.0, 10.0))
sol  = solve(prob, Tsit5())
n    = real.(get_solution(sol, a'*a, eqs).(ts))
```

The differences are local: `@cnumbers` becomes `@variables`, `using ModelingToolkit` becomes `using ModelingToolkitBase`, `to_system(eqs)` becomes `mtkcompile(System(eqs; name))`, and `sol[op]` becomes `get_solution(sol, op, eqs)(t)`. For numeric initial values from *symbolic* level names, pre-translate the `:g`/`:e` labels to integer levels before calling `initial_values(eqs, ψ)` (the v0.4 `level_map` keyword was dropped). Everything else (`meanfield`, `complete`, the Hamiltonian) keeps the same names and meaning.

### Unchanged

These QC names keep their meaning across the migration; code that only uses them ports without changes:

`meanfield`, `complete`/`complete!`, `find_missing`, `find_operators`, `cumulant`, `cumulant_expansion`, `get_order`, `average`, `undo_average`, `CorrelationFunction`, `Spectrum`, `correlation_u0`, `correlation_p0`, `scale`/`scale!`, `evaluate`, `initial_values`.

The entire SecondQuantizedAlgebra algebra surface (Hilbert spaces, operators, indices, sums, `to_numeric`, `numeric_average`) is reexported unchanged; see the [SecondQuantizedAlgebra.jl migration notes](https://qojulia.github.io/SecondQuantizedAlgebra.jl/stable/changelog/) for that layer.


<!-- Links generated by Changelog.jl -->

