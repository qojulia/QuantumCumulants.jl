# Changelog

All notable changes to QuantumCumulants.jl will be documented in this file.

## [0.6.0]

### Added

Direct RHS compilation for the completed moment equations (issue [#294](https://github.com/qojulia/QuantumCumulants.jl/issues/294)): `ODEProblem(eqs, u0, tspan, ps)` and `ODEFunction(eqs, ps)` now build a solver-ready problem straight from a `MeanfieldEquations`, without the `System`/`mtkcompile` detour. Two official backends, selected by the `backend` keyword:

- `KernelBackend()`: the moment-polynomial kernel `du = M * v`. The drifts are lowered once into sparse data (prefix-chained monomials plus one sparse coefficient matrix, stored transposed in CSR layout so the `M * v` is a cache-friendly per-row gather, several times faster than a CSC `mul!` on large systems), so there is no per-model native code and cold construction is fast. `parallel = :auto` additionally threads both RHS passes (the monomial update and the gather) through Polyester's persistent task pool on systems of `≥ 64` equations, for a 2-3x further speedup (overridable; serial and threaded are bit-identical). The kernel keeps one scratch buffer per thread, so a single `prob` is reentrant: an `EnsembleProblem` with `EnsembleThreads()` that shares the problem and only varies `u0` needs no `copy` (per-trajectory parameters still do, to own an independent `Mt`). Supports `jac = true` for an analytic sparse Jacobian (holomorphic systems; conj-folded closures raise `HolomorphicJacobianError`) and `KernelBackend(cache = path)` for persisting the lowered tables in a JLD2 file (digest-keyed, verified byte-exact on load; requires `using JLD2`, provided by a package extension).
- `ShardedBackend()`: chunked native codegen (per-chunk `build_function` compiled as RuntimeGeneratedFunctions behind a FunctionWrapper table). The Expr generation and native-code warm-up are parallelized with OhMyThreads schedulers (greedy load-balancing for the heterogeneous per-chunk compile cost). Handles t-dependent and non-polynomial drifts. `chunk = :auto` sizes the chunks to the thread count, and `parallel = :auto` runs the chunks concurrently at solve time (`Threads.@threads :dynamic` over disjoint `du` ranges) for large systems, roughly halving the RHS time above ~1000 equations on a multi-threaded session; both accept explicit overrides. Supports `jac = true`/`:fd` (a colored finite-difference Jacobian driven through the chunked RHS: sparsity read from the drift expressions, greedily colored, filled in `ncolors + 1` RHS evaluations) and `jac = :analytic` (a symbolic Jacobian codegen'd like the RHS); a finite-difference time-gradient is attached too, so implicit solvers such as `Rodas5P()` run on the complex state with no `autodiff` keyword. Both compute the real-directional linearization, so unlike the kernel's holomorphic Jacobian they also handle conjugate-folded (`get_adjoints = false`) closures.
- `AutoBackend()` (the default) uses the kernel and falls back to sharded exactly when the drift is not polynomial in the moments.

Parameters are required at construction; `update_parameters!(prob, Dict(...))` rewrites the compiled RHS in place for sweeps (array-valued parameters included), and `copy(prob.f)` gives ensembles isolated tables. Both backends resolve drifts through the same treatments-aware machinery as `System(eqs)`, so `scale`d and `evaluate`d systems work identically. Lowering failures are typed (`NonPolynomialDriftError`, `TimeDependentCoefficientError`, `ImParameterCollisionError`, `UnresolvedMomentError`) with guidance in the message. Kernel lowering parallelizes its dominant `polynomial_coeffs` pass with an OhMyThreads greedy scheduler, and a PrecompileTools workload covers the meanfield-to-first-RHS-call path. `CorrelationFunction`, noise/SDE systems, events, and symbolic solution indexing stay on the ModelingToolkit path; see the new "Solving the equations directly (RHS backends)" documentation page.

## [0.5.5]

### Changed

The Lindblad drift builders now collect their per-jump contributions into a vector and sum them in a single pass, instead of folding them in with a running `acc += …`. Every `QAdd + QAdd` rebuilds the accumulator's term dictionary, so the running fold cost O(n²) in the number of jump operators; SQA's `sum(::AbstractArray{<:QAdd})` routes the whole list through one MutableArithmetics accumulator (`_QAddBuilder`), making drift assembly O(n). This affects the forward and backward master-equation drifts (`_lindblad_rhs`, `_master_lindblad_backward`), the forward noise-equation builder, and the ancilla-undo step used by `CorrelationFunction`.

Move to SecondQuantizedAlgebra v0.9

## [0.5.4]

### Changed

Upgraded to **SecondQuantizedAlgebra.jl v0.8** on its public API. No user-facing QuantumCumulants API changed. See the [SecondQuantizedAlgebra.jl changelog](https://qojulia.github.io/SecondQuantizedAlgebra.jl/stable/changelog/).

### Fixed

The correlation conjugate-fold is now collapse-aware. The identity `⟨A†⟩ = ⟨A⟩*` is unsound for the two-time quantum-regression-theorem τ-states, where the ancilla collapse adds a boson commutator (`⟨a'aσ⟩` versus `⟨aa'σ⟩`), so folding could drive the computed power spectrum negative. The fold is now gated by a collapse-aware `foldable` predicate threaded through `MomentMap`, the closure, the state registry and the spectrum; single-time systems keep folding every conjugate pair and are unchanged.

## [0.5.3]

### Changed

Migrated to **SecondQuantizedAlgebra.jl v0.7** on its public API. SQA v0.7 represents indexed sums as a `Term` (rather than the previous `SymbolicUtils` node shape) and exposes the algebra surface QuantumCumulants relies on (tree rewriting, averaging, scaling, cumulant expansion, coefficient handling) through documented entry points, so QuantumCumulants no longer reaches into SQA internals. No user-facing QuantumCumulants API changed. See the [SecondQuantizedAlgebra.jl changelog](https://qojulia.github.io/SecondQuantizedAlgebra.jl/stable/changelog/).

## [0.5.2]

### Fixed

`scale` now folds the conjugate images of a conjugation-reduced system (`complete(...; get_adjoints = false)`), so its equation count is the minimal closed set rather than depending on which member of each conjugate pair the closure happened to keep. `scale` previously reduced only permutation symmetry, leaving moments that coincide only after conjugation as separate equations for a single numerical unknown; which of those survived was `objectid`-seeded and so varied between Julia processes, making the post-`scale` count nondeterministic across sessions (the heterodyne-detection system scaled to 12 or 13 equations depending on the process). `scale` now reduces permutation *and* conjugation symmetry consistently with the generated numerical system, so the count is process-independent. See issue [#295](https://github.com/qojulia/QuantumCumulants.jl/issues/295).

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

