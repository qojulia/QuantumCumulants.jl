# Changelog

All notable changes to QuantumCumulants.jl will be documented in this file.

## [0.7.0]

### Changed

Moved to SecondQuantizedAlgebra v0.10. Load a numeric backend yourself (`using QuantumOpticsBase` or `using QuantumToolbox`), and Hermitian-operator averages now carry the `Real` symtype.
QuantumCumulants now relies only on the public SecondQuantizedAlgebra API. See issue [#300](https://github.com/qojulia/QuantumCumulants.jl/issues/300).

**Breaking.** `complete!` and `complete` now default to `get_adjoints = false`, keeping one representative per conjugate pair instead of a separate state for every conjugate moment. The partner is recovered by `conj` at code generation, so the dynamics are unchanged while the closed system is smaller (the optomechanical-cooling example closes at 8 equations instead of 14). Pass `get_adjoints = true` for the full conjugate-closed state set. See issue [#319](https://github.com/qojulia/QuantumCumulants.jl/issues/319).

### Fixed

`show(::MIME"text/latex", ::MeanfieldEquations)` now wraps its output in `$$ … $$` with the `aligned` environment. The bare `\begin{align}` was mangled by Markdown renderers like Documenter (subscripts parsed as emphasis), so equation systems now render correctly in the docs and notebooks. The example pages are cleaned up accordingly. See issue [#310](https://github.com/qojulia/QuantumCumulants.jl/issues/310).

### Documentation

Added a documentation page, [Convergence and validity of the cumulant expansion](https://qojulia.github.io/QuantumCumulants.jl/stable/convergence/), collecting what is known about when a cumulant truncation can be trusted: the exactness of the expansion, the Marcinkiewicz obstruction and the special role of second order, the absence of a general convergence theorem, the good/bad/ugly benchmark behaviour, the central-spin non-monotonicity lessons, the connection to classical moment closure, the regimes in which the method is provably controlled, and practical guidance for users.


## [0.6.0]

### Fixed

Indexed sums with an index-dependent coefficient no longer lose their scope during the cumulant expansion. When a moment inside a sum `Σ_{i,j} c(i,j)·⟨A_i A_j⟩` exceeded the truncation order, the factorised product `c(i,j)·⟨A_i⟩⟨A_j⟩` was re-scoped by stamping the summation indices onto the first averaged leaf that used them, which left index-dependent coefficients (e.g. a `DoubleIndexedVariable` coupling `J(i,j)`) and any sibling leaves outside the sum with a dangling free index. `evaluate` then could not unroll them and `ODEProblem` construction failed with an `UndefVarError` on the dangling index. The factorised term is now wrapped back in a dedicated moment-layer indexed-sum node covering exactly the coefficient and leaves that share each bound index, so `evaluate` unrolls the whole product in lockstep and the concrete system builds and solves. See issues [#288](https://github.com/qojulia/QuantumCumulants.jl/issues/288) and [#198](https://github.com/qojulia/QuantumCumulants.jl/issues/198).

`scale` now handles a factorised indexed-sum node whose body is a product of several averages (e.g. `Σ_{i≠j} J(i,j)·⟨σ⁻ᵢ⟩⟨σ⁻ⱼ⟩`, from a two-atom coupling on an index-free observable). Previously it reconstructed the operator of the whole node, which fused the moment product `⟨σ⁻ᵢ⟩⟨σ⁻ⱼ⟩` into the operator product `σ⁻ᵢσ⁻ⱼ` and re-averaged it to a spurious higher-order moment `⟨σ⁻ᵢσ⁻ⱼ⟩`, leaving the scaled system above its truncation order and unclosed. `scale` now folds each averaged leaf of the sum body to the symmetry representative and charges the falling-factorial `N(N−1)…` prefactor from the node's scope, agreeing with `evaluate` to numerical precision.

`CorrelationFunction` no longer emits redundant equations for the conjugate of the ancilla operator. For a transition correlation such as the Mollow case ⟨σᵉᵍ(τ) σᵍᵉ(0)⟩, the closure previously also seeded the disconnected conjugate branch (states carrying the ancilla operator's adjoint σᵉᵍ_0) as extra unknowns, giving five equations where three suffice. Because the ancilla operator is inert, that branch never couples to the requested correlation, so it is now dropped: the closure keeps only the states `g(τ)` depends on, matching the documented Mollow example (three equations). The computed correlation and power spectrum are unchanged. See issue [#311](https://github.com/qojulia/QuantumCumulants.jl/issues/311).

Indexed sums whose summation index appears **only** in the coefficient (for example a per-site drive `Ω_l = Σ_k u(l,k)` built from a `DoubleIndexedVariable`) now give the correct populations. When such a sum was multiplied by an operator, the moment-layer re-wrapped the already-split off-diagonal body in `Σ` without carrying its `k ≠ l` constraint, so the diagonal was peeled off a second time and the `k = l` term was double-counted (`_carry_non_equal` now keeps every non-equal pair that references an operator index on the block). On an order-2-exact independent-atom model the excited-state population now matches the exact master equation to numerical precision (previously it read 0.447 against an exact 0.325). This fix also depends on the companion diagonal-split fix in SecondQuantizedAlgebra v0.9.3, which keeps the off-diagonal body free of the diagonal coefficient (`u(l,k)` rather than `u(l,k) + u(l,l)`). See issue [#198](https://github.com/qojulia/QuantumCumulants.jl/issues/198).

### Added

`CorrelationFunction` now accepts a vector of operators for `op1`, seeding one correlation τ-state ⟨opᵢ(τ) op2(0)⟩ per entry. This guarantees the listed operators are all included in the correlation system and fixes their order, which simplifies defining several correlations that share the same `op2`. The first entry is kept as the representative `op1`. See issue [#311](https://github.com/qojulia/QuantumCumulants.jl/issues/311) and PR [#309](https://github.com/qojulia/QuantumCumulants.jl/pull/309).

### Changed

The time-shifted operator in a `CorrelationFunction` now gets a fresh name. `op2`, re-embedded on the ancilla subspace, is renamed `<name>_0` (e.g. `a` becomes `a_0`), and the correlation prints as ⟨op1 op2_0⟩ rather than `op1*op2` collapsed onto a single subspace. The two-time structure ⟨op1(τ) op2(0)⟩ is now visible in the displayed equations, in `show`, and in the LaTeX output. See issue [#311](https://github.com/qojulia/QuantumCumulants.jl/issues/311).

## [0.5.6]

### Fixed

The effective-model comparison in the unique-squeezing example now matches the full model again. After the v0.5 rewrite the effective coupling `gΩ = g²/4Ω` was evaluated once at the global `N_global = 100`, but the effective Hamiltonian still multiplied it by a different `N_ = 69`. The original v0.4 version kept `g` symbolic in `N` (so `g² ∝ 1/N` and the intensive product `N g²` was independent of `N`), which made the mismatched `N` harmless; once `g` became a fixed number the effective model ran at 69% of the intended coupling and its squeezing curves fell short of the full model. The example now binds `gΩ` to `N` (see the new `bindings` support below), so the effective model is driven by the single knob `N` and the intensive `N g²` is restored by construction, bringing the two models back into agreement. See issue [#299](https://github.com/qojulia/QuantumCumulants.jl/issues/299).

`substitute` on a `MeanfieldEquations`/`NoiseMeanfieldEquations` works again. The v0.4 method `substitute(eqs, dict)` was dropped in the v0.5 rewrite, so the same call fell through to the generic `SymbolicUtils.substitute`, which treated the equation set as an opaque object and returned it unchanged. The call succeeded silently, having substituted nothing. This broke the minimal-closure idiom of injecting known-zero moments (e.g. by symmetry) to collapse a system onto exactly its observable subset: the substitution did nothing, `System` then correctly rejected the still-open set, and the natural next step of `complete!` closed the full hierarchy instead, ballooning the equation count. `substitute` now rewrites every RHS (and the noise drift RHSs of a `NoiseMeanfieldEquations`) through the graph, so the reduced closure is recovered.

### Added

`ModelingToolkitBase.System(eqs::MeanfieldEquations; name, kwargs...)` and the `NoiseMeanfieldEquations` method now forward extra keyword arguments to the underlying `ModelingToolkitBase.System` constructor. In particular `bindings = [p => expr, ...]` marks a parameter as *derived* from others, so it is computed from the supplied values rather than set independently. This restores the v0.4 convenience of expressing one parameter as a function of the rest (e.g. an intensive collective coupling in terms of the atom number), which had been lost when the rewrite stopped accepting symbolic expressions as parameter values.

`substitute!(eqs::AbstractMeanfieldEquations, dict)` is exported: the in-place counterpart of `substitute`, matching the `simplify!`/`modify_equations!` family. Both rewrite each RHS via the moment graph; `substitute` returns a fresh struct, `substitute!` mutates in place.

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

