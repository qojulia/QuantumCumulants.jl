# Changelog

All notable changes to QuantumCumulants.jl will be documented in this file.

## [0.7.1]

Move to SecondQuantizedAlgebra v0.11

### Changed

Closure derivation now skips the post-cumulant ground-projector reduction when the operator expression cannot contain ordinary N-level transitions. Required ground reductions are serialized around the shared SymbolicUtils tree traversal, avoiding lazy argument-cache races during threaded closure construction while leaving the derived equations unchanged.

## [0.7.0]

### Changed

Moved to SecondQuantizedAlgebra v0.10. Load a numeric backend yourself (`using QuantumOpticsBase` or `using QuantumToolbox`), and Hermitian-operator averages now carry the `Real` symtype.
QuantumCumulants now relies only on the public SecondQuantizedAlgebra API. See issue [#300](https://github.com/qojulia/QuantumCumulants.jl/issues/300).

**Breaking.** `complete!` and `complete` now default to `get_adjoints = false`, keeping one representative per conjugate pair instead of a separate state for every conjugate moment. The partner is recovered by `conj` at code generation, so the dynamics are unchanged while the closed system is smaller (the optomechanical-cooling example closes at 8 equations instead of 14). Pass `get_adjoints = true` for the full conjugate-closed state set. See issue [#319](https://github.com/qojulia/QuantumCumulants.jl/issues/319).

### Fixed

`show(::MIME"text/latex", ::MeanfieldEquations)` now wraps its output in `$$ … $$` with the `aligned` environment. The bare `\begin{align}` was mangled by Markdown renderers like Documenter (subscripts parsed as emphasis), so equation systems now render correctly in the docs and notebooks. The example pages are cleaned up accordingly. See issue [#310](https://github.com/qojulia/QuantumCumulants.jl/issues/310).

### Documentation

Added a documentation page, [Convergence and validity of the cumulant expansion](https://qojulia.github.io/QuantumCumulants.jl/stable/convergence/), collecting what is known about when a cumulant truncation can be trusted: the exactness of the expansion, the Marcinkiewicz obstruction and the special role of second order, the absence of a general convergence theorem, the good/bad/ugly benchmark behaviour, the central-spin non-monotonicity lessons, the connection to classical moment closure, the regimes in which the method is provably controlled, and practical guidance for users.

Added the official QuantumCumulants logo to the repository landing page and documentation, together with a favicon for bookmarks and social-preview metadata for shared documentation links.


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
