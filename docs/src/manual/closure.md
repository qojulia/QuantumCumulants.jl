```@meta
CurrentModule = QuantumCumulants
```

# [Closing the hierarchy](@id closure)

Truncation bounds the order of the averages, but it does not bound their number. A truncated equation still refers to averages of lower or equal order that have no equation of their own, and the system is closed only once every average on a right-hand side is also a state. Closing is a reachability sweep: derive an equation for each missing average, collect the averages its right-hand side introduces, and repeat until nothing new appears.

```@contents
Pages = ["closure.md"]
Depth = 2:3
```

We use the single-atom laser of [Getting started](@ref) throughout.

```@example closure
using QuantumCumulants

h = FockSpace(:cavity) ⊗ NLevelSpace(:atom, (:g, :e))
@qnumbers a::Destroy(h)
σ(i, j) = Transition(h, :σ, i, j)
@variables Δ::Real g::Real γ::Real κ::Real ν::Real

H = Δ * a' * a + g * (a' * σ(:g, :e) + a * σ(:e, :g))
J = [a, σ(:g, :e), σ(:e, :g)]
eqs = meanfield([a' * a], H, J; rates = [κ, γ, ν], order = 2)
nothing # hide
```

## Finding missing averages

[`find_missing`](@ref) lists the averages that still lack an equation, and [`closure_report`](@ref) summarises the state of the whole set.

```@example closure
find_missing(eqs)
```

```@example closure
closure_report(eqs)
```

```@docs
find_missing
closure_report
```

## Completing a system

[`complete`](@ref) runs the sweep to its end. A moment and its complex conjugate count as one state, since the equation of one is the conjugate of the other.

```@example closure
eqs_c = complete(eqs)
closure_report(eqs_c)
```

```@docs
complete
complete!
```

## Filtering averages

A symmetry often guarantees that some averages vanish at all times. The laser is invariant under a global phase rotation ``a \to a e^{i\phi}``, ``\sigma^{ge} \to \sigma^{ge} e^{i\phi}``, so every average whose operators carry a net phase stays zero from a phase-invariant initial state. Passing a predicate as `filter_func` drops these averages during the sweep and keeps the system minimal.

```@example closure
import QuantumCumulants.SecondQuantizedAlgebra as SQA

phase(op) = SQA.is_create(op) - SQA.is_destroy(op) +
    (SQA.is_transition(op) ? sign(op.l1 - op.l2) : 0)
phase_invariant(avg) = iszero(sum(phase, SQA.get_operators(SQA.undo_average(avg)); init = 0))

eqs_f = complete(eqs; filter_func = phase_invariant)
```

The filter must encode a property of the dynamics, not a preference. Dropping an average that does evolve silently changes the physics.
