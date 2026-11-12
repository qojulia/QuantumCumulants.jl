```@meta
CurrentModule = QuantumCumulants
```

# Working with equations

A [`MeanfieldEquations`](@ref) is more than a list of equations. It holds the cumulant hierarchy as a graph of moments, and every function on this page reads or rewrites that graph. Changes therefore survive later steps: a modified right-hand side is still modified after [`complete`](@ref), [`scale`](@ref) or [`evaluate`](@ref).

```@contents
Pages = ["equations.md"]
Depth = 2:3
```

The examples continue with the phase-invariant laser of [Closing the hierarchy](@ref closure).

```@example equations
using QuantumCumulants
import QuantumCumulants.SecondQuantizedAlgebra as SQA

h = FockSpace(:cavity) ⊗ NLevelSpace(:atom, (:g, :e))
@qnumbers a::Destroy(h)
σ(i, j) = Transition(h, :σ, i, j)
@variables Δ::Real g::Real γ::Real κ::Real ν::Real

H = Δ * a' * a + g * (a' * σ(:g, :e) + a * σ(:e, :g))
J = [a, σ(:g, :e), σ(:e, :g)]

phase(op) = SQA.is_create(op) - SQA.is_destroy(op) +
    (SQA.is_transition(op) ? sign(op.l1 - op.l2) : 0)
phase_invariant(avg) = iszero(sum(phase, SQA.get_operators(SQA.undo_average(avg)); init = 0))

eqs = complete(meanfield([a' * a], H, J; rates = [κ, γ, ν], order = 2); filter_func = phase_invariant)
nothing # hide
```

## Inspecting the moments

Each equation belongs to one tracked moment. The operator products, their averages, and the pairing between the two are available separately.

```@example equations
operators(eqs)
```

```@example equations
moments(eqs)
```

```@docs
states
operators
moments
```

## Simplification

Derivation leaves right-hand sides unsimplified, since simplification is expensive and rarely needed before numerics. Call it when you want to read the equations.

```@example equations
simplify(eqs)
```

```@docs
simplify!
```

## Substitution

Substituting a parameter or a vanishing average rewrites every right-hand side at once, e.g. to go on resonance.

```@example equations
substitute(eqs, Dict(Δ => 0))
```

```@docs
substitute!
```

## Adding terms by hand

Some terms do not come from a Hamiltonian or a Lindblad dissipator, such as a phenomenological drive or a coupling to an external reservoir. [`modify_equations`](@ref) rewrites each right-hand side with a function of the left-hand operator. Here we add a coherent drive of the cavity to the laser.

```@example equations
@variables η::Real
Hd = η * (a + a')
add_drive(lhs, rhs) = rhs + cumulant_expansion(average(1im * commutator(Hd, lhs)), 2)
eqs_d = modify_equations(eqs, add_drive)
nothing # hide
```

The drive breaks the phase invariance, so the modified system is no longer closed. [`find_missing`](@ref) reports the averages that the new terms introduced.

```@example equations
find_missing(eqs_d)
```

The rewrite is recorded on the system, and [`complete`](@ref) applies it to every equation it derives afterwards. The result is the same as putting the drive into the Hamiltonian from the start.

```@example equations
eqs_dc = complete(eqs_d)
closure_report(eqs_dc)
```

[`substitute!`](@ref) is recorded in the same way, so a parameter set to zero stays zero in later equations too.

```@docs
modify_equations
modify_equations!
```
