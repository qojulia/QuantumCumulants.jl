```@meta
CurrentModule = QuantumCumulants
```

# Getting started

QuantumCumulants turns a Hamiltonian and a set of collapse operators into a closed system of ordinary differential equations for expectation values. The workflow has four steps: derive the Heisenberg equations for a few averages, truncate them with the cumulant expansion, close the hierarchy, and integrate the result. This page runs through all four for a single-atom laser, a cavity mode ``a`` coupled to a two-level atom that is incoherently pumped,

```math
H = \Delta a^\dagger a + g \left(a^\dagger \sigma^{ge} + \sigma^{eg} a\right),
```

with cavity loss ``\kappa``, atomic decay ``\gamma`` and pumping ``\nu``. The later pages of the manual take up each step in detail.

```@contents
Pages = ["getting_started.md"]
Depth = 2:3
```

## The model

The Hilbert spaces, operators and symbolic parameters come from [SecondQuantizedAlgebra.jl](https://qojulia.github.io/SecondQuantizedAlgebra.jl/stable/), which QuantumCumulants re-exports. Its manual on [Hilbert spaces and operators](https://qojulia.github.io/SecondQuantizedAlgebra.jl/stable/manual/spaces_and_operators/) covers everything used here.

```@example laser
using Latexify # hide
set_default(double_linebreak = true) # hide
using QuantumCumulants

h = FockSpace(:cavity) ⊗ NLevelSpace(:atom, (:g, :e))
@qnumbers a::Destroy(h)
σ(i, j) = Transition(h, :σ, i, j)
@variables Δ::Real g::Real γ::Real κ::Real ν::Real

H = Δ * a' * a + g * (a' * σ(:g, :e) + a * σ(:e, :g))
J = [a, σ(:g, :e), σ(:e, :g)]
rates = [κ, γ, ν]
nothing # hide
```

## Equations of motion

The photon number, the excited-state population and the atom-field coherence form a natural starting set. [`meanfield`](@ref) derives their equations of motion, and `order = 2` applies the cumulant expansion to second order on the fly.

```@example laser
ops = [a' * a, σ(:e, :e), a' * σ(:g, :e)]
eqs = meanfield(ops, H, J; rates, order = 2)
```

```@docs
meanfield
MeanfieldEquations
AbstractMeanfieldEquations
```

## Closing the system

The right-hand sides above reference averages that have no equation of their own, such as ``\langle a \rangle``. [`complete`](@ref) derives the missing equations until the set is self-contained.

```@example laser
eqs_c = complete(eqs)
nothing # hide
```

The first-order averages in this set stay zero from a zero initial state, because the model is invariant under a global phase rotation. [Closing the hierarchy](@ref closure) shows how to exclude them from the start.

## Solving

The closed system becomes a ModelingToolkit `System`, which is compiled and integrated with OrdinaryDiffEq. [`get_solution`](@ref) then turns the solution into a callable trajectory for any average.

```@example laser
using ModelingToolkitBase, OrdinaryDiffEqTsit5

sys = mtkcompile(System(eqs_c; name = :laser))
ps = Dict(Δ => 0.0, g => 1.5, γ => 0.25, κ => 1.0, ν => 4.0)
prob = ODEProblem(sys, merge(initial_values(eqs_c), ps), (0.0, 10.0))
sol = solve(prob, Tsit5())

using Plots
ts = range(0.0, 10.0; length = 200)
plot(ts, real.(get_solution(sol, a' * a, eqs_c).(ts)); label = "photon number", xlabel = "t")
plot!(ts, real.(get_solution(sol, σ(:e, :e), eqs_c).(ts)); label = "excited-state population")
```

[Numerical solution](@ref) covers initial states, parameter maps, and a direct backend that skips ModelingToolkit altogether.
