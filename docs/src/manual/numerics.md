```@meta
CurrentModule = QuantumCumulants
```

# Numerical solution

A closed [`MeanfieldEquations`](@ref) is a system of ordinary differential equations for complex averages. There are two ways to integrate it. The ModelingToolkit route builds a `System`, which handles time-dependent drives, registered functions and stochastic noise. The direct backend builds an `ODEProblem` without ModelingToolkit, which is faster to set up and cheap to re-parametrise, but accepts only polynomial deterministic drifts.

```@contents
Pages = ["numerics.md"]
Depth = 2:3
```

We integrate the phase-invariant laser of [Closing the hierarchy](@ref closure).

```@example numerics
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
ps = Dict(Δ => 0.0, g => 1.5, γ => 0.25, κ => 1.0, ν => 4.0)
nothing # hide
```

## ModelingToolkit

Each moment becomes one unknown ``u(t)`` of the `System`, and every symbol on a right-hand side that is not a moment becomes a parameter. [`moment_variable_map`](@ref) shows the assignment.

```@example numerics
using ModelingToolkitBase, OrdinaryDiffEqTsit5

sys = mtkcompile(System(eqs; name = :laser))
moment_variable_map(eqs)
```

```@docs
ModelingToolkitBase.System(::MeanfieldEquations)
moment_variable_map
```

### Initial values

ModelingToolkit takes initial values as a map from unknowns to numbers. [`initial_values`](@ref) builds it from a vector aligned with [`states`](@ref), or from a numerical quantum state through SecondQuantizedAlgebra's [numeric conversion](https://qojulia.github.io/SecondQuantizedAlgebra.jl/stable/manual/averages_numerics/). Here the cavity starts in a coherent state and the atom in its ground state.

```@example numerics
using QuantumOpticsBase
ψ0 = coherentstate(FockBasis(10), 0.5) ⊗ nlevelstate(NLevelBasis(2), 1)
u0 = initial_values(eqs, ψ0)
prob = ODEProblem(sys, merge(u0, ps), (0.0, 10.0))
sol = solve(prob, Tsit5())
nothing # hide
```

```@docs
initial_values
parameter_map
```

### Reading the solution

[`get_solution`](@ref) resolves an operator or an average to the matching state, including states stored only through their complex conjugate, and returns a callable trajectory.

```@example numerics
n = get_solution(sol, a' * a, eqs)
n(10.0)
```

```@docs
get_solution
```

## Direct backend

The direct backend lowers the equations once to a polynomial representation and evaluates them with a compact kernel. It is selected explicitly with `backend = KernelBackend()` and takes the same initial values, as a map, a vector aligned with [`states`](@ref), or the quantum state itself.

```@example numerics
prob_k = ODEProblem(eqs, u0, (0.0, 10.0), ps; backend = KernelBackend())
sol_k = solve(prob_k, Tsit5())
get_solution(sol_k, a' * a, eqs)(10.0)
```

Parameters can be changed in place without lowering the equations again, which suits parameter sweeps. Updates may be partial. For independent sweeps on several threads, give each task its own copy of the problem.

```@example numerics
update_parameters!(prob_k, Dict(g => 2.0))
get_solution(solve(prob_k, Tsit5()), a' * a, eqs)(10.0)
```

The backend accepts polynomial drifts with state-independent coefficients. Time-dependent coefficients and other functions of the states have to go through ModelingToolkit. For implicit solvers the backend also provides the exact sparse Jacobian with `jac = true`. This requires an unfolded closure (`get_adjoints = true` in [`complete`](@ref)), since a single complex Jacobian is not the full derivative of ``f(u, u^*)``.

```@docs
KernelBackend
SciMLBase.ODEProblem(::MeanfieldEquations, ::Any, ::Any, ::Any)
update_parameters!
```
