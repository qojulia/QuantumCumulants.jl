# The mean-field pipeline

This page walks through the steps QuantumCumulants takes you through, from a Hamiltonian to a numerical solution: deriving the Heisenberg equations, truncating them with the cumulant expansion, closing the system, and handing it to a numerical solver.

The building blocks of a model (Hilbert spaces, operators, symbolic parameters, commutation relations) belong to [SecondQuantizedAlgebra.jl](https://github.com/qojulia/SecondQuantizedAlgebra.jl) (SQA), which QuantumCumulants re-exports. If you have not met them yet, read SQA's [implementation guide](https://qojulia.github.io/SecondQuantizedAlgebra.jl/stable/implementation/) first; this page assumes you can already build a Hamiltonian. As a one-line reminder, here is the driven cavity used throughout:

```@example meanfield
using QuantumCumulants
h = FockSpace(:fock)
@qnumbers a::Destroy(h)
@variables ω::Real η::Real
H = ω*a'*a + η*(a + a') # driven cavity Hamiltonian
nothing # hide
```

[`FockSpace`](@ref), [`Destroy`](@ref), [`@qnumbers`](@ref) and the symbolic scalar [`@variables`](https://docs.sciml.ai/Symbolics/stable/) all come from SQA / Symbolics. For many-body systems built from [`Index`](@ref) and [`Σ`](@ref), see [Indexed and scaled systems](@ref).

## Cumulant expansion

Equations of motion are written for expectation values. SQA's [`average`](@ref) converts an operator product to a c-number, but a closed set of equations requires truncating at a specified order. The order of an average is the number of operator factors in the product:

```@example cumulant
using QuantumCumulants # hide
h = FockSpace(:fock)
@qnumbers a::Destroy(h)

get_order(average(a))      # 1
get_order(average(a'*a))   # 2
nothing # hide
```

The [`cumulant_expansion`](@ref) expresses any average through averages up to a chosen order (see the [theory section](@ref theory)):

```@example cumulant
cumulant_expansion(average(a'*a), 1)
nothing # hide
```

## Deriving equations of motion

[`meanfield`](@ref) derives the Heisenberg (quantum Langevin) equations for a list of operators under a Hamiltonian `H` and, optionally, a list of collapse operators `J` with their `rates`. Each jump adds the Lindblad term

```math
\sum_i r_i \left( J_i^\dagger \mathcal{O} J_i - \tfrac{1}{2}\{J_i^\dagger J_i, \mathcal{O}\} \right).
```

Pass `order` to apply the [`cumulant_expansion`](@ref) immediately:

```@example meanfield
using Latexify # hide
set_default(double_linebreak=true) # hide
me = meanfield([a, a'*a], H; order=2)
```

The result is a [`MeanfieldEquations`](@ref) holding both the operator-level and the averaged equations. Returned right-hand sides are left unsimplified; call [`simplify!`](@ref) (or `Symbolics.simplify`) on what you want to inspect.

A derived system is generally not closed: right-hand sides reference averages that are not yet on a left-hand side. [`find_missing`](@ref) lists them, and [`complete`](@ref) / [`complete!`](@ref) close the system by deriving an equation for each until the set is self-contained:

```@example meanfield
me_complete = complete(me)
nothing # hide
```

For opt-in measurement-backaction noise (`efficiencies=...`) and retrodiction (`direction=Backward()`), see [Noise & measurement backaction](@ref). For permutation-symmetric many-body systems, see [Indexed and scaled systems](@ref).

## Numerical solution

A [`MeanfieldEquations`](@ref) is converted to a `ModelingToolkitBase.System` of ordinary differential equations, compiled with `mtkcompile`, and solved with [OrdinaryDiffEq](https://github.com/SciML/OrdinaryDiffEq.jl):

```@example meanfield
using ModelingToolkitBase
sys = mtkcompile(System(me_complete; name=:cavity))

using OrdinaryDiffEqTsit5
u0 = zeros(ComplexF64, length(me_complete.states))
p0 = [ω => 1.0, η => 0.1]
prob = ODEProblem(sys, merge(initial_values(me_complete, u0), Dict(p0)), (0.0, 1.0))
sol = solve(prob, Tsit5())
nothing # hide
```

[`get_solution`](@ref) substitutes a symbolic average into the compiled solution and returns a callable in `t`. It works both for averages on the left-hand side of the system and for derived products that are not:

```@example meanfield
ts  = range(0.0, 1.0; length=100)
a_t = get_solution(sol, a, me_complete).(ts)
n_t = real.(get_solution(sol, a'*a, me_complete).(ts))
nothing # hide
```

### Computing the initial state

A trivial `u0 = zeros(...)` is rarely the physical initial condition. [`initial_values`](@ref) computes the expectation values of every state from a numerical `Ket` or density operator, using SQA's [`numeric_average`](@ref) / [`to_numeric`](@ref) bridge to [QuantumOpticsBase.jl](https://github.com/qojulia/QuantumOpticsBase.jl). For example, a coherent initial state of the cavity:

```@example meanfield
using QuantumOpticsBase
b = FockBasis(10)
psi_0 = coherentstate(b, 0.3 + 0.4im)
u0 = initial_values(me_complete, psi_0)
nothing # hide
```

Mixed states work the same way by passing a density operator (`initial_values(me, dm(psi_0))`). The details of the symbolic-to-numeric conversion (`to_numeric`, `numeric_average`, `NLevelSpace` level maps, lazy tensor products for large systems) are documented in SQA's [numeric conversion section](https://qojulia.github.io/SecondQuantizedAlgebra.jl/stable/implementation/#Numeric-conversion).

## Extending the operator algebra

Adding a new operator type (a custom `QSym` with its commutation hooks) is an SQA concern, not a QuantumCumulants one. SQA's [implementation guide](https://qojulia.github.io/SecondQuantizedAlgebra.jl/stable/implementation/) documents the `QSym` interface and walks through a worked example; once a type is defined there, `meanfield` and the rest of this pipeline handle it automatically.
