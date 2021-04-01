# Introduction

**QuantumCumulants.jl** offers a practical approach to the application of the generalized cumulant expansion method in Quantum Optics: operators are often represented by matrices on a Hilbert space, where a suitable basis has been chosen. These matrices can quickly become so large that they can no longer be stored. On a more abstract level, however, operators form a noncommutative alebgra that follows fundamental commutation relations. This is where **QuantumCumulants.jl** comes in. The basic working principle boils down to the following steps:

* The model (Hamiltonian) is specified.
* Equations of motion for average values are derived. This is done by using the fundamental commutation relations of operators. The resulting equations are stored as symbolic equations using the [**Symbolics.jl**](https://github.com/JuliaSymbolics/Symbolics.jl) framework, which is also used for any additional simplification and rewriting.
* Then follows the key step: the equations of motion the averages are truncated at a specified order neglecting higher-order quantum correlations using the generalized cumulant expansion method. This results in a closed set of *c*-number ordinary differential equations.
* Finally, the symbolic system of equations can be turned into an `ODESystem` of the [**ModelingToolkit.jl**](https://github.com/SciML/ModelingToolkit.jl) framework which bridges the gap from symbolics to numerics. This makes it straightforward to obtain a solution of the time dynamics of a system within the  [**DifferentialEquations.jl**](https://diffeq.sciml.ai/latest/) ecosystem.


### Installation

**QuantumCumulants.jl** is not officially released (yet) and so has to be installed directly from GitHub. This can be done with

```julia
|pkg> add https://github.com/david-pl/QuantumCumulants.jl.git
```

For a full list of functions, check out the [API](@ref).
