# Introduction

**Qumulants.jl** offers a practical approach to the application of the generalized cumulant expansion method in Quantum Optics: operators are often represented by matrices on a Hilbert space, where a suitable basis has been chosen. These matrices can quickly become so large that they can no longer be stored. On a more abstract level, however, operators form a noncommutative alebgra that follows fundamental commutation relations. This is where **Qumulants.jl** comes in. The basic working principle boils down to the following steps:

* The model (Hamiltonian) is specified.
* Equations of motion (Heisenberg equations) are symbolically derived for the system operators by using their fundamental commutation relations.
* Then follows the key step: the equations of motion for the noncommutative operators are averaged and truncated at a specified order neglecting higher-order quantum correlations using the generalized cumulant expansion method. This results in a closed set of *c*-number ordinary differential equations.
* Finally, the symbolic system of equations can be turned into a native Julia function which is directly usable in the [**DifferentialEquations.jl**](https://diffeq.sciml.ai/latest/) framework. This makes it straightforward to obtain a solution of the time dynamics of a system.


### Installation

**Qumulants.jl** is not officially released (yet) and so has to be installed directly from GitHub. This can be done with

```julia
|pkg> add https://github.com/david-pl/Qumulants.jl.git
```

For a full list of functions, check out the [API](@ref).
