# Qumulants.jl
**Qumulants.jl** is a package for the symbolic derivation of Heisenberg equations equation using in Julia. Averages over the resulting equations can be automatically expanded in terms of cumulants to an arbitrary order. This procedure yields a system of symbolic *c*-number differential equations stored as [SymPy](https://github.com/JuliaPy/SymPy.jl) expression. Finally, these *c*-number equations can be mapped to a function which can be solved using [DifferentialEquations.jl](http://docs.juliadiffeq.org/latest/).


## Development status

**Qumulants.jl** is still at a very early stage of development. **Expect bugs!**


## Basic usage

The package can be installed with

```julia
|pkg> add git@github.com:david-pl/Qumulants.jl.git
```

The basic usage is probably best illustrated with a brief example. In the following, we solve a simple model for a single-atom laser.

**Note:** It is recommended to use an editor that supports LaTeX output such as [Jupyter](https://jupyter.org/). This enables LaTeX printing of operators and equations.

We start by loading the package, defining some parameters and the photonic annihilation operator `a` as well as the atomic lowering operator `s`. This allows us to quickly write down the Hamiltonian and the collapse operators of the system with their corresponding decay rates.

```julia
using Qumulants

# Define parameters
Δ, g, γ, κ, ν = (0, 1.5, 0.25, 1, 4)

# Define the fundamental operators
a = Destroy(:a) ⊗ Identity()
s = Identity() ⊗ Transition(:σ,:g,:e,(:g,:e))

# Hamiltonian
H = Δ*a'*a + g*(a'*s + a*s')

# Collapse operators
J = [a,s,s']
rates = [κ,γ,ν]
```

Now, we define a list of operators of which we want to compute the Heisenberg equations. We will only consider products of two operators. This is because later we will compute the dynamics of the system up to second order. The first-order contributions are always zero and can therefore be neglected (you can try adding `a` and `s` to the list of operators `ops` in order to see that yourself).

```julia
# Derive a set of Heisenberg equations
ops = [a'*a,s'*s,a'*s]
he = heisenberg(ops,H,J;rates=rates)
```

The equations derived above are differential equations for operators. In order to convert them to *c*-number equations, we need to average over them. To obtain a closed set of equations, we expand higher-order products to second order.

```julia
# Average the above equations and expand to second order
he_avg = average(he,2)
```

Finally, we can generate Julia code from the above set of equations which can be solved directly using the [OrdinaryDiffEq](https://github.com/JuliaDiffEq/OrdinaryDiffEq.jl). The kwarg `set_unknowns_zero=true` removes unknown averages from the equations (in our case first-order averages which are zero anyway).

```julia
# Generate a Julia function that to solve numerically
f = generate_ode(he_avg;set_unknowns_zero=true)

# Solve the system using the OrdinaryDiffEq package
import OrdinaryDiffEq; ODE=OrdinaryDiffEq
u0 = zeros(ComplexF64,length(ops))
prob = ODE.ODEProblem(f,u0,(0.0,10.0))
sol = ODE.solve(prob,ODE.Tsit5())
```

The photon number of our laser and the excited state population of the atom are now stored in the first two fields of `sol.u`.
