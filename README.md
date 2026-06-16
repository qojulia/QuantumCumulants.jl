# QuantumCumulants.jl

[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://qojulia.github.io/QuantumCumulants.jl/dev/)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://qojulia.github.io/QuantumCumulants.jl/stable/)
[![Codecov](https://codecov.io/gh/qojulia/QuantumCumulants.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/qojulia/QuantumCumulants.jl/branch/master/)
[![Benchmarks](https://github.com/qojulia/QuantumCumulants.jl/actions/workflows/Benchmarks.yaml/badge.svg?branch=main)](https://qojulia.github.io/QuantumCumulants.jl/benchmarks/)

[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![jet](https://img.shields.io/badge/%F0%9F%9B%A9%EF%B8%8F_tested_with-JET.jl-233f9a)](https://github.com/aviatesk/JET.jl)

**QuantumCumulants.jl** is a package for the symbolic derivation of mean-field equations for quantum mechanical operators in Julia. The equations are derived using fundamental commutation relations of operators. When averaging these equations they can be automatically expanded in terms of cumulants to an arbitrary order (generalized mean-field approximation). This results in a closed set of symbolic differential equations, which can also be solved numerically.

For the application of commutation relations **QuantumCumulants.jl** implements a simple noncommutative algebra, where any commutation relations are applied immediately. All other symbolic simplification and rewriting is done using the [Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl) package.

The basic working principle boils down to the following steps:

* The model (Hamiltonian) is specified.

* Equations of motion for average values are derived from the fundamental commutation relations of operators. The resulting equations are stored as symbolic equations using the [Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl) framework, which also handles any additional simplification and rewriting.

* The key step: the equations of motion for the averages are truncated at a specified order, neglecting higher-order quantum correlations via the generalized cumulant expansion method. This yields a closed set of *c*-number ordinary differential equations.

* Finally, the symbolic system of equations is turned into a `System` from [ModelingToolkitBase.jl](https://github.com/SciML/ModelingToolkitBase.jl), which bridges the gap from symbolics to numerics. This makes it straightforward to obtain the time dynamics of a system within the [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl) ecosystem.

If you only need the second quantized algebra, you can depend on [SecondQuantizedAlgebra.jl](https://github.com/qojulia/SecondQuantizedAlgebra.jl) directly.


## Relationship to SecondQuantizedAlgebra.jl

QuantumCumulants is the *cumulant* layer. The operator algebra it builds on (Hilbert spaces, operators, indices, symbolic sums, normal ordering, numeric conversion) lives in [SecondQuantizedAlgebra.jl](https://github.com/qojulia/SecondQuantizedAlgebra.jl) (SQA), which QuantumCumulants re-exports, so `using QuantumCumulants` gives you the full algebra surface. When you need the details of *building a model* (defining `FockSpace`/`NLevelSpace`, `Destroy`/`Transition`, `Index`/`Σ`), reach for SQA's [documentation](https://qojulia.github.io/SecondQuantizedAlgebra.jl/stable/); the documentation here focuses on what QuantumCumulants adds on top: deriving mean-field equations, the cumulant expansion, completion, scaling, correlations, noise, and the bridge to a numerical solution.

## Installation

**QuantumCumulants.jl** is a registered Julia package and can be installed with the package manager:

```julia
pkg> add QuantumCumulants
```

For a full list of functions, see the [API documentation](https://qojulia.github.io/QuantumCumulants.jl/stable/api/).

## Short example

To briefly illustrate how **QuantumCumulants.jl** works, here's how you can implement a first-order mean-field model of a laser with a single atom as a gain medium:

```julia
using QuantumCumulants

h_cav = FockSpace(:cavity)
h_atom = NLevelSpace(:atom, (:g,:e))
h = tensor(h_cav, h_atom)

@variables Δ::Real g::Real κ::Real γ::Real ν::Real
@qnumbers a::Destroy(h)
σ(i, j) = Transition(h, :σ, i, j)

H = Δ*a'*a + g*(a'*σ(:g,:e) + a*σ(:e,:g))
J = [a,σ(:g,:e),σ(:e,:g)]
rates = [κ,γ,ν]

eqs = meanfield([a,σ(:g,:e),σ(:e,:e)], H, J; rates=rates, order=1)

using ModelingToolkitBase, OrdinaryDiffEq
sys = mtkcompile(System(eqs; name=:laser))
p0 = [Δ=>0, g=>1.5, κ=>1, γ=>0.25, ν=>4]
u0 = ComplexF64[1e-2, 0, 0]
prob = ODEProblem(sys, merge(initial_values(eqs, u0), Dict(p0)), (0.0,50.0))
sol = solve(prob,RK4())

using Plots
ts = range(0.0, 50.0; length=200)
n = abs2.(get_solution(sol, a, eqs).(ts))
plot(ts, n, xlabel="t", label="n")
```

![photon-number](https://user-images.githubusercontent.com/18166442/114183684-3ae76080-9944-11eb-9d21-94bf4069bb60.png)


The above code implements the Jaynes-Cummings Hamiltonian describing an optical cavity mode that couples to a two-level atom. Additionally, the decay processes are specified. Then, mean-field equations for the average values of the operators `[a,σ(:g,:e),σ(:e,:e)]` are derived and expanded to first order (average values of products are factorized). For the numerical solution a `System` (from [ModelingToolkitBase.jl](https://github.com/SciML/ModelingToolkitBase.jl)) is created and solved with the [OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl) library. Finally, the time dynamics of the photon number `n` is plotted.


## Citing

If you find **QuantumCumulants.jl** useful in your research, please consider citing [this paper](https://arxiv.org/abs/2105.01657):

```bib
@article{plankensteiner2022quantumcumulants,
  doi = {10.22331/q-2022-01-04-617},
  url = {https://doi.org/10.22331/q-2022-01-04-617},
  title = {Quantum{C}umulants.jl: {A} {J}ulia framework for generalized mean-field equations in open quantum systems},
  author = {Plankensteiner, David and Hotter, Christoph and Ritsch, Helmut},
  journal = {{Quantum}},
  issn = {2521-327X},
  publisher = {{Verein zur F{\"{o}}rderung des Open Access Publizierens in den Quantenwissenschaften}},
  volume = {6},
  pages = {617},
  month = jan,
  year = {2022}
}
```
