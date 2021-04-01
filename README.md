# QuantumCumulants.jl
**QuantumCumulants.jl** is a package for the symbolic derivation of equations of motion for average values of quantum mechanical operators in Julia. The equations are derived using fundamental commutation relations of operators and can be automatically expanded in terms of cumulants to an arbitrary order. This results in a closed set of symbolic differential equations. Finally, these equations can be solved on a numeric level.

For the application of commutation relations **QuantumCumulants.jl** implements a simple noncommutative algebra, where any commutation relations are applied immediately. All other symbolic simplification and rewriting is done using the [Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl) package.

To obtain a numerical solution, equations derived with **QuantumCumulants.jl** can be converted to [ModelingToolkit.jl](https://github.com/SciML/ModelingToolkit.jl), which is part of the [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl) ecosystem and bridges the gap between symbolics and numerics.

## Development status

![CI](https://github.com/david-pl/QuantumCumulants.jl/workflows/CI/badge.svg) [![Codecov][codecov-img]][codecov-url] [![Documentation][docs-img]][docs-url]

Note that **QuantumCumulants.jl** is still at a very early stage of development.

## Installation

The package can be installed with

```julia
|pkg> add https://github.com/david-pl/QuantumCumulants.jl.git
```

## Documentation

Please refer to the latest [Documentation][docs-url] for more details and examples.

[codecov-url]: https://codecov.io/gh/david-pl/QuantumCumulants.jl/branch/master/
[codecov-img]: https://codecov.io/gh/david-pl/QuantumCumulants.jl/branch/master/graph/badge.svg

[docs-url]: https://david-pl.github.io/QuantumCumulants.jl/dev/
[docs-img]: https://img.shields.io/badge/docs-dev-blue.svg
