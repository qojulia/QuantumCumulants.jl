# Qumulants.jl
**Qumulants.jl** is a package for the symbolic derivation of Heisenberg equations in Julia. Averages over the resulting equations can be automatically expanded in terms of cumulants to an arbitrary order. This procedure yields a system of symbolic *c*-number differential equations. Finally, these *c*-number equations can be solved on a numeric level.

The noncommutative symbolic algebra, the application of commutation relations as well as general simplification, are implemented using [SymbolicUtils.jl](https://github.com/JuliaSymbolics/SymbolicUtils.jl).

To obtain a numerical solution, equations derived with **Qumulants.jl** can be converted to [ModelingToolkit.jl](https://github.com/SciML/ModelingToolkit.jl), which is part of the [DifferentialEquations.jl] ecosystem.

## Development status

![CI](https://github.com/david-pl/Qumulants.jl/workflows/CI/badge.svg) [![Codecov][codecov-img]][codecov-url] [![Documentation][docs-img]][docs-url]

Note that **Qumulants.jl** is still at a very early stage of development.

## Installation

The package can be installed with

```julia
|pkg> add https://github.com/david-pl/Qumulants.jl.git
```

## Documentation

Please refer to the latest [Documentation][docs-url] for more details and examples.

[codecov-url]: https://codecov.io/gh/david-pl/Qumulants.jl/branch/master/
[codecov-img]: https://codecov.io/gh/david-pl/Qumulants.jl/branch/master/graph/badge.svg

[docs-url]: https://david-pl.github.io/Qumulants.jl/dev/
[docs-img]: https://img.shields.io/badge/docs-dev-blue.svg
