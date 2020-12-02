# Qumulants.jl
**Qumulants.jl** is a package for the symbolic derivation of Heisenberg equations in Julia. Averages over the resulting equations can be automatically expanded in terms of cumulants to an arbitrary order. This procedure yields a system of symbolic *c*-number differential equations. Finally, these *c*-number equations can be mapped to a function which can be solved using [DifferentialEquations.jl](http://docs.juliadiffeq.org/latest/).

For the application of commutation relations and general simplification, **Qumulants.jl** uses [SymbolicUtils.jl](https://github.com/JuliaSymbolics/SymbolicUtils.jl).


## Development status

![CI](https://github.com/david-pl/Qumulants.jl/workflows/CI/badge.svg) [![Codecov][codecov-img]][codecov-url] [![Documentation][docs-img]][docs-url]

**Qumulants.jl** is still at a very early stage of development. **Expect bugs!**


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
