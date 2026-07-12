# CumulantHomotopy.jl

CumulantHomotopy.jl is a proposed Julia package for computing physically selected stationary solutions of truncated quantum cumulant and moment hierarchies.

The package is intended to work closely with QuantumCumulants.jl. QuantumCumulants.jl generates the hierarchy. CumulantHomotopy.jl follows candidate stationary branches across physical parameters and truncation order, then evaluates whether the resulting roots remain consistent with the original open quantum system.

## Motivation

A truncated cumulant hierarchy produces a nonlinear polynomial system

```math
F_n(x; p) = 0,
```

where `n` is the truncation order, `x` contains retained moments or cumulants, and `p` contains physical parameters.

Even when the exact Liouvillian has one stationary density operator, the truncated polynomial system may have many complex or real roots. Some roots are closure artifacts. Some violate adjoint relations or moment positivity. Some approximate metastable states. A generic nonlinear solver has no reliable way to select the physically relevant branch.

CumulantHomotopy.jl will address this selection problem by combining:

1. Continuation from a known physical reference state.
2. Continuation between hierarchy orders using the actual closure map.
3. Automatic realification and symmetry reduction.
4. Physical diagnostics based on stationarity, adjoint consistency, realizability, and cross order convergence.
5. Optional semidefinite validation of candidate moment sequences.

## Central principle

The package will not claim that a finite truncation has one algebraic root. Instead, it will search for a distinguished sequence of roots

```math
x_1, x_2, \ldots, x_n, \ldots
```

whose fixed low order moments converge as the hierarchy order increases and whose physical consistency improves.

The primary output is therefore a stationary branch together with evidence and diagnostics, not an unqualified claim of uniqueness.

## Proposed workflow

```julia
using QuantumCumulants
using CumulantHomotopy

hierarchies = cumulant_hierarchies(problem; orders = 2:8)

result = stationary_sequence(
    hierarchies;
    method = PhysicalContinuation(),
    reference = AutoReference(),
    parameter_path = ResetPath(),
    realizability = MomentMatrixChecks(),
)
```

A result should expose:

```julia
result.solutions
result.shared_moment_errors
result.stationary_residuals
result.order_lift_defects
result.adjoint_residuals
result.positivity_margins
result.stability_abscissae
result.path_ambiguity
result.status
```

## Method summary

### Physical parameter continuation

For a physical path `p(t)`, track

```math
F_n(x; p(t)) = 0
```

from a parameter point where the stationary state is known to the requested target point.

A useful regularized path adds a reset generator

```math
\mathcal L_{\lambda,\epsilon}
= \lambda \mathcal L_{\mathrm{target}}
+ \epsilon \mathcal R_\sigma,
```

where

```math
\mathcal R_\sigma(\rho)
= \operatorname{Tr}(\rho)\sigma - \rho.
```

The reference state is exactly `σ` when `λ = 0`. The reset strength can be removed after the target Liouvillian is activated.

### Continuation in hierarchy order

Moving from order `n` to `n + 1` introduces new retained moments `z`. Their start values must be the factorized values produced by the order `n` closure,

```math
z = \Phi_n(x),
```

not zero unless the numerical variables are connected cumulants.

Define closure deviation coordinates

```math
u = z - \Phi_n(x).
```

The order continuation starts at `u = 0` and gradually activates the full order `n + 1` equations.

### Physical validation

Each candidate stationary root will be evaluated using:

```math
r_n = \|F_n(x_n)\|,
```

```math
\delta_n = \|W F_{n+1}(x_n, \Phi_n(x_n))\|,
```

```math
e_{n,k} = \|\pi_k x_n - \pi_k x_{n-1}\|,
```

and moment matrix tests

```math
M_B(x)_{ij} = \langle B_i^\dagger B_j \rangle_x,
\qquad M_B(x) \succeq 0.
```

These quantities measure stationarity, held out hierarchy consistency, convergence of shared moments, and necessary quantum realizability.

## Scope

The first release should focus on one or a small number of stationary branches. It should not attempt to replace a complete numerical algebraic geometry package.

General root enumeration, singular endgames, monodromy, and certified complex algebraic geometry should initially be delegated to HomotopyContinuation.jl. Real folds and stability analysis may use BifurcationKit.jl.

The package contribution is the hierarchy aware construction, branch management, scaling, physical filtering, and diagnostics.

## Status

This repository currently contains a design specification. Implementation has not started.

See [docs/SPEC.md](docs/SPEC.md) for the technical specification and [docs/ROADMAP.md](docs/ROADMAP.md) for the proposed development sequence.

## License

MIT
