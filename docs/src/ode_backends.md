# Solving `MeanfieldEquations` directly

Completed deterministic `MeanfieldEquations` can be lowered directly to a compact numerical
evaluator and passed to SciML without first constructing a ModelingToolkit system:

```julia
using QuantumCumulants
using SciMLBase: ODEProblem

eqs = complete(meanfield([a], H, [a]; rates = [κ], order = 2))
prob = ODEProblem(eqs, u0, (0.0, 10.0), Dict(κ => 0.5); backend = KernelBackend())
```

The direct path has one execution strategy:

```text
MeanfieldEquations → structured lowering → compact tables → ODEFunction
```

The structured representation is independent of the evaluator. Distinct state monomials are
shared globally and represented as signed state factors. A positive factor reads `u[j]`; a
negative factor reads `conj(u[j])`. The evaluator computes the monomial vector `v` and then
applies the sparse coefficient table `M`, giving `du = M * v`.

## Supported systems

The direct evaluator supports completed, deterministic cumulant equations whose drifts are
polynomial in the tracked moments with parameter-only coefficients. Ordinary unfolded systems,
conjugate-folded systems for RHS evaluation, and systems produced by `scale` or `evaluate` are
supported. The current numeric representation is specialized to `ComplexF64` state vectors and
coefficient values.

Parameters are supplied when the problem is constructed. The lowering and monomial structure
are independent of those values, so a sweep only refreshes the numeric coefficient table:

```julia
update_parameters!(prob, Dict(κ => 0.7))
```

Scalar parameters, evaluated one- and multi-dimensional array parameters, and repeated solves
are supported. The parameter update is equivalent to constructing a fresh problem with the new
values, without re-lowering the equations.

Non-polynomial rewrites, time-dependent coefficients, noise equations, and correlation systems
are outside this first direct path. They should use the existing ModelingToolkit route:

```julia
sys = System(eqs; name = :meanfield)
```

The direct constructors fail with a capability error for unsupported input. There is no silent
generated-code fallback.

## Jacobians

`jac = true` and `jac = :analytic` attach the sparse analytic Jacobian derived from the same
monomial tables:

```julia
prob = ODEProblem(eqs, u0, tspan, ps; backend = KernelBackend(), jac = true)
```

This mode is valid only for holomorphic closures. For a drift
``f(u, \bar u)`` the differential is ``δf = Aδu + Bδ\bar u``. A single complex `n × n`
Jacobian represents the differential only when `B == 0`. A folded system containing conjugate
state factors therefore raises `HolomorphicJacobianError`; it is never silently replaced by a
one-sided finite difference or a real-directional approximation. Use an unfolded holomorphic
closure or the ModelingToolkit path when that Jacobian is required.

## Direct solution access

The QC state registry is used for direct solutions as well as ModelingToolkit solutions:

```julia
sol = solve(prob, Tsit5())
get_solution(sol, a, eqs)(1.0)
get_solution(sol, a', eqs)(1.0)
```

The requested operator is resolved to a registered integer state index. If a folded query is the
conjugate side of the stored representative, `get_solution` returns the conjugated trajectory.
This keeps the user-facing `get_solution(sol, op, eqs)` behavior the same for both numerical
routes.
