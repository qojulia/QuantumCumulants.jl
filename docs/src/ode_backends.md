# Solving the equations directly (RHS backends)

QuantumCumulants can compile the completed moment equations straight into a
`SciMLBase.ODEProblem`, without constructing a `ModelingToolkit.System` first:

```julia
eqs = meanfield(ops, H, J; rates, order = 2)
complete!(eqs)
prob = ODEProblem(eqs, ψ0, (0.0, 10.0), Dict(Δ => 1.0, κ => 0.5))
sol = solve(prob, RK4())
```

`u0` can be a numeric quantum state (a ket or density operator from
QuantumOpticsBase), a `ComplexF64` vector aligned with `eqs.states`, or a `Dict`
keyed by the state averages. Read results with the same
[`get_solution`](@ref) used on the ModelingToolkit path.

## Backend selection

The `backend` keyword picks the RHS compilation strategy:

| Backend | What it is | When it wins |
|---|---|---|
| `AutoBackend()` (default) | tries the kernel, falls back to sharded exactly when the drift is not polynomial in the moments | almost always the right choice |
| `KernelBackend(; cache, parallel)` | the moment-polynomial kernel `du = M * v`: the drifts are lowered once to sparse data (CSR layout; both RHS passes optionally threaded via Polyester), no per-model native code | cold construction everywhere; warm runtime beyond roughly a thousand equations; analytic Jacobian; caching |
| `ShardedBackend(; chunk, threads, parallel)` | chunked native codegen (`build_function` per chunk of equations, codegen parallelized with OhMyThreads; the chunks also run concurrently at solve time for large systems) | warm RHS runtime on small and mid-size systems; t-dependent or rewritten non-polynomial drifts |

Both backends resolve the drifts through the same treatments-aware machinery as
`System(eqs)`, so scaled (`scale`) and evaluated (`evaluate`) systems, indexed
operators, and array-valued parameters work identically on either.

## Parameters

Parameters are required at construction and baked into the compiled RHS; there is no
symbolic parameter object at solve time. For sweeps and ensembles, rewrite the
compiled problem in place:

```julia
update_parameters!(prob, Dict(κ => 2.0))       # unmentioned parameters keep their values
```

`remake(prob, p = ...)` with a plain vector does not work on the kernel backend (its
`prob.p` is a `KernelParameters` object) and raises a typed error pointing back to
`update_parameters!`.

For ensembles, `copy` the `ODEFunction` so each trajectory owns its tables:

```julia
f2 = copy(prob.f)
prob2 = remake(prob; f = f2, p = f2.f.kp)
update_parameters!(prob2, Dict(κ => 3.0))       # prob is unaffected
```

## Analytic Jacobian (`jac = true`)

The kernel backend can attach the analytic sparse Jacobian, enabling implicit
(stiff) solvers on the complex-valued state without ForwardDiff:

```julia
prob = ODEProblem(eqs, ψ0, tspan, ps; jac = true)
sol = solve(prob, Rodas5P(autodiff = AutoFiniteDiff()))
```

The Jacobian is holomorphic-only: a system closed with `get_adjoints = false`
(conjugate-folded) raises `HolomorphicJacobianError`, because its true derivative
needs the Wirtinger pair. `ShardedBackend` has no analytic Jacobian and rejects
`jac = true`.

## Kernel cache

Lowering large systems takes seconds; the kernel tables can be persisted:

```julia
using JLD2                                       # activates the cache extension
kb = KernelBackend(cache = "my_kernels.jld2")
prob = ODEProblem(eqs, u0, tspan, ps; backend = kb)   # first run lowers and stores
prob = ODEProblem(eqs, u0, tspan, ps; backend = kb)   # later runs load the tables
```

The cache key is a SHA-256 digest of a canonical text of the equations plus version
stamps, and every hit re-verifies the stored text byte-for-byte, so a stale or
colliding entry silently degrades to a fresh lowering. Corrupt or unreadable files
warn and fall back; they never error. A loaded kernel is fully sweepable with
`update_parameters!`. With `jac = true` the cache is bypassed (the stored tables do
not carry the Jacobian's complement monomials).

## Error taxonomy

All lowering failures are typed (subtypes of `KernelLoweringError`) with guidance in
the message:

| Error | Meaning | Fix |
|---|---|---|
| `NonPolynomialDriftError` | a rewritten drift is not polynomial in the moments | use `ShardedBackend` (what `AutoBackend` does automatically) |
| `TimeDependentCoefficientError` | a coefficient depends on `t` | use `ShardedBackend()` or the MTK path |
| `ImParameterCollisionError` | a user parameter named `im` | rename the parameter |
| `UnresolvedMomentError` | an RHS average is not among the states | `complete!(eqs)` first |
| `HolomorphicJacobianError` | `jac = true` on a conjugate-folded system | close with `get_adjoints = true` or drop `jac` |

## What stays on the ModelingToolkit path

Events and callbacks defined at the `System` level, observed variables, symbolic
indexing of solutions, [`CorrelationFunction`](@ref)/[`Spectrum`](@ref), and
noise/SDE systems (`NoiseMeanfieldEquations`) are not handled by the direct
backends in this version. Build those through `System(eqs; name = ...)` as before;
the two routes coexist and agree.
