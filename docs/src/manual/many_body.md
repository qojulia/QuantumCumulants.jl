```@meta
CurrentModule = QuantumCumulants
```

# Many-body systems

An ensemble of ``N`` emitters would naively need ``N`` copies of every equation. With indexed operators the Hamiltonian is written once as a symbolic sum, the equations are derived for a running index, and ``N`` stays a symbol until the numerical stage. Building indexed operators and sums is covered by the [many-body manual](https://qojulia.github.io/SecondQuantizedAlgebra.jl/stable/manual/many_body/) of SecondQuantizedAlgebra. This page covers the two ways QuantumCumulants turns the resulting symbolic system into a finite one.

```@contents
Pages = ["many_body.md"]
Depth = 2:3
```

## Indexed equations

We take ``N`` two-level atoms with individual couplings ``g_i`` in a single-mode cavity, the Tavis-Cummings model with incoherent pumping.

```@example many
using QuantumCumulants

ha = NLevelSpace(:atoms, 2)
hc = FockSpace(:cavity)
h = hc ⊗ ha

@variables N::Int Δ::Real κ::Real γ::Real ν::Real
i = Index(h, :i, N, ha)

@qnumbers b::Destroy(h)
σ(x, y, z) = IndexedOperator(Transition(h, :σ, x, y), z)
gi = IndexedVariable(:g, i)

H = Δ * b' * b + ∑(gi * (b * σ(2, 1, i) + b' * σ(1, 2, i)), i)
J = [b, σ(1, 2, i), σ(2, 1, i)]
rates = [κ, γ, ν]

eqs = complete(meanfield(b' * b, H, J; rates, order = 2))
nothing # hide
```

An indexed jump with a scalar rate damps every atom independently. The closed equations still contain sums over `i` and the symbolic bound `N`.

## Unrolling or collapsing

There are two routes to a finite system, and the physics decides between them.

- [`evaluate`](@ref) unrolls every sum into ``N`` per-site equations. It is exact for any ``N`` and allows the atoms to differ, at the price of a system that grows with ``N``.
- [`scale`](@ref) assumes the atoms are identical and keeps one representative per symmetry class. The number of equations is independent of ``N``, so ``N`` can stay symbolic and be large.

Here the couplings ``g_i`` differ, so we unroll three atoms.

```@example many
evaled = evaluate(eqs; limits = (N => 3))
length(states(evaled))
```

Both functions take an `h` keyword that selects the Hilbert subspaces to act on. A hybrid system can thus unroll one ensemble and collapse another.

```@docs
evaluate
scale
scale!
```

## Parameters of an unrolled system

After [`evaluate`](@ref), the indexed coupling ``g_i`` is an array parameter of length ``N``. [`parameter_map`](@ref) translates the user-facing key into it, so values are given per atom.

```@example many
using ModelingToolkitBase, OrdinaryDiffEqTsit5

sys = mtkcompile(System(evaled; name = :tc))
ps = parameter_map(evaled, Dict(Δ => 0.0, gi => [0.75, 1.2, 1.5], γ => 0.25, κ => 1.0, ν => 1.5))
prob = ODEProblem(sys, merge(initial_values(evaled), ps), (0.0, 10.0))
sol = solve(prob, Tsit5())

using Plots
ts = range(0.0, 10.0; length = 200)
plot(ts, real.(get_solution(sol, b' * b, evaled).(ts)); label = "photon number", xlabel = "t")
```
