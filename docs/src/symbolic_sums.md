# Symbolic Sums and Indices

Many physical systems contain multiple elements of the same kind, which basically do the same thing just with different rates. For these systems it is convenient to describe the Hamiltonian and the dissipative processes with indexed objects and sums. A well-known example is the Tavis-Cummings Hamiltonian, which describes the interaction of $N$ two-level atoms with a cavity mode according to the Hamiltonian

```math
\begin{equation}
H_\mathrm{TC} = \omega_c a^† a + \sum_i^N \omega_i \sigma_i^{22} + \sum_i^N g_i (a^\dagger \sigma_i^{12} + a \sigma_i^{21}).
\end{equation}
```

In principle we can write down and derive the equations for all $N$ atoms explicitly, but this can take a long time for large $N$. The more practical and elegant approach is to derive the equations for averages of indexed operators and insert all possible number combinations afterwards. The implementation of symbolic sums and indices allows for exactly this.

## Implementation

### Index

The main tool to use symbolic summations is the [`Index`](@ref) object. An index is constructed from the full [`HilbertSpace`](@ref SecondQuantizedAlgebra.HilbertSpace) `h`, a `name` (a `Symbol`), a `range` (a `Symbol`, a concrete `Int`, or a symbolic `Num`), and the specific subspace `space` (or its integer `space_index`) the index acts on. This means an index for a [`NLevelSpace`](@ref) can only be used by [`Transition`](@ref) operators. In the example below, two indices are defined equivalently, and a third one is defined on the [`FockSpace`](@ref) of the same [`ProductSpace`](@ref) `h`.


```@example symbolic_sums
using QuantumCumulants

@variables N::Int

ha = NLevelSpace(:atoms, 2)
hc = FockSpace(:cavity)
h = hc ⊗ ha

i  = Index(h, :i, N, ha)
i2 = Index(h, :i, N, 2) # equivalent: ha is the second subspace
n  = Index(h, :n, 5, hc)
```


### IndexedOperators

Operators like [`Destroy`](@ref) or [`Transition`](@ref) can be tied to an [`Index`](@ref) of the corresponding subspace by wrapping them in an [`IndexedOperator`](@ref). The constructor takes the underlying operator and the index. Below, indexed callables are built for the two subspaces.


```@example symbolic_sums
σ(x, y, z) = IndexedOperator(Transition(h, :σ, x, y), z)
a(z) = IndexedOperator(Destroy(h, :a), z)
```


We can now form products that already carry their per-site information:


```@example symbolic_sums
a(n) * σ(2, 2, i)
nothing #hide
```

```math
{a}_{n} {\sigma}_{i}^{{22}}
```

Symbolic per-site coefficients are expressed via [`IndexedVariable`](@ref):


```@example symbolic_sums
gi = IndexedVariable(:g, i)
nothing #hide
```

```math
{g}_{i}
```

### Summations

Indexed operators (and indexed variables) compose into symbolic summations. A sum needs the `term` to be summed and the [`Index`](@ref) it runs over. The simplest case is a single-index sum:

```@example symbolic_sums
∑(σ(2, 2, i), i)
nothing #hide
```

```math
\underset{i}{\overset{N}{\sum}} {σ}_{i}^{{22}}
```

The constructor is available as both `∑` (`\sum`) and `Σ` (`\Sigma`). A trailing vector of indices marks them as non-equal to the running index:


```@example symbolic_sums
j = Index(h, :j, N, ha)
∑(σ(2, 2, i), i, [j])
nothing #hide
```

```math
\underset{i ≠j }{\overset{N}{\sum}} {σ}_{i}^{{22}}
```

Multi-index sums are written by nesting the constructor:

```@example symbolic_sums
∑(∑(a(n) * σ(2, 1, i), i), n)
nothing #hide
```

```math
\underset{i}{\overset{N}{\sum}} \underset{n}{\overset{5}{\sum}} {a}_{n}  {σ}_{i}^{{21}}
```


When two running indices act on the same subspace, the `i = j` diagonal slice is split out automatically:


```@example symbolic_sums
k = Index(h, :k, N, ha)
l = Index(h, :l, N, ha)

∑(∑(σ(2, 1, k) * σ(1, 2, l), k), l)
nothing #hide
```

```math
\underset{k{\ne}l}{\overset{N}{\sum}} \underset{l{\ne}k}{\overset{N}{\sum}} {\sigma}_{l}^{{12}}  {\sigma}_{k}^{{21}} + \underset{k}{\overset{N}{\sum}} {\sigma}_{k}^{{22}}
```


The same diagonal split fires when a sum is multiplied by an [`IndexedOperator`](@ref) on the same subspace:


```@example symbolic_sums
∑(σ(2, 2, k), k) * σ(2, 1, l)
nothing #hide
```

```math
\underset{k{\ne}l}{\overset{N}{\sum}} {\sigma}_{k}^{{22}}  {\sigma}_{l}^{{21}} + {\sigma}_{l}^{{21}}
```

## Short Example

We walk through the full pipeline (Hamiltonian to equations to numeric solution) for `N` two-level atoms in a single-mode cavity.

Set up indices and operators, then write the Hamiltonian:


```@example symbolic_sums
using QuantumCumulants

ha = NLevelSpace(:atoms, 2)
hc = FockSpace(:cavity)
h = hc ⊗ ha

@variables N::Int Δ::Real κ::Real γ::Real ν::Real

i = Index(h, :i, N, ha)
j = Index(h, :j, N, ha)

@qnumbers b::Destroy(h)
σ(x, y, z) = IndexedOperator(Transition(h, :σ, x, y), z)
gi = IndexedVariable(:g, i)

H = Δ*b'*b + ∑(gi*(b*σ(2, 1, i) + b'*σ(1, 2, i)), i)
nothing #hide
```

```math
\underset{i}{\overset{N}{\sum}} {g}_{i}  b  {\sigma}_{i}^{{21}} + \underset{i}{\overset{N}{\sum}} {g}_{i}  b^\dagger  {\sigma}_{i}^{{12}} + \Delta b^\dagger b
```

The operators we derive equations for and the jump operators are specified next. Indexed operators used as the derivation targets must carry an index distinct from any already bound in the Hamiltonian or jump terms. An indexed jump operator $J_i$ contributes the dissipator

```math
\begin{equation}
\frac{1}{2} \sum_{i} R_{i} \bigg( 2 J_i^\dagger \mathcal{O} J_i - \mathcal{O} J_i^\dagger J_i -  J_i^\dagger J_i \mathcal{O} \bigg).
\end{equation}
```

The rate may be scalar or an indexed variable; if indexed, its index must match the operator's.


```@example symbolic_sums
J     = [b, σ(1, 2, i), σ(2, 1, i)]
rates = [κ, γ, ν]

eqs = meanfield(b'b, H, J; rates=rates, order=2)
nothing #hide
```

Closing the system with [`complete`](@ref) derives equations for any averages appearing on the RHS that are missing from the LHS:

```@example symbolic_sums
eqs_comp = complete(eqs)
nothing #hide
```

### Evaluate and Scale

The closed equations still contain symbolic sums and the symbolic upper bound `N`. There are two routes to concrete numeric equations: [`evaluate`](@ref) unrolls each sum into `N` per-site equations, and [`scale`](@ref) collapses permutation-equivalent terms by assuming the atoms are identical. `scale` typically produces far fewer equations than `evaluate` because it exploits the permutation symmetry; `evaluate` is the right call when the atoms differ (different couplings, different rates).

Here we unroll for $N=3$ atoms with `evaluate`. The numeric value of $N$ is passed via the `limits` keyword:

```@example symbolic_sums
evaled = evaluate(eqs_comp; limits=(N => 3))
nothing #hide
```

The unrolled equations are now ready to feed into `System`, which builds a [`ModelingToolkitBase.System`](https://github.com/SciML/ModelingToolkitBase.jl) that can be solved with [OrdinaryDiffEq](https://github.com/SciML/OrdinaryDiffEq.jl). The per-atom coupling `g_i` is a vector parameter; [`parameter_map`](@ref) takes a dict of symbolic parameters (scalar or array-valued) and produces a parameter substitution dict that matches the compiled system's array-shaped parameters.


```@example symbolic_sums
using ModelingToolkitBase
sys = mtkcompile(System(evaled; name=:tc))

using OrdinaryDiffEq
u0 = zeros(ComplexF64, length(evaled.states))
p  = parameter_map(evaled, Dict(
    Δ  => 0.0,
    gi => [0.75, 1.2, 1.5],
    γ  => 0.25,
    κ  => 1.0,
    ν  => 1.5,
))
prob = ODEProblem(sys, merge(initial_values(evaled, u0), p), (0.0, 10.0))
sol = solve(prob, Tsit5())
nothing #hide
```

Use [`get_solution`](@ref) to evaluate any operator-average trajectory. After `evaluate`, the per-atom excited-state averages live as concrete entries in `evaled.states`; we filter them out and plot them alongside the photon number:


```@example symbolic_sums
using Plots
using SecondQuantizedAlgebra: undo_average

ts = range(0.0, 10.0; length=200)
n  = real.(get_solution(sol, b'*b, evaled).(ts))

# Per-atom excited-state averages produced by `evaluate(eqs; limits=(N => 3))`.
pe_states = [s for s in evaled.states if string(undo_average(s)) |> contains("σ") &&
                                          string(undo_average(s)) |> contains("22")]

pl = plot(ts, n, label="Photon number", xlabel="t")
for (k, s) in enumerate(pe_states)
    plot!(pl, ts, real.(get_solution(sol, undo_average(s), evaled).(ts)),
        label="Excited state population of atom $(k)")
end
pl # hide
```
