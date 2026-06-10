# Indexed and scaled systems

Many physical systems contain multiple identical elements: ``N`` atoms in a cavity, a chain of emitters, a register of qubits. Rather than writing ``N`` copies of every equation, you write the Hamiltonian once with *indexed* operators and symbolic sums, derive the equations in terms of a running index, and only commit to a concrete ``N`` at the numerical stage.

The indexed-operator machinery itself ([`Index`](@ref), [`IndexedOperator`](@ref), [`IndexedVariable`](@ref), the summation constructor [`Σ`](@ref) (also written `∑`), and the automatic diagonal splitting of products) is provided by [SecondQuantizedAlgebra.jl](https://github.com/qojulia/SecondQuantizedAlgebra.jl); see its [Symbolic Sums and Indices](https://qojulia.github.io/SecondQuantizedAlgebra.jl/stable/symbolic_sums/) guide for how to build indexed expressions. This page covers what QuantumCumulants adds on top: deriving, closing, and **collapsing or unrolling** an indexed mean-field system.

## A worked example: Tavis-Cummings

We take ``N`` two-level atoms in a single-mode cavity. The atom subspace carries an [`Index`](@ref) `i`, the per-atom coupling is an [`IndexedVariable`](@ref), and the Hamiltonian sums over the atoms with [`Σ`](@ref):

```@example sums
using QuantumCumulants

ha = NLevelSpace(:atoms, 2)
hc = FockSpace(:cavity)
h = hc ⊗ ha

@variables N::Int Δ::Real κ::Real γ::Real ν::Real

i = Index(h, :i, N, ha)

@qnumbers b::Destroy(h)
σ(x, y, z) = IndexedOperator(Transition(h, :σ, x, y), z)
gi = IndexedVariable(:g, i)

H = Δ*b'*b + ∑(gi*(b*σ(2, 1, i) + b'*σ(1, 2, i)), i)
nothing # hide
```

```math
\underset{i}{\overset{N}{\sum}} {g}_{i}  b  {\sigma}_{i}^{{21}} + \underset{i}{\overset{N}{\sum}} {g}_{i}  b^\dagger  {\sigma}_{i}^{{12}} + \Delta b^\dagger b
```

[`meanfield`](@ref) and [`complete`](@ref) work exactly as in the scalar case; an indexed jump operator ``J_i`` with rate ``R_i`` contributes the dissipator

```math
\frac{1}{2} \sum_{i} R_{i} \left( 2 J_i^\dagger \mathcal{O} J_i - \mathcal{O} J_i^\dagger J_i -  J_i^\dagger J_i \mathcal{O} \right),
```

with a scalar or matching-index rate:

```@example sums
J     = [b, σ(1, 2, i), σ(2, 1, i)]
rates = [κ, γ, ν]

eqs = meanfield(b'b, H, J; rates=rates, order=2)
eqs_comp = complete(eqs)
nothing # hide
```

## Evaluate vs. scale

The closed equations still contain symbolic sums and the symbolic bound `N`. There are two routes to concrete numeric equations:

- [`evaluate`](@ref) unrolls each sum into `N` per-site equations. Use it when the atoms *differ* (different couplings or rates).

- [`scale`](@ref) collapses permutation-equivalent terms by assuming the atoms are *identical*, typically yielding far fewer equations.

Both accept an `h::Vector{Int}` of subspace `space_index` values to target specific Hilbert factors, so a hybrid system can unroll some subspaces and collapse others.

Here we unroll for ``N = 3`` atoms; the numeric value of `N` is passed via `limits`:

```@example sums
evaled = evaluate(eqs_comp; limits=(N => 3))
nothing # hide
```

## Numerical solution

The unrolled (or scaled) equations feed into `System` like any other. A per-atom coupling `g_i` is a vector parameter; [`parameter_map`](@ref) turns a dict of symbolic parameters (scalar or array-valued) into the substitution dict the compiled system expects:

```@example sums
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
nothing # hide
```

[`get_solution`](@ref) evaluates any operator-average trajectory. After `evaluate`, the per-atom excited-state averages are concrete entries in `evaled.states`:

```@example sums
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
