```@meta
CurrentModule = QuantumCumulants
```

# Cumulant expansion

The Heisenberg equation for a product of ``n`` operators generally involves products of ``n+1`` operators, so the equations for averages never close on their own. The cumulant expansion cuts this tower: it assumes that all joint cumulants above a chosen order vanish and rewrites every higher moment through lower ones. The [theory page](@ref theory) derives the expansion, and [Convergence and validity](@ref convergence) discusses when the truncation can be trusted.

```@contents
Pages = ["cumulant_expansion.md"]
Depth = 2:3
```

## Order of an average

The order of an average is the number of operator factors it contains. It decides which averages the expansion rewrites and which it keeps.

```@example cumulant
using QuantumCumulants

h = FockSpace(:cavity) ⊗ NLevelSpace(:atom, (:g, :e))
@qnumbers a::Destroy(h)
σ(i, j) = Transition(h, :σ, i, j)

get_order(average(a' * a * σ(:e, :e)))
```

```@docs
get_order
```

## Joint cumulants

The joint cumulant of a product measures the correlation that is not already contained in its factors. For three operators it is the third moment minus every way of splitting it into smaller moments.

```@example cumulant
cumulant(a' * a * σ(:e, :e))
```

Setting this cumulant to zero yields the second-order expansion of the third moment.

```@example cumulant
cumulant_expansion(average(a' * a * σ(:e, :e)), 2)
```

```@docs
cumulant
```

## Expanding equations

[`cumulant_expansion`](@ref) acts on single averages, on arbitrary expressions of averages, and on a whole [`MeanfieldEquations`](@ref). The `order` keyword of [`meanfield`](@ref) applies it during the derivation, which is the usual route. Expanding afterwards is useful to compare several orders on the same raw equations.

```@example cumulant
@variables Δ::Real g::Real κ::Real
H = Δ * a' * a + g * (a' * σ(:g, :e) + a * σ(:e, :g))
eqs = meanfield([a' * σ(:g, :e)], H, [a]; rates = [κ])
cumulant_expansion(eqs, 1)
```

A hybrid system may need different orders on different subspaces, e.g. a classical-like cavity coupled to strongly correlated atoms. A vector `order` sets one cap per Hilbert subspace, and `mix_choice` decides the cap of a term that spans several.

```@docs
cumulant_expansion
```
