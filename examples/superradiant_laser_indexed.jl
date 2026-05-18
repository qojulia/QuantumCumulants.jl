# # Superradiant Laser

# A superradiant laser of $N$ identical two-level atoms coupled to a single cavity
# mode (cf. Meiser et al. PRL 102, 163601 (2009)). The Hamiltonian is
#
# $H = -\Delta a^\dagger a + \sum_{j=1}^{N} g_j (a^\dagger \sigma_j^{12} + a \sigma_j^{21}).$

using QuantumCumulants
using OrdinaryDiffEq, ModelingToolkitBase
using Plots

hc = FockSpace(:cavity)
ha = NLevelSpace(:atom, 2)
h = hc ⊗ ha

@qnumbers a::Destroy(h)
σ(α, β, i) = IndexedOperator(Transition(h, :σ, α, β), i)

@variables N Δ κ Γ R ν
g(i) = IndexedVariable(:g, i)

i = Index(h, :i, N, ha)
j = Index(h, :j, N, ha)

H = -Δ * a' * a + Σ(g(i) * (a' * σ(1, 2, i) + a * σ(2, 1, i)), i)

J = [a, σ(1, 2, i), σ(2, 1, i), σ(2, 2, i)]
rates = [κ, Γ, R, ν]

ops = [a' * a, σ(2, 2, j)]
eqs = meanfield(ops, H, J; rates = rates, order = 2, simplify = false)

# Note: in v1.0, the deeper `complete!` for this system needs phase-invariant
# filtering and additional indexed-correlation infrastructure that is still in
# development. This example currently demonstrates the equation set and the
# `to_system` boundary; downstream solver integration will follow.

println("Generated ", length(eqs.equations), " base meanfield equations.")
println("First equation: ", eqs.equations[1])
