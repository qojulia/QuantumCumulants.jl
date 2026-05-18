# # Cavity Antiresonance (meanfield setup)

# Driven cavity with $N$ dipole-coupled atoms (cf. Plankensteiner et al. PRL 119,
# 093601 (2017)). The transmission antiresonance arises from collective decay.
# The full example also unrolls the indexed sums for a fixed $N$ via `evaluate`;
# in v1.0 the evaluate-with-limits pass is still under development, so this
# port stops at the symbolic equations.

using QuantumCumulants
using OrdinaryDiffEq, ModelingToolkitBase
using Plots

hc = FockSpace(:cavity)
ha = NLevelSpace(:atom, 2)
h = hc ⊗ ha

@variables N Δc η Δa κ Γ Ω
g(i) = IndexedVariable(:g, i)

i = Index(h, :i, N, ha)
j = Index(h, :j, N, ha)

@qnumbers a::Destroy(h)
σ(x, y, k) = IndexedOperator(Transition(h, :σ, x, y), k)

Hc = Δc * a' * a + η * (a' + a)
Ha = Δa * Σ(σ(2, 2, i), i)
Hi = Σ(g(i) * (a' * σ(1, 2, i) + a * σ(2, 1, i)), i)
H = Hc + Ha + Hi

J = [a, σ(1, 2, i)]
rates = [κ, Γ]

eqs = meanfield(a, H, J; rates = rates, order = 1)

println("Number of base equations: ", length(eqs.equations))
println("First equation: ", eqs.equations[1])
