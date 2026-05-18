# # Unique Steady-State Squeezing (meanfield setup)

# Driven Dicke model (Gietka et al. PRL 131, 223604 (2023)).
# The full example requires sums whose summation index also appears in a
# multiplied operator at the same site; SQA v0.5 currently rejects this with
# "Summation index appears in both factors". The example below shows the bare
# meanfield setup for the simplest single-spin case while that algebraic
# limitation is being addressed upstream.

using QuantumCumulants
using OrdinaryDiffEq, ModelingToolkitBase
using Plots

hf = FockSpace(:harmonic)
ha = NLevelSpace(:spin, 2)
h = hf ⊗ ha

@variables ω Ω η κ g γ ξ

@qnumbers a::Destroy(h)
σ(x, y) = Transition(h, :σ, x, y)

b = a * cosh(ξ) + a' * sinh(ξ)

H = ω * a' * a + Ω / 2 * (σ(2, 2) - σ(1, 1)) +
    g * (σ(1, 2) + σ(2, 1)) * (a + a') / 2

J = [b, σ(1, 2)]
rates = [κ, γ]

eqs = meanfield(
    [a, a' * a, σ(2, 2)], H, J;
    rates = rates, order = 2, simplify = false
)
println("Generated ", length(eqs.equations), " equations.")
println("First eq: ", eqs.equations[1])
