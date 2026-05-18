# # Heterodyne detection (forward-noise meanfield setup)

# Stochastic-master-equation framework with measurement backaction
# (Yu et al. PRL 133, 073601 (2024)). v1.0 supports forward-noise meanfield
# equations via `meanfield(..., efficiencies = ..., direction = Forward())`;
# the full `SDESystem` integration is still in development, so this example
# stops at constructing the `NoiseMeanFieldEquations`.

using QuantumCumulants
using ModelingToolkitBase

hc = FockSpace(:resonator)
ha = NLevelSpace(:atom, 2)
h = hc ⊗ ha

@variables ωa γ η χ ωc κ g ξ

@qnumbers a::Destroy(h, 1)
σ(α, β) = Transition(h, :σ, α, β, 2)

H = ωc * a' * a + ωa * σ(2, 2) + g * a' * σ(1, 2) + g * a * σ(2, 1)

J = [a, σ(1, 2), σ(2, 1), σ(2, 2)]
rates = [κ, γ, η, 2 * χ]
efficiencies = [ξ, 0, 0, 0]

ops = [a', a' * a, σ(2, 2), σ(1, 2), a * a]
eqs = meanfield(
    ops, H, J;
    rates = rates, efficiencies = efficiencies,
    direction = Forward(), order = 2, simplify = false
)

@assert eqs isa NoiseMeanFieldEquations
@assert eqs.direction isa Forward
println(
    "Built NoiseMeanFieldEquations: ", length(eqs.equations), " drift / ",
    length(eqs.noise_equations), " noise equations."
)
