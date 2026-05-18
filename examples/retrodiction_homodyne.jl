# # Retrodiction with Homodyne Detection (backward-noise meanfield setup)

# Backward (retrodicted) noise meanfield equations using
# `direction = Backward()`. Full SDE integration with `SDESystem` is still in
# development; this example shows the equation construction.

using QuantumCumulants
using ModelingToolkitBase

hc = FockSpace(:cavity)
@qnumbers a::Destroy(hc)

@variables ω κ η

H = ω * a' * a

eqs = meanfield(
    [a, a' * a], H, [a];
    rates = [κ], efficiencies = [η],
    direction = Backward(), order = 2
)

@assert eqs isa NoiseMeanFieldEquations
@assert eqs.direction isa Backward
println(
    "Built backward NoiseMeanFieldEquations: ",
    length(eqs.equations), " drift / ", length(eqs.noise_equations), " noise."
)
