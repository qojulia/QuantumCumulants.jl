using QuantumCumulants
using Symbolics: Symbolics, @variables
using SymbolicUtils
using Test

# 6th-order cumulant closure of the damped JC system.

@testset "higher-order: 6th-order closure runs (no spectrum assertion)" begin
    @variables Δ::Real g::Real γ::Real κ::Real ν::Real
    hf = FockSpace(:cavity)
    ha = NLevelSpace(:atom, 2)
    h = hf ⊗ ha
    a = Destroy(h, :a)
    s = Transition(h, :σ, 1, 2)

    H = Δ * a' * a + g * (a' * s + a * s')
    J = [a, s, s']
    rates = [κ, γ, ν]

    # `phase_invariant` is the shared U(1) filter from runtests.jl `init_code`.
    he6 = complete(
        meanfield(a' * a, H, J; rates = rates);
        order = 6, filter_func = phase_invariant
    )
    @test length(he6.equations) > 0
    @test isempty(find_missing(he6; filter_func = phase_invariant))
end
