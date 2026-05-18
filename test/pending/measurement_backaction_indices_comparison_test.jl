# PENDING: port of master test/test_measurement_backaction_indices_comparison.jl
#
# Status: blocked on the same SQA index-conflict issue
# (`commutator(im*H, op)` with collective Σ-Hamiltonian and indexed ops),
# AND on `scale(::NoiseMeanFieldEquations)` (TODO.md: not defined in v1).
# Also depends on `cumulant_expansion(::NoiseMeanFieldEquations, k)`
# rebuilding the cumulant-truncated noise equations, which doesn't have a
# v1 surface yet.
#
# To enable: implement `scale(::NoiseMeanFieldEquations)` (drift + noise
# both scaled), `cumulant_expansion(::NoiseMeanFieldEquations, _)`, then
# resolve the SQA index-conflict guard for collective Σ-Hamiltonians.

#=
using QuantumCumulants
using Symbolics: Symbolics, @variables, simplify, expand, @syms, @register_symbolic
using SymbolicUtils
using Test

@testset "measurement backaction indexed: scaled stoch vs det at order 4" begin
    @variables N::Real ωa::Real γ::Real η::Real χ::Real ωc::Real κ::Real
    @variables g::Real ξ::Real ωl::Real
    @syms t::Real
    @register_symbolic pulse(t)

    hc = FockSpace(:resonator); ha = NLevelSpace(:atom, 2)
    h = hc ⊗ ha
    j = Index(h, :j, N, ha)
    k = Index(h, :k, N, ha)
    @qnumbers a::Destroy(h, 1)
    σ(α, β, kk) = IndexedOperator(Transition(h, :σ, α, β, 2), kk)

    H = ωc * a' * a + ωa * Σ(σ(2, 2, j), j) +
        g * a' * Σ(σ(1, 2, j), j) + g * a * Σ(σ(2, 1, j), j)
    J = [a * exp(1.0im * ωl * t), σ(1, 2, j), σ(2, 1, j), σ(2, 2, j)]
    rates = [κ, γ, η * pulse(t), 2 * χ]
    efficiencies = [ξ, 0, 0, 0]
    ops = [a, a' * a, σ(2, 2, k), σ(1, 2, k)]

    stoch = meanfield(ops, H, J;
                      rates = rates, efficiencies = efficiencies, order = 2)
    det = meanfield(ops, H, J; rates = rates, order = 2)

    # Stochastic drift must equal deterministic drift (efficiencies enter
    # the noise channel, not the drift).
    for (s, d) in zip(stoch.equations, det.equations)
        @test isequal(simplify(s.lhs - d.lhs), 0)
        @test isequal(simplify(s.rhs - d.rhs), 0)
    end

    stoch_sc = scale(stoch)
    det_sc = scale(det)
    for (s, d) in zip(stoch_sc.equations, det_sc.equations)
        @test isequal(simplify(s.lhs - d.lhs), 0)
        @test isequal(simplify(expand(s.rhs - d.rhs)), 0)
    end

    # Order-4 vs order-2 cumulant-collapse on the scaled stochastic system
    # should reproduce the scaled order-2 noise equations.
    full_stoch = meanfield(stoch.operators, H, J;
                            rates = rates, efficiencies = efficiencies, order = 4)
    stoch_complete = complete(stoch)
    full_stoch_2 = cumulant_expansion(full_stoch, 2)
    for (lhs, rhs) in zip(stoch_complete.noise_equations, full_stoch_2.equations)
        @test isequal(simplify(lhs.lhs - rhs.lhs), 0)
        @test isequal(simplify(lhs.rhs - rhs.rhs), 0)
    end
end
=#
