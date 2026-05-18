# PENDING: port of master test/test_measurement_backaction_indices.jl
#
# Status: blocked on SQA v0.5 strict index-conflict guard during
# `commutator(im * H, op)` for collective JC Hamiltonians (the same
# "Summation index appears in both factors" error noted in TODO.md). Once
# the indexed `meanfield` pipeline accepts collective Σ-Hamiltonians, this
# port will exercise the indexed measurement-backaction equality:
#       eqs_c (deterministic).equations == eqs_c_noise (with efficiencies).equations
# (the noise enters `noise_equations` only; drift must agree term-for-term).
#
# The non-indexed bare-`a` portion is already covered by
# `test/noise_test.jl::"measurement backaction: full drift + noise match
# analytic"`. The skeleton below covers only the indexed part.

#=
using QuantumCumulants
using Symbolics: Symbolics, @variables, simplify
using SymbolicUtils
using Test

@testset "measurement backaction (indexed): drift agrees det vs stoch" begin
    @variables κ::Real g::Real gf::Real κf::Real R::Real Γ::Real
    @variables Δ::Real ν::Real N::Real M::Real η::Real
    δ(i) = IndexedVariable(:δ, i)

    hc = FockSpace(:cavity); hf = FockSpace(:filter); ha = NLevelSpace(:atom, 2)
    h = hc ⊗ hf ⊗ ha
    i = Index(h, :i, M, hf)
    j = Index(h, :j, N, ha)

    @qnumbers a::Destroy(h, 1)
    b(k) = IndexedOperator(Destroy(h, :b, 2), k)
    σ(α, β, k) = IndexedOperator(Transition(h, :σ, α, β, 3), k)

    H = Δ * Σ(σ(2, 2, j), j) +
        Σ(δ(i) * b(i)' * b(i), i) +
        gf * Σ(a' * b(i) + a * b(i)', i) +
        g * Σ(a' * σ(1, 2, j) + a * σ(2, 1, j), j)

    J = [a, b(i), σ(1, 2, j), σ(2, 1, j), σ(2, 2, j)]
    rates = [κ, κf, Γ, R, ν]
    efficiencies = [η, 0, 0, 0, 0]

    eqs_det = meanfield(a' * a, H, J; rates = rates, order = 2)
    eqs_det_c = complete(eqs_det)
    eqs_noise = meanfield(a' * a, H, J;
                          rates = rates, efficiencies = efficiencies, order = 2)
    eqs_noise_c = complete(eqs_noise)

    @test length(eqs_det_c.equations) == length(eqs_noise_c.equations)
    for (eq, eq_noise) in zip(eqs_det_c.equations, eqs_noise_c.equations)
        @test isequal(eq.lhs, eq_noise.lhs)
        @test isequal(simplify(eq.rhs - eq_noise.rhs; expand = true), 0)
    end
end
=#
