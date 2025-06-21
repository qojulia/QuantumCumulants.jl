using QuantumCumulants
using SymbolicUtils
using Test
using Symbolics

@testset "test_measurement_backaction" begin

    @cnumbers ω κ η

    hc = FockSpace(:cavity)
    a = Destroy(hc, :a)

    eqs = meanfield(a, ω * a' * a, [a]; rates = [κ], efficiencies = [η])
    test_eqn = sqrt(κ*η)*(average(a'*a+a * a)-average(a')*average(a)-average(a) ^ 2)

    @test isequal(simplify(test_eqn-eqs.noise_equations[1].rhs, expand = true), 0)
    @test isequal(eqs[1].lhs, average(a))

    @cnumbers ω κ η γ ωa

    ha = NLevelSpace(:atom, 2)
    hc = FockSpace(:cavity)
    h = ha ⊗ hc
    a = Destroy(h, :a)
    σ = Transition(h, :σ, 2, 1)

    eqs = meanfield(
        [a, σ],
        ω * a' * a + ωa * σ' * σ,
        [a];
        rates = [κ, γ],
        efficiencies = [0, 0],
    )

    @test isequal(eqs.noise_equations[1].rhs, 0)
    @test isequal(eqs.noise_equations[2].rhs, 0)

end
