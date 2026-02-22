using QuantumCumulants
using Test
using SymbolicUtils
using Symbolics

@testset "test_measurement_backaction_indices" begin

    # Parameters
    @cnumbers κ g gf κf R Γ Δ ν N M η
    δ(i) = IndexedVariable(:δ, i)

    # Hilbertspace
    hc = FockSpace(:cavity)
    hf = FockSpace(:filter)
    ha = NLevelSpace(:atom, 2)
    h = hc ⊗ hf ⊗ ha

    # Indices and Operators
    i = Index(h, :i, M, hf)
    j = Index(h, :j, N, ha)

    @qnumbers a::Destroy(h, 1)
    b(k) = IndexedOperator(Destroy(h, :b, 2), k)
    σ(α, β, k) = IndexedOperator(Transition(h, :σ, α, β, 3), k)

    # test meanfield - measurement backaction
    test_eq_a = -1im*Δ*average(a) - 0.5κ*average(a)
    test_noise_eq_a =
        √(η*κ)*(average(a'a) + average(a*a) - average(a)^2 - average(a)*average(a'))
    me_a = meanfield(a, Δ*a'a, [a]; rates = [κ], efficiencies = [η])
    @test iszero(me_a.equations[1].rhs - test_eq_a)
    @test iszero(me_a.noise_equations[1].rhs - test_noise_eq_a)

    # Hamiltonian
    H =
        Δ*Σ(σ(2, 2, j), j) +
        Σ(δ(i)*b(i)'b(i), i) +
        gf*(Σ(a'*b(i) + a*b(i)', i)) +
        g*(Σ(a'*σ(1, 2, j) + a*σ(2, 1, j), j))

    # Jumps & rates
    J = [a, b(i), σ(1, 2, j), σ(2, 1, j), σ(2, 2, j)]
    rates = [κ, κf, Γ, R, ν]

    eqs = meanfield(a'a, H, J; rates = rates, order = 2)
    eqs_c = complete(eqs);

    efficiencies = [η, 0, 0, 0, 0]
    eqs_noise = meanfield(a'a, H, J; rates = rates, efficiencies = efficiencies, order = 2)
    eqs_c_noise = indexed_complete(eqs);

    @test length(eqs_c.equations) == length(eqs_c_noise.equations)
    for (eq, eq_noise) in zip(eqs_c.equations, eqs_c_noise.equations)
        @test isequal(eq.lhs, eq_noise.lhs)
        @test isequal(simplify(eq.rhs-eq_noise.rhs), 0)
    end

end
