using QuantumCumulants
using SymbolicUtils: SymbolicUtils
using Test

@testset "find_missing on JC" begin
    hc = FockSpace(:cavity); ha = NLevelSpace(:atom, 2)
    h = hc ⊗ ha
    @qnumbers a::Destroy(h, 1)
    σ(i, j) = Transition(h, :σ, i, j, 2)
    @variables Δ g κ γ
    H = Δ * a' * a + g * (a * σ(2, 1) + a' * σ(1, 2))
    eqs = meanfield([a, σ(2, 2)], H, [a, σ(1, 2)]; rates = [κ, γ], order = 2)
    missing_states = find_missing(eqs)
    @test !isempty(missing_states)
    @test all(
        SymbolicUtils.iscall(m) && QuantumCumulants.SecondQuantizedAlgebra.is_average(m)
            for m in missing_states
    )
end

@testset "complete! closes JC at order 2" begin
    hc = FockSpace(:cavity); ha = NLevelSpace(:atom, 2)
    h = hc ⊗ ha
    @qnumbers a::Destroy(h, 1)
    σ(i, j) = Transition(h, :σ, i, j, 2)
    @variables Δ g κ γ
    H = Δ * a' * a + g * (a * σ(2, 1) + a' * σ(1, 2))
    eqs = meanfield([a, σ(2, 2)], H, [a, σ(1, 2)]; rates = [κ, γ], order = 2)
    complete!(eqs)
    @test isempty(find_missing(eqs))
    @test length(eqs.equations) >= 2
end

@testset "complete is non-mutating" begin
    hc = FockSpace(:cavity)
    @qnumbers a::Destroy(hc)
    @variables ω κ U
    H = ω * a' * a + U * a' * a' * a * a
    eqs = meanfield([a], H, [a]; rates = [κ], order = 2)
    n_before = length(eqs.equations)
    eqs2 = complete(eqs)
    @test length(eqs.equations) == n_before
    @test isempty(find_missing(eqs2))
end
