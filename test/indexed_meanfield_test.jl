using QuantumCumulants
using Symbolics: Symbolics, @variables
using Test

@testset "indexed meanfield: collective atomic emission" begin
    hc = FockSpace(:cavity); ha = NLevelSpace(:atom, 2); h = hc ⊗ ha
    @qnumbers a::Destroy(h, 1)
    σ(x, y, k) = IndexedOperator(Transition(h, :σ, x, y, 2), k)
    @variables N Δ g κ γ
    i = Index(h, :i, N, ha)
    H = Δ * a' * a + Σ(g * (a * σ(2, 1, i) + a' * σ(1, 2, i)), i)
    eqs = meanfield([a], H, [a, σ(1, 2, i)]; rates = [κ, γ])
    @test eqs isa MeanFieldEquations
    @test length(eqs.equations) == 1
end

@testset "indexed average appears in RHS" begin
    hc = FockSpace(:cavity); ha = NLevelSpace(:atom, 2); h = hc ⊗ ha
    @qnumbers a::Destroy(h, 1)
    σ(x, y, k) = IndexedOperator(Transition(h, :σ, x, y, 2), k)
    @variables N Δ g κ γ
    i = Index(h, :i, N, ha)
    H = Δ * a' * a + Σ(g * (a * σ(2, 1, i) + a' * σ(1, 2, i)), i)
    eqs = meanfield([a], H, [a, σ(1, 2, i)]; rates = [κ, γ])
    @test occursin("σ", repr(eqs.equations[1].rhs))
end

@testset "indexed meanfield: two-atom collective JC drift" begin
    @variables Δc::Real η::Real Δa::Real κ::Real Γ::Real
    N = 2
    hc = FockSpace(:cavity)
    ha = NLevelSpace(Symbol(:atom), 2)
    h = hc ⊗ ha

    i_ind = Index(h, :i, N, ha)
    k_ind = Index(h, :k, N, ha)

    g(k) = IndexedVariable(:g, k)

    @qnumbers a::Destroy(h)
    σ(i, j, k) = IndexedOperator(Transition(h, :σ, i, j), k)

    Hc = Δc * a' * a + η * (a' + a)
    Ha = Δa * Σ(σ(2, 2, i_ind), i_ind)
    Hi = Σ(g(i_ind) * (a' * σ(1, 2, i_ind) + a * σ(2, 1, i_ind)), i_ind)
    H = Hc + Ha + Hi

    J = [a, σ(1, 2, i_ind)]
    rates = [κ, Γ]

    ops = [a, σ(2, 2, k_ind), σ(1, 2, k_ind)]
    eqs = meanfield(ops, H, J; rates = rates, order = 2)
    @test length(eqs.equations) == 3
    @test eqs isa MeanFieldEquations
    ceqs = complete(eqs)
    @test ceqs isa MeanFieldEquations
    @test length(ceqs.equations) >= 3
end
