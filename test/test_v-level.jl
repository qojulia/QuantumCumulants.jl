using QuantumCumulants
using OrdinaryDiffEq
using ModelingToolkit
using Test
using SymbolicUtils

@testset "v-level" begin

    # Hilbertspace
    hf = FockSpace(:cavity)
    ha = NLevelSpace(:atom, 3)
    h = hf⊗ha

    # Parameters
    @cnumbers g
    Δ2, Δ3, Ω2, Ω3 = cnumbers("Δ_2 Δ_3 Ω_2 Ω_3")
    Δc, Γ2 = rnumbers("Δ_c Γ_2")
    Γ3 = rnumber(:Γ3)
    Γ3 = rnumber("Γ3")
    @rnumbers κ

    @test isequal(κ', κ)
    @test isequal(conj(Γ3), Γ3)
    @test isequal(adjoint(κ), κ)

    @test typeof(Δ2*(Γ3 + 1)) == SymbolicUtils.BasicSymbolic{Complex{Real}}
    @test typeof(κ*(Γ3 + 1)) == SymbolicUtils.BasicSymbolic{Real}
    @test typeof(1im*κ*(Γ3 + 1)) == SymbolicUtils.BasicSymbolic{Complex{Real}}
    @test isequal((Δ2*(Γ3 + 1))', (conj(Δ2*(Γ3 + 1))))
    @test isequal(κ*(Γ3 + 1)', κ*(Γ3 + 1))
    @test isequal((1im*κ*(Γ3 + 1))', conj(1im*κ*(Γ3 + 1)))
    @test_broken isequal(simplify(exp(1im*κ)*(exp(1im*κ))'), 1)

    # Operators
    @qnumbers a::Destroy(h) σ::Transition(h)

    # Hamiltonian
    H_atom = -Δ2*σ(2, 2) - Δ3*σ(3, 3) + Ω2*(σ(2, 1) + σ(1, 2)) + Ω3*(σ(3, 1) + σ(1, 3))
    H_cav = -Δc*a'*a + g*(a'*σ(1, 2) + a*σ(2, 1))
    H = H_atom + H_cav
end # testset
