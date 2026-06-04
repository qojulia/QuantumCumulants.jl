using QCNew
using QuantumOpticsBase
using Symbolics: @variables
using Random
using Test

Random.seed!(0)

@testset "initial_values: Fock + 2-level state to u0 vector (integer levels)" begin
    h = FockSpace(:cavity) ⊗ NLevelSpace(:atom, 2)
    a = Destroy(h, :a)
    s(i, j) = Transition(h, :σ, i, j)

    @variables Δ::Real g::Real κ::Real η::Real
    H = Δ * a' * a + g * (a' * s(1, 2) + a * s(2, 1)) + η * (a + a')
    ops = [a, s(1, 2), a' * a, s(2, 2), a' * s(1, 2)]
    eqs = meanfield(ops, H, [a]; rates = [κ], order = 2)

    bcav = FockBasis(10)
    batom = NLevelBasis(2)
    b = bcav ⊗ batom
    ψ0 = randstate(b)

    u0 = initial_values(eqs, ψ0)

    @test u0[1] ≈ expect(destroy(bcav) ⊗ one(batom), ψ0)
    @test u0[2] ≈ expect(one(bcav) ⊗ QuantumOpticsBase.transition(batom, 1, 2), ψ0)
    @test u0[3] ≈ expect(number(bcav) ⊗ one(batom), ψ0)
    @test u0[4] ≈ expect(one(bcav) ⊗ QuantumOpticsBase.transition(batom, 2, 2), ψ0)
    @test u0[5] ≈ expect(create(bcav) ⊗ QuantumOpticsBase.transition(batom, 1, 2), ψ0)
end

@testset "initial_values: LazyKet route matches Ket route" begin
    if isdefined(QuantumOpticsBase, :LazyKet)
        h = FockSpace(:cavity) ⊗ NLevelSpace(:atom, 2)
        a = Destroy(h, :a)
        s(i, j) = Transition(h, :σ, i, j)

        @variables Δ::Real g::Real κ::Real
        H = Δ * a' * a + g * (a' * s(1, 2) + a * s(2, 1))
        ops = [a, s(1, 2), a' * a, s(2, 2)]
        eqs = meanfield(ops, H, [a]; rates = [κ], order = 2)

        bcav = FockBasis(8)
        batom = NLevelBasis(2)
        b = bcav ⊗ batom
        ψlazy = LazyKet(b, (randstate(bcav), randstate(batom)))
        ψfull = Ket(ψlazy)

        u0_lazy = initial_values(eqs, ψlazy)
        u0_full = initial_values(eqs, ψfull)
        @test u0_lazy ≈ u0_full
    else
        @info "QuantumOpticsBase.LazyKet not available; skipping"
    end
end
