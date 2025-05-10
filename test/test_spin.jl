using QuantumCumulants
using ModelingToolkit
using OrdinaryDiffEq
using QuantumOpticsBase
using Test

@testset "spin" begin
    hs1 = PauliSpace(:Spin1)
    hs2 = PauliSpace(:Spin2)
    h = hs1 ⊗ hs2

    s(axis) = Pauli(hs1, :σ, axis) # axis ∈ [1,2,3] → [x,y,z]
    @test isequal(s(1)*s(2), 1im*s(3))
    @test !isequal(s(1)*s(2), 1im*s(2))
    @test isequal(s(1)*s(3), -1im*s(2))
    @test isequal(s(3)*s(3), 1)
    @test isequal(s(3)*s(1), 1im*s(2))
    @test isequal(s(1)*s(2)*s(3), 1im)

    σ(i, axis) = Pauli(h, Symbol(:σ_, i), axis, i)
    σx(i) = σ(i, 1)
    σy(i) = σ(i, 2)
    σz(i) = σ(i, 3)

    sx(i) = σ(i, :x)
    sy(i) = σ(i, :y)
    sz(i) = σ(i, :z)
    sx(1) == σx(1)

    @test isequal(σx(2)*σx(1), σx(1)*σx(2))
    @test isequal(σy(2)*σx(1), σx(1)*σy(2))
    @test isequal(σz(2)*σz(1), σz(1)*σz(2))

    @cnumbers J
    Δi(i) = cnumber(Symbol(:Δ_, i))
    H = Δi(1)*σz(1) + Δi(2)*σz(2) + J*σx(1)*σx(2)

    ops = [σz(1)]

end # testset
