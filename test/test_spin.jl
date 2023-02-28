using QuantumCumulants
using Test

@testset "spin" begin

hs1 = SpinSpace(:Spin1)
hs2 = SpinSpace(:Spin2)
h = hs1 ⊗ hs2

s(axis) = Sigma(hs1, :σ, axis) # axis ∈ [1,2,3] → [x,y,z]
@test isequal(s(1)*s(2),1im*s(3))
@test !isequal(s(1)*s(2),1im*s(2))
@test isequal(s(1)*s(3),-1im*s(2))
@test isequal(s(3)*s(3),1)
@test isequal(s(3)*s(1),1im*s(2))
@test isequal(s(1)*s(2)*s(3),1im)

σ(i, axis) = Sigma(h,Symbol(:σ_,i), axis, i)
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
Δ(i) = cnumber(Symbol(:Δ_,i))
H = Δ(1)*σz(1) + Δ(2)*σz(2) + J*σx(1)*σx(2)

ops = [σz(1)]
eqs = meanfield(ops, H)
eqs_c = complete(eqs, order=2)
@test length(eqs_c) == 6

ops2 = [σx(1), σy(1), σz(1), σx(2), σy(2), σz(2)]
eqs2 = meanfield(ops2, H)
eqs2_c = complete(eqs2, order=2)
@test length(eqs2_c) == 14

# TODO - Add more test: timeevolution, comparison NLevelSpace, ...
end # testset
