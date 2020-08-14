using Qumulants
using Test

@testset "spin" begin

h = SpinSpace(:S,1)
sx = SigmaX(h, :σ)
sy = SigmaY(h, :σ)
sz = SigmaZ(h, :σ)

@test simplify_operators(sy*sx) == sx*sy + -2im*sz
@test simplify_operators(sz*sy) == -2im*sx + sy*sz
@test simplify_operators(sz*sx) == sx*sz + 2im*sy

@test simplify_operators(sx^2) != 1
@test simplify_operators(sx^5) != sx
@parameters p
@test simplify_operators(sx*p*sx) == p*sx^2


h = SpinSpace(:S,1//2)
sx = SigmaX(h, :σ)
sy = SigmaY(h, :σ)
sz = SigmaZ(h, :σ)

@test simplify_operators(sy*sx) == -im*sz
@test simplify_operators(sz*sy) == -im*sx
@test simplify_operators(sz*sx) == im*sy

@test simplify_operators(sx^2) == 1
@test simplify_operators(sx^5) == sx
@parameters p
@test simplify_operators(sx*p*sx) == p


end # testset
