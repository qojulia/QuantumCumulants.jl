using Qumulants
using Test

@testset "spin" begin

h = SpinSpace(:S,1//2)
sx = SigmaX(h, :σ)
sy = SigmaY(h, :σ)
sz = SigmaZ(h, :σ)



end # testset
