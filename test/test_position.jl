using Qumulants
using Test

@testset "position" begin

h = PositionSpace(:X)

x = Position(h,:x)
p = Momentum(h,:p)

@test commutator(x,p) == im
@test commutator(p,x) == -im
@test simplify_operators(p*x + x*p) == p*x + x*p

# Test quantum Harmonic oscillator
@parameters m ω
H = 1/2*m*ω^2*x^2 + p^2/(2m)

@test iszero(commutator(H,H))

ops = [x,p]
he = heisenberg(ops,H)

@test he.rhs[1] == simplify_operators(p/m)
@test he.rhs[2] == simplify_operators(-m*ω*x)

r = symmetrize(x)
q = symmetrize(p)

@test commutator(H,r*q) == commutator(H,x*p)
@test simplify_operators(r*p) == simplify_operators(x*q)

ops2 = [x,p,x^2,p^2,r*q]
he2 = heisenberg(ops2,H)
he2_avg = average(he2,2)

ps = (m,ω)
meta_f = build_ode(he2_avg,ps)

end # testset
