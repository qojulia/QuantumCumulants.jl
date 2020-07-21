using Qumulants
using Test

@testset "position" begin

h = PositionSpace(:X)

x = Position(h,:x)
p = Momentum(h,:p)

@test commutator(x,p) == im
@test commutator(p,x) == -im

@test commutator(cos(x),p) == sin(x)
@test commutator(sin(x),p) == -1*cos(x)

@parameters k q r
@test commutator(cos(k*x),p) == k*sin(k*x)
@test commutator(sin(k*x),p) == simplify_operators(-k*cos(k*x))
@test commutator(cos(q*k*x),r*p) == simplify_operators(r*q*k*sin(q*k*x))
@test commutator(sin(q*k*x),r*p) == simplify_operators(-r*q*k*cos(q*k*x))
@test iszero(commutator(cos(k*x),q*x))

# Test two-level atom coupled to standing wave
hx = PositionSpace(:motion)
ha = NLevelSpace(:internal, (:g,:e))
h = ha⊗hx

x = Position(h,:x)
p = Momentum(h,:p)
σ(i,j) = Transition(h,:σ,i,j)

@parameters Δ Ω k γ m
H = -Δ*σ(:e,:e) + Ω*cos(k*x)*(σ(:e,:g) + σ(:g,:e)) + p^2/(2m)
J = [σ(:g,:e)]

dx = heisenberg(x,H,J;rates=[γ])
@test dx.rhs[1] == simplify_operators(p/m)

dp = heisenberg(p,H,J;rates=[γ])
@test dp.rhs[1] == simplify_operators(im*k*Ω*sin(k*x)*(σ(:e,:g) + σ(:g,:e)))

ops = [σ(:g,:e),σ(:e,:e),x,p]
he = heisenberg(ops,H,J;rates=[γ])

end # testset
