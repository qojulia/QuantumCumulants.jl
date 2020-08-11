using Qumulants
using OrdinaryDiffEq
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
@test he.rhs[2] == simplify_operators(-m*ω^2*x)

r = symmetrize(x)
q = symmetrize(p)

@test commutator(H,r*q) == commutator(H,x*p)
@test simplify_operators(r*p) == simplify_operators(x*q)

ops2 = [x,p,x^2,p^2,r*q]
he2 = heisenberg(ops2,H)
he2_avg = average(he2,2)

ps = (m,ω)
f = generate_ode(he2_avg,ps)

# Initial state
α = 5.0im + 4 # coherent amplitude
u0 = [sqrt(2)*real(α), sqrt(2)*imag(α)]
p0 = (1,1,0.1)
u2 = zeros(ComplexF64, length(he2))
u2[1:2] = u0
u2[3] = real(α^2) + abs2(α) + 0.5
u2[4] = -real(α^2) + abs2(α) + 0.5
u2[5] = imag(α^2) + 0.5im

prob2 = ODEProblem(f,u2,(0.0,10.0),p0)
sol2 = solve(prob2,Tsit5(), abstol=1e-11, reltol=1e-9);

Δx = real.([u[3] - u[1]^2 for u in sol2.u])
Δp = real.([u[4] - u[2]^2 for u in sol2.u])
@test isapprox(minimum(Δx .* Δp), 0.25)


end # testset
