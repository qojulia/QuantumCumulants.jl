using Qumulants
using Test
import OrdinaryDiffEq
ODE = OrdinaryDiffEq

@testset "diffeq" begin

# Test single-atom laser
σ(i,j) = Transition(:σ,i,j,2;GS=1)
a = Destroy(:a) ⊗ Identity()
s = Identity() ⊗ σ(1,2)
sz = simplify_operators(2*s'*s - one(s))
sps = simplify_operators(s'*s)
ωc, ωa, g, γ, κ, ν = (0.1, 3, 0.5, 0.25, 1.0, 4)
H = ωc*a'*a + ωa*s'*s + g*(a'*s + a*s')
J = [sqrt(κ)*a,sqrt(γ)*s,sqrt(ν)*s']

ops = [a'*a,sps,a'*s]
eqs_op = heisenberg(ops,H,J)
eqs_avg = average(eqs_op,2)

meta_f = build_ode(eqs_avg,set_unknowns_zero=true)
f = generate_ode(eqs_avg,set_unknowns_zero=true)

u0 = zeros(ComplexF64,3)
prob = ODE.ODEProblem{true}(f,u0,(0.0,10.0))
sol = ODE.solve(prob,ODE.Tsit5())

n = getindex.(sol.u,1)
p = getindex.(sol.u,2)
@test imag.(n) == zeros(length(n))
@test all(real.(n) .>= 0.0)
@test imag.(p) == zeros(length(p))
@test all(0.0 .<= real.(p) .<= 1.0)

end # testset
