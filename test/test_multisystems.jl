using Qumulants
using Test
import OrdinaryDiffEq
import SymPy
ODE = OrdinaryDiffEq

@testset "multisystems" begin

# Parameters
N = 2
g = 1.5
γ = 0.25
κ = 1.0
ν = 4.0
ωc = ωa = 0.0

a = embed(Destroy(:a), 1, N+1)
s = [embed(Transition(Symbol("σ$i"),:g,:e,(:g,:e)), i+1, N+1) for i=1:N]
sps = [simplify_operators(s1'*s1) for s1=s]
# TODO: off-diagonal products for atoms

H = ωc*a'*a + ωa*sum(sps) + g*sum(a'*s1 + a*s1' for s1=s)
J = [sqrt(κ)*a]
for s1=s
    push!(J,sqrt(γ)*s1)
end
for s1=s
    push!(J,sqrt(ν)*s1')
end
ops = [a'*a; sps; [a'*s1 for s1=s]]

eqs_op = heisenberg(ops,H,J)
eqs_avg = average(eqs_op,2)

meta_f = build_ode(eqs_avg;set_unknowns_zero=true)
f = generate_ode(eqs_avg;set_unknowns_zero=true)

u0 = zeros(ComplexF64,length(ops))
prob = ODE.ODEProblem{true}(f,u0,(0.0,10.0))

sol = ODE.solve(prob,ODE.Tsit5());

n = getindex.(sol.u,1)
pop1 = getindex.(sol.u,2)

@test imag.(n) == zeros(length(n))
@test imag.(pop1) == zeros(length(n))
@test all(0.0 .<= real.(n))
@test all(0.0 .<= real.(pop1) .<= 1.0)

end # testset
