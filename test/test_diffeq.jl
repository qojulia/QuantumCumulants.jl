using Qumulants
using OrdinaryDiffEq
using Test

@testset "diffeq" begin

hf = FockSpace(:cavity)
ha = NLevelSpace(:atom,(:g,:e))
h = hf⊗ha

a = Destroy(h,:a)
σ = Transition(h,:σ,:g,:e)

# Single-atom laser
@cnumbers Δ g κ γ ν

H = Δ*a'*a + g*(a'*σ + σ'*a)
J = [a,σ,σ']
he_laser = heisenberg([a'*a,σ'*σ,a*σ'],H,J;rates=[κ,γ,ν],multithread=true,simplify_input=true)

he_avg = average(he_laser;multithread=true)
he_exp = cumulant_expansion(he_avg,2;multithread=true)
@test isequal(he_exp, average(he_laser,2))

ps = [Δ,g,κ,γ,ν]
missed = find_missing(he_exp)
@test !any(Qumulants._in(p, missed) for p=ps)

# Exploit phase invariance
subs = Dict([missed; Qumulants._conj.(missed)] .=> 0)
he_nophase = substitute(he_exp, subs)
@test isempty(find_missing(he_nophase))

f = generate_ode(he_nophase,ps)

# Numerical solution
p0 = [0.0,0.5,1.0,0.1,0.9]
u0 = zeros(ComplexF64,3)
tmax = 10.0

prob = ODEProblem(f,u0,(0.0,tmax),p0)
sol = solve(prob,RK4());
n = getindex.(sol.u,1)
pe = getindex.(sol.u,2)

@test all(iszero.(imag.(n)))
@test all(iszero.(imag.(pe)))
@test all(real.(n) .>= 0.0)
@test all(1.0 .>= real.(pe) .>= 0.0)

end # testset
