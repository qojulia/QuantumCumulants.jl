using Qumulants
using OrdinaryDiffEq
using Test

@testset "correlation" begin

hf = FockSpace(:cavity)
ha = NLevelSpace(:atom,(:g,:e))
h = hf⊗ha

a = Destroy(h,:a)
σ = Transition(h,:σ,:g,:e)

# Single-atom laser
@parameters Δ g κ γ ν

H = Δ*a'*a + g*(a'*σ + σ'*a)
J = [a,σ,σ']
he_laser = heisenberg([a'*a,σ'*σ,a*σ'],H,J;rates=[κ,γ,ν],multithread=true)

he_avg = average(he_laser,2;multithread=true)
he_comp = complete(he_avg)

ps = (Δ, g, γ, κ, ν)
f = generate_ode(he_comp,ps)

# Numerical solution
# p0 = [0.0,0.5,1.0,0.1,0.9]
p0 = (1, 1.5, 0.25, 1, 4)
u0 = zeros(ComplexF64,length(he_comp))
tmax = 10.0

prob = ODEProblem(f,u0,(0.0,tmax),p0)
sol = solve(prob,RK4());
n = getindex.(sol.u,1)
pe = getindex.(sol.u,2)

@test all(iszero.(imag.(n)))
@test all(iszero.(imag.(pe)))
@test all(real.(n) .>= 0.0)
@test all(1.0 .>= real.(pe) .>= 0.0)

# Correlation function
c = CorrelationFunction(a', a, he_comp)


c_steady = CorrelationFunction(a', a, he_comp; steady_state=true)

# Spectrum
S = Spectrum(c_steady, ps)

usteady = sol.u[end]
ω_ls = range(-pi, pi, length=501)

using PyPlot; pygui(true)
plot(ω_ls, S(ω_ls,usteady,p0))

S_ = S(ω_ls,usteady,p0)
max, ind = findmax(S_)
hm_idx = findmin(abs.(S_ .- 0.5max))[2]
fwhm = 2*abs(ω_ls[ind] - ω_ls[hm_idx])

end # testset
