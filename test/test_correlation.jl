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
@cnumbers Δ g κ γ ν

H = Δ*a'*a + g*(a'*σ + σ'*a)
J = [a,σ,σ']
he_laser = heisenberg([a'*a,σ'*σ,a*σ'],H,J;rates=[κ,γ,ν],multithread=true,simplify_input=true)

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
sol = solve(prob,RK4())
n = getindex.(sol.u,1)
pe = getindex.(sol.u,2)

@test all(iszero.(imag.(n)))
@test all(iszero.(imag.(pe)))
@test all(real.(n) .>= 0.0)
@test all(1.0 .>= real.(pe) .>= 0.0)

# Correlation function
c_steady = CorrelationFunction(a', a, he_comp; steady_state=true, multithread=true)
cf = generate_ode(c_steady, ps; check_bounds=true)

u0_c = initial_values(c_steady, sol.u[end])
p0_c = (p0..., sol.u[end]...)
τ = range(0.0, 10tmax; length=1001)
prob_c = ODEProblem(cf,u0_c,(0.0,τ[end]),p0_c)
sol_c = solve(prob_c,RK4(),saveat=τ,save_idxs=1)

# Spectrum via FFT of g(τ)
using QuantumOptics.timecorrelations: correlation2spectrum
ω, S1 = correlation2spectrum(sol_c.t, sol_c.u)
# using PyPlot; pygui(true)
# plot(ω, S1 ./ maximum(S1), label="FFT of g(τ)")

# Spectrum via Laplace transform
S = Spectrum(c_steady, ps)
usteady = sol.u[end]

S2 = S(ω,usteady,p0)
# plot(ω, S2 ./ maximum(S2), label="Laplace transform")
# xlim(-2pi,2pi)

# S_ = S(ω_ls,usteady,p0)
# max, ind = findmax(S_)
# hm_idx = findmin(abs.(S_ .- 0.5max))[2]
# fwhm = 2*abs(ω_ls[ind] - ω_ls[hm_idx])

S1_ = S1 ./ maximum(S1)
S1_ .-= minimum(S1_)
S_check = abs.(S2 ./ maximum(S2) .- S1_)
@test maximum(S_check) < 0.02

# Phase invariant case
has_phase(x) = !iszero(phase(x))
import SymbolicUtils
phase(avg::Average) = phase(avg.arguments[1])
phase(op::Destroy) = -1
phase(op::Create) = 1
function phase(t::Transition)
    lvls = Qumulants.levels(t.hilbert, t.aon)
    i = findfirst(isequal(t.i), lvls)
    j = findfirst(isequal(t.j), lvls)
    if i < j
        -1
    elseif i > j
        1
    else
        0
    end
end
phase(op::QTerm{<:typeof(*)}) = sum(phase(arg) for arg in op.arguments)
c_nophase = CorrelationFunction(a', a, he_avg; steady_state=true, filter_func=!has_phase)

S_nophase = Spectrum(c_nophase, ps)
S3 = S_nophase(ω,usteady,p0)
# plot(ω, S3 ./ maximum(S3), label="Laplace transform (phase invariant)")
# legend()

# When not in steady state -- cavity that decays
h = FockSpace(:fock)
a = Destroy(h,:a)
@cnumbers ωc κ
H = ωc*a'*a
he = heisenberg(a'*a,H,[a];rates=[κ])
he_avg = average(he)
ps = (ωc,κ)
f = generate_ode(he_avg,ps)
n0 = 20.0
u0 = [n0]
p0 = (1,1)
prob = ODEProblem(f,u0,(0.0,10.0),p0)
sol = solve(prob,RK4())

c = CorrelationFunction(a', a, he_avg)
cf = generate_ode(c, ps)
idx = 5
u0_c = initial_values(c, sol.u[idx])
prob_c = ODEProblem(cf,u0_c,(0.0,10.0),p0)
sol_c = solve(prob_c,RK4(),save_idxs=1)
# plot(sol_c.t,real.(sol_c.u), label="Re(g) -- numeric")
# plot(sol_c.t,imag.(sol_c.u), label="Im(g) -- numeric")

gfunc(τ) = @. sol.u[idx] * exp((im*p0[1]-0.5p0[2])*τ)
# plot(sol_c.t, real.(gfunc(sol_c.t)), ls="dashed", label="Re(g) -- analytic")
# plot(sol_c.t, imag.(gfunc(sol_c.t)), ls="dashed", label="Re(g) -- analyitc")
# legend()

@test isapprox(sol_c.u, gfunc(sol_c.t), rtol=1e-4)

end # testset
