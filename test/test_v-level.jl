using QuantumCumulants
using OrdinaryDiffEq
using ModelingToolkit
using Test

@testset "v-level" begin

# Hilbertspace
hf = FockSpace(:cavity)
ha = NLevelSpace(:atom, 3)
h = hf⊗ha

# Parameters
@cnumbers κ g
Δc, Γ2, Γ3, Δ2, Δ3, Ω2, Ω3 = cnumbers("Δ_c Γ_2 Γ_3 Δ_2 Δ_3 Ω_2 Ω_3")

# Operators
@qnumbers a::Destroy(h) σ::Transition(h)

# Hamiltonian
H_atom = -Δ2*σ(2,2) - Δ3*σ(3,3) + Ω2*(σ(2,1) + σ(1,2)) + Ω3*(σ(3,1) + σ(1,3))
H_cav = -Δc*a'*a + g*(a'*σ(1,2) + a*σ(2,1))
H = H_atom + H_cav

# Jump operators
J = [a,σ(1,2),σ(1,3)]
rates = [κ,Γ2,Γ3]
ops = [a'*a]

# Derive first equation
he_n = meanfield([a'*a],H,J;rates=rates)

# Average and expand to 2nd order
hn_avg = cumulant_expansion(he_n,2)

# Complete the system
he_comp = complete(hn_avg)

# Compare to finding operators from start
ops = find_operators(h,2)

@test length(ops)==length(he_comp)==16
he = meanfield(ops,H,J;rates=rates)
he_avg = cumulant_expansion(he,2)

@test all(((l ∈ Set(he_comp.states) || QuantumCumulants._adjoint(l) ∈ Set(he_comp.states))) for l in he_avg.states)
@test all(((l ∈ Set(he_comp.states) || QuantumCumulants._adjoint(l) ∈ Set(he_avg.states))) for l in he_comp.states)


p = [κ, g, Δc, Γ2, Γ3, Δ2, Δ3, Ω2, Ω3]
sys = ODESystem(he_avg)

u0 = zeros(ComplexF64,length(he_avg))

# Parameters - numerical values
Γ3n = 1.0
Γ2n = Γ3n/4266
κn = Γ3n/32.0
Δ2n = 0.5 * Γ3n
Δ3n = -0.2 .* Γ3n
Δcn = Δ2n
Ω2n = 0.2 * Γ3n
Ω3n = 0.05 * Γ3n
gn = 10.0 * Γ2n
tmax = 5/Γ2n

p0 = p .=> [κn, gn, Δcn, Γ2n, Γ3n, Δ2n, Δ3n, Ω2n, Ω3n]
prob = ODEProblem(sys,u0,(0.0,tmax),p0)
sol = solve(prob,RK4(),jac=true,sparse=true)

avg = average(a'*σ(2,1))
@test sol[avg] == map(conj, getindex.(sol.u, 7)) == map(conj, sol[QuantumCumulants._conj(avg)])

# Filter cavity equations to compute spectrum
# Hilbert space
hfilter = FockSpace(:filter)
h_tot = h⊗hfilter

# New cnumbers
ωf, κf, gf = cnumbers("ω_f κ_f g_f")

# Operators
c = Destroy(h_tot,:c, 3)
af = Destroy(h_tot,:a, 1)
σf(i,j) = Transition(h_tot, :σ, i, j)

# Hamiltonian and jumps
H_atom_f = -Δ2*σf(2,2) - Δ3*σf(3,3) + Ω2*(σf(2,1) + σf(1,2)) + Ω3*(σf(3,1) + σf(1,3))
H_cav_f = -Δc*af'*af + g*(af'*σf(1,2) + af*σf(2,1))
Hf = ωf*c'*c + gf*(af'*c + af*c') + H_atom_f + H_cav_f

Jf = [af,σf(1,2),σf(1,3),c]
rates_f = [rates;κf]

# Gather operators
ops_f = find_operators(h_tot,2;names=[:a,:σ,:c])
filter!(x->(3 in acts_on(x)), ops_f)

# Compute equations
he_f = meanfield(ops_f,Hf,Jf;rates=rates_f)
he_f_avg = cumulant_expansion(he_f,2)

# Find missing averages and them as parameter
missing_avgs = find_missing(he_f_avg)
avg_ps = QuantumCumulants._make_parameter.(missing_avgs)

he_f_avg = substitute(he_f_avg, Dict(missing_avgs .=> avg_ps))

# Gather all new cnumbers
pf = [ωf; gf; κf; avg_ps; p]

# Generate function for the filter cavities
sys_f = ODESystem(he_f_avg)

end # testset
