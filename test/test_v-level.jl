using Qumulants
using OrdinaryDiffEq
using Test

@testset "v-level" begin

# Hilbertspace
hf = FockSpace(:cavity)
ha = NLevelSpace(:atom, 3)
h = hf⊗ha

# Parameters
@parameters κ g
Δc, Γ2, Γ3, Δ2, Δ3, Ω2, Ω3 = parameters("Δ_c Γ_2 Γ_3 Δ_2 Δ_3 Ω_2 Ω_3")

# Operators
a = Destroy(h,:a)
σ(i,j) = Transition(h,:σ,i,j)

# Hamiltonian
H_atom = -Δ2*σ(2,2) - Δ3*σ(3,3) + Ω2*(σ(2,1) + σ(1,2)) + Ω3*(σ(3,1) + σ(1,3))
H_cav = -Δc*a'*a + g*(a'*σ(1,2) + a*σ(2,1))
H = H_atom + H_cav

# Jump operators
J = [a,σ(1,2),σ(1,3)]
rates = [κ,Γ2,Γ3]
ops = [a'*a];

# Derive first equation
he_n = heisenberg([a'*a],H,J;rates=rates)

# Average and expand to 2nd order
hn_avg = average(he_n,2)

# Complete the system
he_comp = complete(hn_avg)

# Compare to finding operators from start
ops = find_operators(h,2)

@test length(ops)==length(he_comp.lhs)
he = heisenberg(ops,H,J;rates=rates)
he_avg = average(he,2)

@test all((l in he_comp.lhs || l' in he_comp.lhs) for l in he_avg.lhs)
@test all((l in he_comp.lhs || l' in he_avg.lhs) for l in he_comp.lhs)


p = [κ, g, Δc, Γ2, Γ3, Δ2, Δ3, Ω2, Ω3]
meta_f = build_ode(he_avg,p)
f = Meta.eval(meta_f)

u0 = zeros(ComplexF64,length(he_avg.lhs))

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

p0 = [κn, gn, Δcn, Γ2n, Γ3n, Δ2n, Δ3n, Ω2n, Ω3n]
prob = ODEProblem(f,u0,(0.0,tmax),p0)
sol = solve(prob,RK4());

# Filter cavity equations to compute spectrum
# Hilbert space
hfilter = FockSpace(:filter)
h_tot = h⊗hfilter

# New parameters
ωf, κf, gf = parameters("ω_f κ_f g_f")

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
he_f = heisenberg(ops_f,Hf,Jf;rates=rates_f)
he_f_avg = average(he_f,2)

# Find missing averages and them as parameter
missing_avgs = filter(x->isa(x,Average), find_missing(he_f_avg));

# Gather all new parameters
pf = [ωf; gf; κf; missing_avgs; p]

# Generate function for the filter cavities
meta_ff = build_ode(he_f_avg,pf);

# Filter cavity parameters
ωfn = 0.0
κfn = 0.05κn
gfn = 0.1κfn
tf = 5/Γ2n

# Numerical parameters - get steady state values
steady_vals = ComplexF64[]
avg_exprs = Qumulants._to_expression.(he_avg.lhs)
for m in missing_avgs
    m_ex = Qumulants._to_expression(m)
    m_adj_ex = Qumulants._to_expression(m')
    i = findfirst(isequal(m_ex), avg_exprs)
    j = findfirst(isequal(m_adj_ex), avg_exprs)
    if !(i isa Nothing)
        push!(steady_vals, sol.u[end][i])
    else
        push!(steady_vals, conj(sol.u[end][j]))
    end
end
pf0 = [ωfn;κfn;gfn;steady_vals;p0]

# Initial state
u0f = zeros(ComplexF64,length(he_f_avg.lhs))

ff = Meta.eval(meta_ff)
prob_f = ODEProblem(ff,u0f,(0.0,tf),pf0);

# Solve for different frequencies of the filters; the spectrum is then equal to ⟨fᵗf⟩
ω = range(-0.8,-0.2,length=81)
spec = zeros(length(ω))

freq_ind = findfirst(isequal(ωf),pf)
nf_idx = findfirst(isequal(c'*c),ops_f)
for i=1:length(ω)
    prob_f.p[freq_ind] = ω[i]
    sol_f = solve(prob_f,RK4(),save_idxs=nf_idx)
    spec[i] = real(sol_f.u[end])
end

@test all(spec .> 0.0)

end # testset
