using Qumulants
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
ops = [a'*a];

# Derive first equation
he_n = heisenberg([a'*a],H,J;rates=rates)

# Average and expand to 2nd order
hn_avg = average(he_n,2)

# Complete the system
he_comp = complete(hn_avg)

# Compare to finding operators from start
ops = find_operators(h,2)

@test length(ops)==length(he_comp.lhs)==16
he = heisenberg(ops,H,J;rates=rates)
he_avg = average(he,2)

# @test all((Qumulants._in(l, he_comp.lhs) || l' in he_comp.lhs) for l in he_avg.lhs)
@test all((Qumulants._in(l, he_comp.lhs) || Qumulants._in(Qumulants._adjoint(l), he_comp.lhs) for l in he_avg.lhs))
@test all((Qumulants._in(l, he_comp.lhs) || Qumulants._in(Qumulants._adjoint(l), he_avg.lhs)) for l in he_comp.lhs)


p = [κ, g, Δc, Γ2, Γ3, Δ2, Δ3, Ω2, Ω3]
sys = ODESystem(he_avg)

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

p0 = p .=> [κn, gn, Δcn, Γ2n, Γ3n, Δ2n, Δ3n, Ω2n, Ω3n]
prob = ODEProblem(sys,u0,(0.0,tmax),p0)
sol = solve(prob,RK4(),jac=true,sparse=true)

avg = average(a'*σ(2,1))
@test get_solution(avg,sol,he_avg) == get_solution(avg,sol.u,he_avg) == map(conj, getindex.(sol.u, 7))
@test get_solution(avg,sol,he_avg)[end] == get_solution(avg,sol.u[end],he_avg) == conj(sol.u[end][7])

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
he_f = heisenberg(ops_f,Hf,Jf;rates=rates_f)
he_f_avg = average(he_f,2)

# Find missing averages and them as parameter
import SymbolicUtils
missing_avgs = filter(SymbolicUtils.sym_isa(Average), find_missing(he_f_avg))
avg_ps = Qumulants._make_parameter.(missing_avgs)

he_f_avg = substitute(he_f_avg, Dict(missing_avgs .=> avg_ps))

# Gather all new cnumbers
pf = [ωf; gf; κf; avg_ps; p]

# Generate function for the filter cavities
sys = ODESystem(he_f_avg)

# Filter cavity cnumbers
ωfn = 0.0
κfn = 0.05κn
gfn = 0.1κfn
tf = 5/Γ2n

# Numerical cnumbers - get steady state values
steady_vals = ComplexF64[]
avg_exprs = Qumulants._to_expression.(he_avg.lhs)
for m in missing_avgs
    m_ex = Qumulants._to_expression(m)
    m_adj_ex = Qumulants._to_expression(Qumulants._adjoint(m))
    i = findfirst(isequal(m_ex), avg_exprs)
    j = findfirst(isequal(m_adj_ex), avg_exprs)
    if !(i isa Nothing)
        push!(steady_vals, sol.u[end][i])
    else
        push!(steady_vals, conj(sol.u[end][j]))
    end
end
pf0 = pf .=> [ωfn;κfn;gfn;steady_vals;getindex.(p0, 2)]

# Initial state
u0f = zeros(ComplexF64,length(he_f_avg.lhs))

prob_f = ODEProblem(sys,u0f,(0.0,tf),pf0,jac=true,sparse=false)

# Solve for different frequencies of the filters; the spectrum is then equal to ⟨fᵗf⟩
ω = range(-0.8,-0.2,length=81)
spec = zeros(length(ω))

freq_ind = findfirst(isequal(ωf),pf)
n_avg = average(c'*c)
for i=1:length(ω)
    # TODO fix indexing here
    prob_f.p[freq_ind] = ω[i]
    sol_f = solve(prob_f,RK4())
    spec[i] = real(get_solution(n_avg, sol_f.u[end], he_f_avg))
end

@test all(spec .> 0.0)

end # testset
