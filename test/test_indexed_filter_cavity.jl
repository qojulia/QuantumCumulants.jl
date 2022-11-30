using QuantumCumulants
using OrdinaryDiffEq, SteadyStateDiffEq, ModelingToolkit
using Test

@testset "filter-cavity-indexing" begin
### indexing results
order = 2 #order of the cumulant expansion
@cnumbers κ g gf κf R Γ Δ ν N M
δ(i) = IndexedVariable(:δ, i)

# Hilbertspace
hc = FockSpace(:cavity)
hf = FockSpace(:filter)
ha = NLevelSpace(:atom, 2)
h = hc ⊗ hf ⊗ ha

# Indices and Operators
i = Index(h,:i,M,hf)
i1 = Index(h,:i1,M,hf)
i2 = Index(h,:i2,M,hf)

j = Index(h,:j,N,ha)
j1 = Index(h,:j1,N,ha)
j2 = Index(h,:j2,N,ha)

@qnumbers a::Destroy(h,1)
b(k) = IndexedOperator(Destroy(h,:b,2), k)
σ(α,β,k) = IndexedOperator(Transition(h,:σ,α,β,3), k)

# Hamiltonian
H = Δ*Σ(σ(2,2,j),j) + Σ(δ(i)*b(i)'b(i),i) +
    gf*(Σ(a'*b(i) + a*b(i)',i)) + g*(Σ(a'*σ(1,2,j) + a*σ(2,1,j),j))

# Jumps & rates
J = [a, b(i), σ(1,2,j), σ(2,1,j), σ(2,2,j)]
rates = [κ, κf, Γ, R, ν]

ops = [a'a]
eqs = meanfield(ops,H,J;rates=rates,order=order)
eqs_c = complete(eqs)
@test length(eqs_c) == 23

# first evaluate, then scale
M_ = 3 #number of filter cavities
eqs_eval_ = evaluate(eqs_c; limits=(M=>M_), h=[2])
@test length(eqs_eval_) == 43
eqs_sc_ = scale(eqs_eval_;h=[3])
eqs_sc_test = scale(eqs_eval_;h=[ha])
@test all([eqs_sc_[i] == eqs_sc_test[i] for i=1:length(eqs_sc_)])
@test length(eqs_sc_) == 42

# first scale, then evaluate
eqs_sc = scale(eqs_c;h=[3])
@test length(eqs_sc) == 22
eqs_eval = evaluate(eqs_sc; limits=(M=>M_))
@test length(eqs_eval) == 42
@test length(eqs_sc_) == length(eqs_eval)

# Initial state
u0 = zeros(ComplexF64, length(eqs_eval))

# System parameters
N_ = 1e4
Γ_ = 1.0
Δ_ = 0Γ_
g_ = 2Γ_
κ_ = 1e3Γ_
R_ = 1e2Γ_
ν_ = 10Γ_

gf_ = 0.1Γ_
κf_ = 0.1Γ_
δ_ls = [0:1/M_:1-1/M_;]*10Γ_

ps = [Γ, κ, g, κf, gf, R, [δ(i) for i=1:M_]..., Δ, ν, N]
p0 = [Γ_, κ_, g_, κf_, gf_, R_, δ_ls..., Δ_, ν_, N_]

@named sys = ODESystem(eqs_eval)
prob = ODEProblem(sys,u0,(0.0, 1.0/κf_), ps.=>p0)
sol = solve(prob, Tsit5(); maxiters=1e7, abstol=1e-12, reltol=1e-12)
n_ind = sol[a'a][end]
nf_ind = [sol[b(k)'b(k)][end] for k=1:M_]
s22_ind = sol[σ(2,2,1)][end]

@named sys2 = ODESystem(eqs_sc_)
prob2 = ODEProblem(sys2,u0,(0.0, 1.0/κf_), ps.=>p0)
sol2 = solve(prob2, Tsit5(); maxiters=1e7, abstol=1e-12, reltol=1e-12)
n_ind2 = sol2[a'a][end]
nf_ind2 = [sol2[b(k)'b(k)][end] for k=1:M_]
s22_ind2 = sol2[σ(2,2,1)][end]

@test n_ind == n_ind2
@test nf_ind == nf_ind2
@test s22_ind == s22_ind2

### brute-force for M filter cavities
hc_ = FockSpace(:cavity)
hf_ = [FockSpace(Symbol(:filter,i)) for i=1:M_]
ha__ = NLevelSpace(:atom,2)
ha_ = ClusterSpace(ha__, N, order)
h_ = tensor(hc_, hf_..., ha_)

@qnumbers c::Destroy(h_,1)
f(i) = Destroy(h_,Symbol(:f,i),1+i)
s(i,j) = Transition(h_, :s, i, j, 2+M_)
s(1,2)

d(i) = cnumbers(Symbol(:d,i))[1]

H_ = Δ*sum(s(2,2)) + sum(d(i)*f(i)'f(i) for i=1:M_) +
    gf*(sum(c'*f(i) + c*f(i)' for i=1:M_)) + g*(c'*sum(s(1,2)) + c*sum(s(2,1)))

# Jumps & rates
J_ = [c, [f(i) for i=1:M_]..., s(1,2), s(2,1), s(2,2)]
rates_ = [κ, [κf for i=1:M_]..., Γ, R, ν]

ops_ = [c'c]
eqs_ = meanfield(ops_,H_,J_;rates=rates_,order=order)
eqs_c_ = complete(eqs_)
@test length(eqs_c_) == length(eqs_eval) == length(eqs_sc_)

@named sys_ = ODESystem(eqs_c_)

ps_ = [Γ, κ, g, κf, gf, R, [d(i) for i=1:M_]..., Δ, ν, N]
p0_ = [Γ_, κ_, g_, κf_, gf_, R_, δ_ls..., Δ_, ν_, N_]
u0_ = zeros(ComplexF64, length(eqs_c_))

prob_ = ODEProblem(sys_,u0_,(0.0, 1.0/κf_), ps_.=>p0_)
sol_ = solve(prob_, Tsit5(); maxiters=1e7, abstol=1e-12, reltol=1e-12)

n_bf = sol_[c'c][end]
nf_bf = [sol_[f(k)'f(k)][end] for k=1:M_]
s22_bf = sol_[s(2,2)[1]][end]

@test n_ind ≈ n_bf
@test nf_ind ≈ nf_bf
@test s22_ind ≈ s22_bf

end #testset
