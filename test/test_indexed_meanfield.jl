using QuantumCumulants
using QuantumOpticsBase
using ModelingToolkit
using OrdinaryDiffEq
using Test
using Random
using SteadyStateDiffEq

const qc = QuantumCumulants

@testset "meanfield" begin

order = 2
@cnumbers Δc η Δa κ

N = 2 #number of atoms
hc = FockSpace(:cavity)
ha = NLevelSpace(Symbol(:atom),2)
h = hc ⊗ ha

#define indices
i_ind = Index(h,:i,N,ha)
j_ind = Index(h,:j,N,ha)
k_ind = Index(h,:k,N,ha)

#define indexed variables
g(k) = IndexedVariable(:g,k)
Γ_ij = DoubleIndexedVariable(:Γ,i_ind,j_ind)
Ω_ij = DoubleIndexedVariable(:Ω,i_ind,j_ind;identical=false)

@qnumbers a::Destroy(h)
σ(i,j,k) = IndexedOperator(Transition(h,:σ,i,j),k)

# Hamiltonian

DSum = Σ(Ω_ij*σ(2,1,i_ind)*σ(1,2,j_ind),j_ind,i_ind;non_equal=true)

@test DSum isa DoubleSum
@test isequal(Σ(Σ(Ω_ij*σ(2,1,i_ind)*σ(1,2,j_ind),i_ind,[j_ind]),j_ind),DSum)

Hc = Δc*a'a + η*(a' + a)
Ha = Δa*Σ(σ(2,2,i_ind),i_ind) + DSum
Hi = Σ(g(i_ind)*(a'*σ(1,2,i_ind) + a*σ(2,1,i_ind)),i_ind)
H = Hc + Ha + Hi

J = [a, [σ(1,2,i_ind),σ(1,2,j_ind)] ] 
rates = [κ,Γ_ij]

ops = [a, σ(2,2,k_ind), σ(1,2,k_ind)]
eqs = meanfield(ops,H,J;rates=rates,order=order)

@test isequal([i_ind,j_ind,k_ind],sort(qc.get_indices_equations(eqs)))
@test isequal([:i,:j,:k],sort(qc.getIndName.(qc.get_indices_equations(eqs))))

@test length(eqs) == 3

ind1 = Index(h,:q,N,ha)
ind2 = Index(h,:r,N,ha)
ind3 = Index(h,:s,N,ha)

eqs_comp = qc.complete(eqs;extra_indices=[ind1,ind2,ind3])
eqs_comp2 = qc.complete(eqs)

@test length(eqs_comp.equations) == length(eqs_comp2.equations)

eqs_ = evaluate(eqs_comp)
eqs_2 = evaluate(eqs_comp2)

@test length(eqs_2) == length(eqs_)

@test length(eqs_) == 18

@named sys = ODESystem(eqs_);

u0 = zeros(ComplexF64, length(eqs_))
# parameter
Γ_ = 1.0
d = 2π*0.08 #0.08λ
θ = π/2

Ωij_(i,j) = Γ_*(-3/4)*( (1-(cos(θ))^2)*cos(d)/d-(1-3*(cos(θ))^2)*(sin(d)/(d^2)+(cos(d)/(d^3))) )
function Γij_(i,j)
    i==j ? Γ_ : Γ_*(3/2)*( (1-(cos(θ))^2)*sin(d)/d+(1-3*(cos(θ))^2)*((cos(d)/(d^2))-(sin(d)/(d^3))))
end

ΓMatrix = [Γij_(i,j) for i = 1:2, j=1:2]
ΩMatrix = [Ωij_(i,j) for i = 1:2, j=1:2]


g_ = 2Γ_
κ_ = 20Γ_
Δa_ = 0Γ_
Δc_ = 0Γ_
η_ = κ_/100

g_v = [g_*(-1)^j for j=1:2]
ps = [Δc, η, Δa, κ, g(i_ind), Γ_ij, Ω_ij];

eqs_4 = meanfield(ops,H,J;rates=rates,order=4)

Δc_i = -10*Γ_
Δa_i = Δc_i + Ωij_(1,2) #cavity on resonace with the shifted collective emitter
p0_i = [Δc_i, η_, Δa_i, κ_, g_v, ΓMatrix, ΩMatrix]

ps_ = value_map(ps,p0_i) #Combine all the parameters + values to one list for solving
prob = ODEProblem(sys,u0,(0.0, 20Γ_), ps_);
prob_ss = SteadyStateProblem(prob);

sol_ss = solve(prob_ss, DynamicSS(Tsit5(); abstol=1e-8, reltol=1e-8),
        reltol=1e-14, abstol=1e-14, maxiters=5e7)

@test length(eqs_4) == length(eqs)

order = 1

@cnumbers g N κ

# Hilbertspace
hc = FockSpace(:cavity)
hf = FockSpace(:filter)

h = hc ⊗ hf

i = Index(h,:i,N,hf)
j = Index(h,:j,N,hf)
k = Index(h,:k,N,hf)

xij = IndexedVariable(:x,i,j)


@qnumbers a_::Destroy(h,1)
b(k) = IndexedOperator(Destroy(h,:b,2), k)

H = g*a_'a_ + Σ(xij*a_'a_*b(i)'b(j),i,j)
J = [a_]
rates = [κ]


eqs1 = meanfield(a_,H,J;rates=rates,order=order) 
eqs2 = meanfield([a_],H,J;rates=rates,order=order) 
@test isequal(eqs1.equations,eqs2.equations)

@test isequal(sort([i,j]),sort(qc.get_all_indices(eqs1)))


#example for testing evaluation of individual hilbertspaces
@cnumbers N N2 Δ g κ Γ R ν M

hc = FockSpace(:cavity)
ha = NLevelSpace(:atom,2)

h = hc ⊗ ha

k = Index(h,:k,N,ha)
l = Index(h,:l,N,ha)

m = Index(h,:m,N2,hc)
n = Index(h,:n,N2,hc)

order = 2

σ(i,j,k) = IndexedOperator(Transition(h,:σ,i,j),k)
ai(k) = IndexedOperator(Destroy(h,:a),k)

H_2 = -Δ*∑(ai(m)'ai(m),m) + g*(∑(Σ(ai(m)'*σ(1,2,k),k),m) + ∑(Σ(ai(m)*σ(2,1,k),k),m))

J_2 = [ai(m),σ(1,2,k),σ(2,1,k),σ(2,2,k)]
rates_2 = [κ, Γ, R, ν]
ops_2 = [ai(n)'*ai(n),σ(2,2,l)]
eqs_2 = meanfield(ops_2,H_2,J_2;rates=rates_2,order=order)

q = Index(h,:q,N,ha)
r = Index(h,:r,N2,hc)

extra_indices = [q,r]

eqs_com = qc.complete(eqs_2;extra_indices=extra_indices);
eqs_com2 = qc.complete(eqs_2)

@test length(eqs_com) == 15
@test length(eqs_com) == length(eqs_com2)

@test isequal(sort([k,m,n,l,q,r]),sort(qc.get_all_indices(eqs_com)))

e_1 = evaluate(eqs_com; h=2,limits=(N=>5))
e_2 = evaluate(eqs_com; h=1,limits=(N2=>6))

@test length(e_1) != length(e_2)
@test !(e_1.equations == e_2.equations)
@test !(e_1.states == e_2.states)

@test sort(qc.get_indices_equations(e_1)) == sort([m,n,r])
@test sort(qc.get_indices_equations(e_2)) == sort([k,l,q])

limits = Dict(N=>5,N2=>6)
s1 = evaluate(eqs_com; h=[1,2],limits=limits)
s2 = evaluate(eqs_com; limits=limits)

s1_ = evaluate(eqs_com; h=[ha],limits=(N=>5))
s2_ = evaluate(eqs_com; h=[2],limits=(N=>5))

@test s1_.equations == s2_.equations

@test length(s1) == length(s2)
@test s1.equations == s2.equations

@test qc.get_indices_equations(s1) == []

@test s1.states == s2.states


end