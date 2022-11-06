using Test
using QuantumCumulants
using SymbolicUtils
using Symbolics

const qc = QuantumCumulants

@testset "indexed_meanfield" begin

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

@test DSum isa IndexedDoubleSum
@test isequal(Σ(Σ(Ω_ij*σ(2,1,i_ind)*σ(1,2,j_ind),i_ind,[j_ind]),j_ind),DSum)

Hc = Δc*a'a + η*(a' + a)
Ha = Δa*Σ(σ(2,2,i_ind),i_ind) + DSum
Hi = Σ(g(i_ind)*(a'*σ(1,2,i_ind) + a*σ(2,1,i_ind)),i_ind)
H = Hc + Ha + Hi

J = [a, [σ(1,2,i_ind),σ(1,2,j_ind)] ] 
rates = [κ,Γ_ij]

ops = [a, σ(2,2,k_ind), σ(1,2,k_ind)]
eqs = indexed_meanfield(ops,H,J;rates=rates,order=order)

@test isequal([i_ind,j_ind,k_ind],sort(qc.get_indices_equations(eqs)))
@test isequal([:i,:j,:k],sort(qc.getIndName.(qc.get_indices_equations(eqs))))

@test length(eqs) == 3

ind1 = Index(h,:q,N,ha)
ind2 = Index(h,:r,N,ha)
ind3 = Index(h,:s,N,ha)

eqs_comp = complete(eqs;extra_indices=[ind1,ind2,ind3])
eqs_comp2 = complete(eqs)

@test length(eqs_comp.equations) == length(eqs_comp2.equations)

eqs_ = evaluate(eqs_comp)
eqs_2 = evaluate(eqs_comp2)

@test length(eqs_2) == length(eqs_)

for i = 1:length(eqs_)
    @test length(arguments(eqs_[i].rhs)) == length(arguments(eqs_2[i].rhs))
end

@test length(eqs_) == 18

eqs_4 = indexed_meanfield(ops,H,J;rates=rates,order=4)

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


eqs1 = indexed_meanfield(a_,H,J;rates=rates,order=order) 
eqs2 = indexed_meanfield([a_],H,J;rates=rates,order=order) 
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
eqs_2 = indexed_meanfield(ops_2,H_2,J_2;rates=rates_2,order=order)

q = Index(h,:q,N,ha)
r = Index(h,:r,N2,hc)

extra_indices = [q,r]

eqs_com = complete(eqs_2;extra_indices=extra_indices);
@test length(eqs_com) == 15

@test isequal(sort([k,m,n,l,q,r]),sort(qc.get_all_indices(eqs_com)))

e_1 = evaluate(eqs_com; h=ha,mapping=(N=>5))
e_2 = evaluate(eqs_com; h=hc,mapping=(N2=>6))

@test length(e_1) != length(e_2)
@test !(e_1.equations == e_2.equations)
@test !(e_1.states == e_2.states)

@test sort(qc.get_indices_equations(e_1)) == sort([m,n,r])
@test sort(qc.get_indices_equations(e_2)) == sort([k,l,q])

mapping = Dict(N=>5,N2=>6)
s1 = evaluate(eqs_com; h=[hc,ha],mapping=mapping)
s2 = evaluate(eqs_com; mapping=mapping)

@test length(s1) == length(s2)
@test s1.equations == s2.equations

@test qc.get_indices_equations(s1) == []

@test s1.states == s2.states



end