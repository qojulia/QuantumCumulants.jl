using Test
using QuantumCumulants
using SymbolicUtils
using Symbolics

const qc = QuantumCumulants

@testset "indexed_scale" begin

#test a system, where scaling is done individually per hilbertspace

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
eqs_com2 = qc.complete(eqs_2);
@test length(eqs_com) == 15
@test length(eqs_com) == length(eqs_com2)

s_1 = scale(eqs_com; h=2)
s_2 = scale(eqs_com; h=1)

@test length(s_1) == length(s_2)
@test !(s_1.equations == s_2.equations)
@test !(s_1.states == s_2.states)

@test sort(qc.get_indices_equations(s_1)) == sort([m,n,r])
@test sort(qc.get_indices_equations(s_2)) == sort([k,l,q])

s1 = scale(eqs_com; h=[1,2])
s2 = scale(eqs_com)

@test length(s1) == length(s2)
@test s1.equations == s2.equations

@test qc.get_indices_equations(s1) == []

@test s1.states == s2.states

end