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

# test scaling keyword in complete
eqs_com_2 = qc.complete(eqs_2;scaling=true,h=1)
eqs_sc_2 = scale(eqs_com2;h=1)

st = scale_term.(eqs_com_2.states;h=1)

@test isequal(length(eqs_com_2.states),length(eqs_sc_2.states))
@test isequal(st,eqs_sc_2.states)

eqs_com_2_2 = qc.complete(eqs_2;scaling=true,h=2)
eqs_sc_2_2 = scale(eqs_com2;h=2)

st_2 = scale_term.(eqs_com_2_2.states;h=2)

@test isequal(length(eqs_com_2_2.states),length(eqs_sc_2_2.states))
@test isequal(st_2,eqs_sc_2_2.states)

eqs_com_2_3 = qc.complete(eqs_2;scaling=true)
eqs_sc_2_3 = scale(eqs_com2)

st_3 = scale_term.(eqs_com_2_3.states)

@test isequal(length(eqs_com_2_3.states),length(eqs_sc_2_3.states))
@test isequal(st_3,eqs_sc_2_3.states)
#

# @test length(s_1) == length(s_2)
@test !(s_1.equations == s_2.equations)
@test !(s_1.states == s_2.states)

@test sort(qc.get_indices_equations(s_1)) == sort([m,n,r])
@test sort(qc.get_indices_equations(s_2)) == sort([k,l,q])

s1 = scale(eqs_com; h=[1,2])
s2 = scale(eqs_com)

s1_ = scale(eqs_com; h=[hc])
s2_ = scale(eqs_com; h=[1])

@test s1_.equations == s2_.equations

se = scale(eqs_com;h=[1])
se2 = qc.evaluate(se;h=[2],limits=(N=>2))

es = qc.evaluate(eqs_com;h=[2],limits=(N=>2))
es2 = scale(es;h=[1])

@test length(se2.equations) == length(es2.equations)
@test isequal(se2.equations[1],es2.equations[1])

@test length(s1) == length(s2)
@test s1.equations == s2.equations

avg = average(σ(2,2,k))
c_avg = conj(avg)

@test operation(qc.inorder!(c_avg)) == conj
@test operation(qc.inorder!(avg)) == qc.sym_average

@test operation(qc.insert_index(avg,k,1)) == qc.sym_average
@test operation(qc.insert_index(c_avg,k,1)) == conj

@test qc.get_indices_equations(s1) == []

@test s1.states == s2.states

@test isequal(scale(∑(average(σ(2,2,k)),k)), N*average(σ(2,2,1)))
@test isequal(scale(∑(average(ai(m)),m)), N2*average(ai(1)))
@test isequal(scale(average(∑(σ(2,1,k)*σ(1,2,l),k))), (N-1)*average(σ(2,1,1)*σ(1,2,2)) + average(σ(2,2,1)))


hc_ = FockSpace(:cavity) 
ha_ = NLevelSpace(:atom, 3)
h_ = hc_ ⊗ ha_
i2 = Index(h_,:i,N,ha_);
j2 = Index(h_,:j,N,ha_);

σ2(i,j,k) = IndexedOperator(Transition(h_,:σ,i,j),k)

Sz_(i) = ∑(σ2(2,2,i) - σ2(3,3,i),i)
Sz2 = average(Sz_(i2)*Sz_(j2))

@test QuantumCumulants.isscaleequal(scale(Sz2),average(N*average(σ2(3,3,1)) + N*σ2(2,2,1)) + (-1+N)^2*average(σ2(2,2,1)*σ2(2,2,2)) +
    (-1+N)^2*average(σ2(3,3,1)*σ2(3,3,2)) -2*(-1+N)^2*average(σ2(2,2,1)*σ2(3,3,2)))

# issue #169 #interacting two level systems
@cnumbers N
h= NLevelSpace(:s,2)
i1 = Index(h,:i1,N,h)
i2 = Index(h,:i2,N,h)
s(α,β,i) = IndexedOperator(Transition(h, :σ, α, β),i)

@test isequal(scale(average(Σ(s(1,2,i1)*s(2,2,i2),i1,[i2]))), (N-1)*average(s(1,2,1)*s(2,2,2)))
@test isequal(scale(average(Σ(s(1,2,i1),i1,[i2]))), (N-1)*average(s(1,2,1))) 
@test isequal(scale(average(Σ(s(1,2,i1),i1))), N*average(s(1,2,1)))

h = NLevelSpace(:spin,2)
@cnumbers N V Ω
order = 1
i1 = Index(h,:i1,N,h)
i2 = Index(h,:i2,N,h)
i = Index(h,:i,N,h)
s(α,β,i) = IndexedOperator(Transition(h, :S, α, β, 1),i)
sp(i1) = (s(2,1,i1))
sm(i1) = (s(1,2,i1))
int_sum = Σ(sp(i1) * sm(i2) + sm(i1) * sp(i2),i1,i2) - Σ(sp(i1) * sm(i1) + sm(i1) * sp(i1),i1)
Hint = V * int_sum + Ω*Σ(sp(i1) + sm(i1),i1)

eqs = meanfield([s(1,2,i)],Hint;order)
eqs_c = complete(eqs)

N_ = 3
V_ = 4.79 / 2
Ω_ = 1.0;

# scale
eqs_sc = scale(eqs_c)
@named sys_sc = ODESystem(eqs_sc)
u0_sc = zeros(ComplexF64, length(eqs_sc))
prob_sc = ODEProblem(sys_sc, u0_sc, (0,2.0), [N,V,Ω].=>[N_,V_,Ω_])
sol_sc = solve(prob_sc,RK4())

# evaluate
eqs_ev = evaluate(eqs_c; limits=(N=>N_))
@named sys_ev = ODESystem(eqs_ev)
u0_ev = zeros(ComplexF64, length(eqs_ev))
prob_ev = ODEProblem(sys_ev, u0_ev, (0,2.0), [V,Ω].=>[V_,Ω_])
sol_ev = solve(prob_ev,RK4())

#### straight forward (no indexing) ###
h_ = tensor([h for i=1:N_]...)
σ(α,β,i) = Transition(h_, Symbol(:σ_,i), α, β, i)
σp(i1) = (σ(2,1,i1))
σm(i1) = (σ(1,2,i1))

int_sum_ = sum(σp(i1) * σm(i2) + σm(i1) * σp(i2) for i1=1:N_ for i2=1:N_) - sum(σp(i1) * σm(i1) + σm(i1) * σp(i1) for i1=1:N_)
Hint_ = V * int_sum_ + Ω*sum(σp(i1) + σm(i1) for i1=1:N_)

eqs_ = meanfield([σ(1,2,1)],Hint_;order)
eqs_c_ = complete(eqs_)

@named sys_ = ODESystem(eqs_c_)
u0_ = zeros(ComplexF64, length(eqs_c_))
prob_ = ODEProblem(sys_, u0_, (0,2.0), [N,V,Ω].=>[N_,V_,Ω_])
sol_ = solve(prob_,RK4())

@test sol_ev[s(1,2,1)][end] ≈ sol_[σ(1,2,1)][end] ≈ sol_sc[s(1,2,1)][end]
@test sol_ev[s(2,2,1)][end] ≈ sol_[σ(2,2,1)][end] ≈ sol_sc[s(2,2,1)][end]

end