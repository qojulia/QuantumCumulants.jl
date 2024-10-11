using Test
using QuantumCumulants
using SymbolicUtils
using Symbolics

const qc = QuantumCumulants

@testset "double_sums" begin


N = 10
ha = NLevelSpace(Symbol(:atom),2)
hf = FockSpace(:cavity)
h = hf⊗ha

ind(i) = Index(h,i,N,ha)

σ(i,j,k) = IndexedOperator(Transition(h,:σ,i,j),k)
g(k) = IndexedVariable(:g,k)

innerSum = SingleSum(σ(2,1,ind(:i))*σ(1,2,ind(:j)),ind(:i))
Dsum = DoubleSum(innerSum,ind(:j),[ind(:i)])
@test(isequal(
    DoubleSum(innerSum,ind(:j)), DoubleSum(SingleSum(σ(2,1,ind(:i))*σ(1,2,ind(:j)),ind(:i),[ind(:j)]),ind(:j)) + SingleSum(σ(2,2,ind(:j)),ind(:j))
))

@test(isequal(
    DoubleSum(innerSum,ind(:j)), DoubleSum(SingleSum(σ(2,1,ind(:i)),ind(:i))*σ(1,2,ind(:j)),ind(:j))
))

N_atoms = 4
N_modes = 2

hf = FockSpace(:cavity)
ha = NLevelSpace(Symbol(:atom),2)
h = hf ⊗ ha

i_ind = Index(h,:i,N_atoms,ha)
j_ind = Index(h,:j,N_atoms,ha)
k_ind = Index(h,:k,N_modes,hf)
l_ind = Index(h,:l,N_modes,hf)


g_ik = IndexedVariable(:g,i_ind,k_ind)

a(k) = IndexedOperator(Destroy(h,:a),k)
σ(i,j,k) = IndexedOperator(Transition(h,Symbol("σ"),i,j),k)

Ssum1 = Σ(g_ik*a(k_ind)*σ(2,1,i_ind),i_ind)
Ssum2 = Σ(g_ik*a(k_ind)'*σ(1,2,i_ind),i_ind)

### issue 221 (DoubleSum)
@cnumbers c1 N1
i_ind2 = Index(h,:i,N1,ha)
j_ind2 = Index(h,:j,N1,ha)
@test isequal(simplify(Σ(-σ(2,2,i_ind),i_ind,j_ind)),simplify(Σ(-σ(2,2,i_ind),i_ind)*4))
@test isequal(simplify(Σ(3*σ(2,2,i_ind),i_ind,j_ind)),simplify(Σ(3*σ(2,2,i_ind),i_ind)*4))
@test isequal(simplify(Σ(c1*σ(2,2,i_ind),i_ind,j_ind)),simplify(Σ(c1*σ(2,2,i_ind),i_ind)*4))
@test isequal(simplify(Σ(-σ(2,2,i_ind2),i_ind2,j_ind2)),simplify(Σ( (1-N1)*σ(2,2,i_ind2) ,i_ind2)) - Σ( σ(2,2,i_ind2), i_ind2))
@test isequal(simplify(Σ(3*σ(2,2,i_ind2),i_ind2,j_ind2)), simplify(Σ( 3*(N1-1)*σ(2,2,i_ind2) ,i_ind2)) + 3*Σ( σ(2,2,i_ind2), i_ind2))
@test isequal(simplify(Σ(c1*σ(2,2,i_ind2),i_ind2,j_ind2)),simplify(c1*Σ( σ(2,2,i_ind2), i_ind2) + Σ( c1*(N1-1)*σ(2,2,i_ind2) ,i_ind2)) )
###

@test isequal(Σ(conj(g_ik)*a(k_ind)'*σ(1,2,i_ind),i_ind),Ssum1')

@test isequal(Σ(Σ(g_ik*(a(k_ind)*σ(2,1,i_ind) + a(k_ind)'*σ(1,2,i_ind)),i_ind),k_ind),
    Σ(Ssum1+Ssum2,k_ind))

@test Ssum1 isa SingleSum
@test Σ(Ssum1,k_ind) isa DoubleSum

H = Σ(Σ(g_ik*(a(k_ind)*σ(2,1,i_ind) + a(k_ind)'*σ(1,2,i_ind)),i_ind),k_ind)

for arg in H.arguments
    @test arg isa DoubleSum
end

H2 = H*a(l_ind)

@test Ssum1*a(l_ind) isa SingleSum

DSum1 = H.arguments[1]
Dsum2 = H.arguments[2]

@test isequal(DSum1*a(l_ind),Σ(Ssum1*a(l_ind),k_ind))

@test isequal(Σ(Σ(σ(2,1,i_ind)*σ(1,2,j_ind),i_ind,[j_ind]),j_ind,[i_ind]),Σ(Σ(σ(2,1,i_ind)*σ(1,2,j_ind),i_ind,[j_ind]),j_ind))
@test isequal(Σ(Σ(σ(2,1,i_ind)*σ(1,2,j_ind),i_ind,[j_ind]),j_ind,[i_ind]),Σ(Σ(σ(2,1,i_ind),i_ind,[j_ind])*σ(1,2,j_ind),j_ind))

ADsum1 = average(DSum1)
ADsum2 = average(Dsum2)

split0 = split_sums(ADsum1,j_ind,15)
@test isequal(split0,ADsum1)

split1 = split_sums(ADsum1,k_ind,5)
@test split1 isa SymbolicUtils.BasicSymbolic && operation(split1) === *
@test isequal(5,arguments(split1)[1])
@test isequal(N_modes/5,arguments(split1)[2].metadata[qc.IndexedAverageDoubleSum].sum_index.range)

split2 = split_sums(ADsum1,i_ind,5)
@test split2 isa SymbolicUtils.BasicSymbolic && operation(split2) === *
@test isequal(5,arguments(split2)[1])
@test isequal(N_atoms/5,arguments(split2)[2].metadata[qc.IndexedAverageDoubleSum].innerSum.metadata[qc.IndexedAverageSum].sum_index.range)

@test isequal(N_modes,arguments(split2)[2].metadata[qc.IndexedAverageDoubleSum].sum_index.range)

innerSum = Σ(σ(1,2,i_ind)*σ(2,1,j_ind),i_ind)
DSum = Σ(innerSum,j_ind)
@test innerSum isa qc.QAdd
@test Dsum isa qc.QAdd

@test DSum isa qc.QAdd

@test Dsum.arguments[1] isa qc.DoubleSum
@test Dsum.arguments[2] isa qc.SingleSum

@test isequal(Σ(Σ(σ(1,2,i_ind)*σ(2,1,j_ind),i_ind,[j_ind]),j_ind) + Σ(σ(1,2,j_ind)*σ(2,1,j_ind),j_ind),DSum)

sum1 = Σ(σ(1,2,i_ind),i_ind)
sum2 = Σ(σ(2,1,i_ind),i_ind)

mul1 = *(sum1,sum2;ind=j_ind)

sum2_ = Σ(σ(2,1,j_ind),j_ind)
mul2 = sum1*sum2_

@test mul2 isa qc.QAdd
@test mul1 isa qc.QAdd

@test isequal(mul1,mul2)

# Double indexed variable
@cnumbers N
i = Index(h,:i,N,ha)
j = Index(h,:j,N,ha)
Γ(i,j) = IndexedVariable(:Γ,i,j)
Ω(i,j) = IndexedVariable(:Ω,i,j;identical=false)
@test iszero(Ω(i,i))
@test iszero(Ω(2,2))
@test !isequal(Ω(2,3),0)
@test !isequal(Γ(i,i),0)
@test !isequal(Γ(2,2),0)

# complete and get_indices
h=NLevelSpace(:spin,2)
@cnumbers N
i1 = Index(h,:i1,N,h)
i2 = Index(h,:i2,N,h)
i = Index(h,:i,N,h)
s(α,β,i) = IndexedOperator(Transition(h, :S, α, β, 1),i)
Hint = Σ(s(2,1,i1) * s(1,2,i2) + s(1,2,i1) * s(2,1,i2),i1,i2) - Σ(s(2,1,i1) * s(1,2,i1) + s(1,2,i1) * s(2,1,i1),i1)
#
mf = meanfield([s(1,2,i)],Hint,order=2)
mf_c = complete(mf)
length(mf_c.states) == 7
#
Ω(i,j) = IndexedVariable(:Ω,i,j)
Hint2 = Σ(Ω(i1,i2)*(s(2,1,i1) * s(1,2,i1) + s(1,2,i1) * s(2,1,i2)),i1,i2)
mf2 = meanfield([s(1,2,i)],Hint,order=2)
mf2_c = complete(mf)
isequal(mf2_c.states, mf_c.states)
#
s1 = Σ(Ω(i1,i2),i1,i2)
isequal(sort(qc.get_indices(Hint)), sort([i2, i1]))
g(i) = IndexedVariable(:g, i)
s2 = Σ(g(i1),i1)
isequal(qc.get_indices(s2),[i1])


### issue 223
ha = NLevelSpace(:atom,2)
σ(α,β,i) = IndexedOperator(Transition(ha, :σ, α, β),i)
@cnumbers N g
i = Index(ha,:i,N,1)
j = Index(ha,:j,N,1)
k = Index(ha,:k,N,1)
#
H = Σ(2*σ(2,2,j),i,j)
H_ji = Σ(2*σ(2,2,j),j,i)
H_s = simplify(H)
H_g = Σ(g*σ(2,2,j),i,j)
H_ji_g = Σ(g*σ(2,2,j),j,i)
H_ji_g_s = simplify(Σ(g*σ(2,2,j),j,i))

dict_N = Dict(N => 10)
sub_dict(x) = simplify(substitute(x, dict_N))
#
@test isequal(sub_dict(simplify(commutator(H,σ(2,1,k)))), 20*σ(2,1,k))
@test isequal(sub_dict(simplify(commutator(H_s,σ(2,1,k)))), 20*σ(2,1,k))
#
@test isequal(sub_dict(simplify(commutator(H,σ(1,2,k)))), -20*σ(1,2,k))
@test isequal(sub_dict(simplify(commutator(H_s,σ(1,2,k)))), -20*σ(1,2,k))

@test isequal(sub_dict(simplify(commutator(H_g,σ(2,1,k)))), 10*g*σ(2,1,k))
@test isequal(sub_dict(simplify(commutator(H_ji_g,σ(2,1,k)))), 10*g*σ(2,1,k))
@test isequal(sub_dict(simplify(commutator(H_ji_g_s,σ(2,1,k)))), 10*g*σ(2,1,k))
#
@test isequal(sub_dict(simplify(commutator(H_g,σ(1,2,k)))), -10g*σ(1,2,k))
@test isequal(sub_dict(simplify(commutator(H_ji_g,σ(1,2,k)))), -10g*σ(1,2,k))
@test isequal(sub_dict(simplify(commutator(H_ji_g_s,σ(1,2,k)))), -10g*σ(1,2,k))

end
