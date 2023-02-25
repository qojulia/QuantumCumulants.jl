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

end
