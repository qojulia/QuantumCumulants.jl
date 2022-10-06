using Test
using QuantumCumulants
using SymbolicUtils
using Symbolics


@testset "double_sums" begin


N = 10
ha = NLevelSpace(Symbol(:atom),2)
hf = FockSpace(:cavity)
h = hf⊗ha

ind(i) = Index(h,i,N,ha)

σ(i,j,k) = IndexedOperator(Transition(h,:σ,i,j),k)
g(k) = IndexedVariable(:g,k)

innerSum = IndexedSingleSum(σ(2,1,ind(:i))*σ(1,2,ind(:j)),ind(:i))
Dsum = IndexedDoubleSum(innerSum,ind(:j),[ind(:i)])
@test(isequal(
    IndexedDoubleSum(innerSum,ind(:j)), IndexedDoubleSum(IndexedSingleSum(σ(2,1,ind(:i))*σ(1,2,ind(:j)),ind(:i),[ind(:j)]),ind(:j)) + IndexedSingleSum(σ(2,2,ind(:j)),ind(:j))
))

@test(isequal(
    IndexedDoubleSum(innerSum,ind(:j)), IndexedDoubleSum(IndexedSingleSum(σ(2,1,ind(:i)),ind(:i))*σ(1,2,ind(:j)),ind(:j))
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


g_ik = DoubleIndexedVariable(:g,i_ind,k_ind,true)

a(k) = IndexedOperator(Destroy(h,:a),k)
σ(i,j,k) = IndexedOperator(Transition(h,Symbol("σ"),i,j),k)

Ssum1 = Σ(g_ik*a(k_ind)*σ(2,1,i_ind),i_ind)
Ssum2 = Σ(g_ik*a(k_ind)'*σ(1,2,i_ind),i_ind)

@test isequal(Σ(conj(g_ik)*a(k_ind)'*σ(1,2,i_ind),i_ind),Ssum1')

@test isequal(Σ(Σ(g_ik*(a(k_ind)*σ(2,1,i_ind) + a(k_ind)'*σ(1,2,i_ind)),i_ind),k_ind),
    Σ(Ssum1+Ssum2,k_ind))

@test Ssum1 isa IndexedSingleSum
@test Σ(Ssum1,k_ind) isa IndexedDoubleSum

H = Σ(Σ(g_ik*(a(k_ind)*σ(2,1,i_ind) + a(k_ind)'*σ(1,2,i_ind)),i_ind),k_ind)

for arg in H.arguments
    @test arg isa IndexedDoubleSum
end

H2 = H*a(l_ind)

@test Ssum1*a(l_ind) isa IndexedSingleSum

DSum1 = H.arguments[1]
Dsum2 = H.arguments[2]

@test isequal(DSum1*a(l_ind),Σ(Ssum1*a(l_ind),k_ind))

@test isequal(Σ(Σ(σ(2,1,i_ind)*σ(1,2,j_ind),i_ind,[j_ind]),j_ind,[i_ind]),Σ(Σ(σ(2,1,i_ind)*σ(1,2,j_ind),i_ind,[j_ind]),j_ind))
@test isequal(Σ(Σ(σ(2,1,i_ind)*σ(1,2,j_ind),i_ind,[j_ind]),j_ind,[i_ind]),Σ(Σ(σ(2,1,i_ind),i_ind,[j_ind])*σ(1,2,j_ind),j_ind))
 
ADsum1 = average(DSum1)
ADsum2 = average(Dsum2)

split0 = splitSums(ADsum1,j_ind,15)
@test isequal(split0,ADsum1)

split1 = splitSums(ADsum1,k_ind,5)
@test split1 isa SymbolicUtils.Mul
@test isequal(5,arguments(split1)[1])
@test isequal(N_modes/5,arguments(split1)[2].metadata.sumIndex.rangeN)

split2 = splitSums(ADsum1,i_ind,5)
@test split2 isa SymbolicUtils.Mul
@test isequal(5,arguments(split2)[1])
@test isequal(N_atoms/5,arguments(split2)[2].metadata.innerSum.metadata.sumIndex.rangeN)

@test isequal(N_modes,arguments(split2)[2].metadata.sumIndex.rangeN)


#for arg in H.arguments
#    @test arg isa IndexedDoubleSum
#end





end