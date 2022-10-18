using Test
using QuantumCumulants
using SymbolicUtils
using Symbolics

const qc=QuantumCumulants

@testset "index_basic" begin

N = 10
ha = NLevelSpace(Symbol(:atom),2)
hf = FockSpace(:cavity)
h = hf⊗ha

indT(i) = Index(h,i,N,ha) #transition index
indF(i) = Index(h,i,N,hf) #fock index
i_ind = indT(:i)
j_ind = indT(:j)

ind(a) = indT(a)

@test(!isequal(indT(:i),indT(:j)))
@test(!isequal(indT(:i),indF(:j)))
@test(!isequal(indT(:i),indF(:i)))

@test(isequal(indT(:i),Index(h,:i,10,ha)))

g(k) = IndexedVariable(:g,k)
@test(!isequal(g(indT(:i)),g(indT(:j))))
@test(isequal(g(indT(:i)),g(Index(h,:i,10,ha))))

σ(i,j,k) = IndexedOperator(Transition(h,:σ,i,j),k)
σ12i = σ(1,2,indT(:i))
@test(isequal(σ12i,σ(1,2,i_ind)))
@test(!isequal(σ12i,σ(2,2,i_ind)))
@test(!isequal(σ12i,σ(1,2,j_ind)))

@test(isequal(0,σ12i*σ(1,2,i_ind)))
@test(isequal(σ(2,2,i_ind),σ(2,1,i_ind)*σ12i))

#@test(isequal(σ(2,2,i_ind)+σ(1,2,j_ind),σ(1,2,j_ind)+σ(2,2,i_ind)))
#apperently QAdd isequal function is dependant in order of terms inside the addition (?)

@test(isequal(adjoint(σ(1,2,i_ind)),σ(2,1,i_ind)))


@qnumbers a::Destroy(h)
sum1 = IndexedSingleSum(σ(1,2,i_ind)*a',i_ind)
sum2 = IndexedSingleSum(σ(2,1,i_ind)*a,i_ind)
@test(isequal(adjoint(sum1),sum2))

sum3 = IndexedSingleSum(a'*σ(1,2,i_ind) + a*σ(2,1,i_ind),i_ind)
@test(isequal(sum3,(sum1+sum2)))
@test(isequal(acts_on(σ12i),2))
@test(i_ind < j_ind)

@test isequal(0,Σ(0,i_ind))
@test isequal(0,Σ(σ(2,1,i_ind)*σ(2,1,i_ind),i_ind))

k_ind = indT(:k)
Γij = DoubleIndexedVariable(:Γ,i_ind,j_ind)

@test(isequal(change_index(Γij,j_ind,k_ind), DoubleIndexedVariable(:Γ,i_ind,k_ind)))
@test(isequal(change_index(σ(1,2,j_ind)*σ(1,2,i_ind),j_ind,i_ind),0))
@test(isequal(change_index(g(k_ind),k_ind,j_ind),g(j_ind)))

@test(isequal(
    order_by_index(σ(1,2,k_ind)*σ(1,2,j_ind)*σ(1,2,i_ind),[i_ind]), σ(1,2,i_ind)*σ(1,2,k_ind)*σ(1,2,j_ind)
    ))

@test(isequal(
    reorder(σ(1,2,k_ind)*σ(1,2,j_ind)*σ(1,2,i_ind),[(i_ind,j_ind)]), 
    SpecialIndexedTerm(σ(1,2,k_ind)*σ(1,2,i_ind)*σ(1,2,j_ind),[(i_ind,j_ind)])
))
@test(isequal(
    σ(1,2,k_ind) * sum1, simplify(IndexedSingleSum(σ(1,2,k_ind)*σ(1,2,i_ind)*a',i_ind))
))
@test(isequal(
    simplify(σ(2,1,k_ind) * sum1), simplify(IndexedSingleSum(σ(2,1,k_ind)*σ(1,2,i_ind)*a',i_ind,[k_ind]) + a'*σ(2,2,k_ind))
))
innerSum = IndexedSingleSum(σ(2,1,i_ind)*σ(1,2,j_ind),i_ind)
@test(isequal(
    IndexedDoubleSum(innerSum,j_ind), IndexedDoubleSum(IndexedSingleSum(σ(2,1,i_ind)*σ(1,2,j_ind),i_ind,[j_ind]),j_ind) + IndexedSingleSum(σ(2,2,j_ind),j_ind)
))
@test(isequal(SymbolicUtils.arguments(σ(1,2,indT(:i))*a'),SymbolicUtils.arguments(sum1)))

@test isequal(N*g(ind(:j)),Σ(g(ind(:j)),ind(:i)))
@test Σ(g(ind(:j)),ind(:j)) isa qc.IndexedSingleSum

@test isequal(N*Γij,Σ(Γij,ind(:k)))
@test Σ(Γij,ind(:i)) isa qc.IndexedSingleSum

@test (sum1 + a') isa qc.QAdd
@test (sum1 + σ(1,2,i_ind)) isa qc.QAdd

@test (σ(2,1,j_ind) + σ(1,2,i_ind)) isa qc.QAdd
@test (sum1 + g(i_ind)) isa qc.QAdd
@test isequal(sum1 + g(i_ind), g(i_ind) + sum1)
@test (a + σ(1,2,i_ind)) isa qc.QAdd
@test (σ(1,2,i_ind)+a) isa qc.QAdd

qadd = a + a'
@test length((qadd + sum1).arguments) == 3
@test isequal((sum1+qadd),(qadd + sum1))

@test isequal(sum1 + g(i_ind),g(i_ind) + sum1)

@test length((qadd + σ(1,2,i_ind)).arguments) == 3
@test isequal((σ(1,2,i_ind)+qadd),(qadd + σ(1,2,i_ind)))

@test length((qadd + g(i_ind)).arguments) == 3
@test isequal((g(i_ind)+qadd),(qadd + g(i_ind)))

qmul = a'*a
@test sum1+qmul isa qc.QAdd
@test isequal(sum1+qmul,qmul+sum1)
@test isequal((σ(1,2,i_ind)+qmul),(qmul + σ(1,2,i_ind)))
@test isequal((g(i_ind)+qmul),(qmul + g(i_ind)))
@test isequal(g(i_ind) + σ(1,2,j_ind),σ(1,2,j_ind) + g(i_ind))

specTerm = qc.SpecialIndexedTerm(σ(1,2,i_ind)*σ(1,2,j_ind),[(i_ind,j_ind)])
@test isequal((sum1+specTerm),(specTerm + sum1))
@test isequal((σ(1,2,i_ind)+specTerm),(specTerm + σ(1,2,i_ind)))
@test isequal((g(i_ind)+specTerm),(specTerm + g(i_ind)))
@test isequal((specTerm+qadd),(qadd + specTerm))
@test isequal((specTerm+qmul),(qmul + specTerm))
@test isequal((specTerm+2),(2+ specTerm))

@test isequal(-σ(1,2,i_ind),-1*σ(1,2,i_ind))
@test isequal(-g(i_ind),-1*g(i_ind))

@test isequal(g(i_ind) + a,a + g(i_ind))

@test isequal(qadd+g(i_ind),g(i_ind)+qadd)

@test g(i_ind)*a isa qc.QMul
@test (g(i_ind)*a).args_nc == [a]
@test g(i_ind)*a' isa qc.QMul
@test (g(i_ind)*a').args_nc == [a']

@test a*g(i_ind) isa qc.QMul
@test (a*g(i_ind)).args_nc == [a]
@test a'*g(i_ind) isa qc.QMul
@test (a'*g(i_ind)).args_nc == [a']

@test σ(1,2,i_ind)*g(i_ind) isa qc.QMul
@test (σ(1,2,i_ind)*g(i_ind)).args_nc == [σ(1,2,i_ind)]
@test g(i_ind)*σ(1,2,i_ind) isa qc.QMul
@test (g(i_ind)*σ(1,2,i_ind)).args_nc == [σ(1,2,i_ind)]

@test length((qmul*g(i_ind)).args_nc) == 2
@test isequal((qmul*g(i_ind)).arg_c,g(i_ind))
@test length((g(i_ind)*qmul).args_nc) == 2
@test isequal((g(i_ind)*qmul).arg_c,g(i_ind))

@test isequal(acts_on(σ(1,2,i_ind)),2)

ai(k) = IndexedOperator(Destroy(h,:a),k)

@test isequal((ai(indF(:m))*ai(indF(:m))'),ai(indF(:m))'*ai(indF(:m)) + 1)


specTerm = qc.SpecialIndexedTerm(σ(1,2,i_ind)*σ(1,2,j_ind),[(i_ind,j_ind)])
asdf = specTerm*σ(1,2,k_ind)
asdf2 = σ(1,2,k_ind)*specTerm

@test isequal(asdf,asdf2)
@test isequal(specTerm*qmul,qmul*specTerm)
@test isequal(qadd*specTerm,specTerm*qadd)
@test isequal(2*specTerm,specTerm*2)

@test isequal(commutator(σ(1,2,i_ind),σ(2,1,i_ind)),σ(1,2,i_ind)*σ(2,1,i_ind) - σ(2,1,i_ind)*σ(1,2,i_ind))
@test isequal(simplify(commutator(σ(1,2,i_ind),qadd)),0)
@test isequal(simplify(commutator(σ(1,2,i_ind),qmul)),0)

Ωij = DoubleIndexedVariable(:Ω,i_ind,j_ind;identical=false)

@test change_index(Ωij,i_ind,j_ind) == 0
@test reorder(qc.QAdd([]),[(i_ind,j_ind)]) == 0
@test reorder(qc.QAdd([0]),[(i_ind,j_ind)]) == 0
@test reorder(qc.QAdd([σ(1,2,i_ind),σ(2,1,j_ind)]),[(i_ind,j_ind)]) isa qc.QAdd

@test reorder(average(qc.QAdd([0])),[(i_ind,j_ind)]) == 0

@test isequal(NumberedOperator(Transition(h,:σ,1,2),1),σ(1,2,1))

@test isequal(∑(σ(1,2,i_ind),i_ind),Σ(σ(1,2,i_ind),i_ind))
@test isequal(∑(σ(1,2,i_ind)*σ(2,1,j_ind),i_ind,j_ind),Σ(σ(1,2,i_ind)*σ(2,1,j_ind),i_ind,j_ind))

@test isequal([i_ind,j_ind],qc.getIndices(σ(1,2,i_ind) + σ(2,1,j_ind)))
@test isequal([i_ind,j_ind],sort(qc.getIndices(average(σ(1,2,i_ind)) + 3 + average(σ(2,1,j_ind)))))

@test isequal(IndexedVariable(:Ω,1,2),qc.DoubleNumberedVariable(:Ω,1,2))
@test isequal(IndexedVariable(:Ω,2),qc.SingleNumberedVariable(:Ω,2))

@test isequal(σ(1,2,1.0),σ(1,2,1))
@test isequal(σ(1,2,1.9),σ(1,2,2))

@test isequal(g(1.1),qc.SingleNumberedVariable(:g,1))
@test isequal(IndexedVariable(:Ω,1.1,2.1),qc.DoubleNumberedVariable(:Ω,1,2))

hc = FockSpace(:cavity)
hf = FockSpace(:filter)

h = hc ⊗ hf

i = Index(h,:i,N,hf)
j = Index(h,:j,N,hf)
k = Index(h,:k,N,hf)

xij = IndexedVariable(:x,i,j)


@qnumbers a_::Destroy(h,1)
b(k) = IndexedOperator(Destroy(h,:b,2), k)

@test reorder(b(i)*b(k)*b(i)',[(i,k)]) isa qc.QAdd
@test reorder(b(i)'*b(i)*b(k),[(i,k)]) isa qc.SpecialIndexedTerm
@test isequal(reorder(b(i)*b(k)*b(i)',[(i,k)]),reorder(b(i)'*b(i)*b(k),[(i,k)]) + reorder(b(k),[(i,k)]))

end

