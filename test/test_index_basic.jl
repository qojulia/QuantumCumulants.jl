using Test
using QuantumCumulants
using QuantumOpticsBase
using SymbolicUtils
using Symbolics
using OrdinaryDiffEq
using UUIDs

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

    # @test isequal(simplify(σ(2,2,i_ind)+σ(1,2,j_ind)),simplify(σ(1,2,j_ind)+σ(2,2,i_ind)))
    @test(isequal(adjoint(σ(1,2,i_ind)),σ(2,1,i_ind)))

    @test isequal(σ(2,1,i_ind)*σ(1,2,i_ind), σ(2,2, i_ind))

    ex = σ(2,1,i_ind)*σ(1,2,j_ind)
    s1 = σ(2,1,i_ind)
    s2 = σ(1,2,j_ind)
    id = uuid4()
    push!(s1.merge_events, id)
    push!(s2.merge_events, id)
    @test isequal(ex, σ(2,2,i_ind)*(i_ind == j_ind) + (1 - (i_ind==j_ind)) * s1*s2)

    a = Destroy(h,:a)
    a_indexed(i) = IndexedOperator(a, i)
    r_ind = indF(:r)
    s_ind = indF(:s)
    @test isequal(a_indexed(r_ind) * a_indexed(r_ind)', a_indexed(r_ind)' * a_indexed(r_ind) + 1)
    @test isequal(a_indexed(r_ind) * a_indexed(s_ind)', a_indexed(r_ind)' * a_indexed(r_ind) + r_ind == j_ind)

end


# @testset "sums" begin

    N = 10
    ha = NLevelSpace(Symbol(:atom),2)
    hf = FockSpace(:cavity)
    h = hf⊗ha

    indT(i) = Index(h,i,N,ha) #transition index
    indF(i) = Index(h,i,N,hf) #fock index
    i_ind = indT(:i)
    j_ind = indT(:j)

    ind(a) = indT(a)

    a = Destroy(h,:a)
    σ(i,j,k) = IndexedOperator(Transition(h,:σ,i,j),k)
    σ12i = σ(1,2,indT(:i))

sum1 = Sum(σ(1,2,i_ind)*a',i_ind)
sum2 = Sum(σ(2,1,i_ind)*a,i_ind)
@test(isequal(adjoint(sum1),sum2))

sum3 = Sum(a'*σ(1,2,i_ind) + a*σ(2,1,i_ind),i_ind)
@test(isequal(sum3,(sum1+sum2)))
@test(isequal(acts_on(σ12i),2))
# @test(i_ind < j_ind)

@test isequal(0,Σ(0,i_ind))
@test isequal(0,Σ(σ(2,1,i_ind)*σ(2,1,i_ind),i_ind))

k_ind = indT(:k)
Γij = IndexedParameter(:Γ,i_ind,j_ind)
g = IndexedParameter(:g)

@test(isequal(change_index(Γij,j_ind,k_ind), IndexedParameter(:Γ,i_ind,k_ind)))
@test(isequal(change_index(σ(1,2,j_ind)*σ(1,2,i_ind),j_ind,i_ind),0))
@test(isequal(change_index(g(k_ind),k_ind,j_ind),g(j_ind)))
@test isequal(change_index(Σ(2g(i_ind),i_ind), i_ind, j_ind), Σ(2g(j_ind),j_ind))

@test(isequal(
    order_by_index(σ(1,2,k_ind)*σ(1,2,j_ind)*σ(1,2,i_ind),[i_ind]), σ(1,2,i_ind)*σ(1,2,k_ind)*σ(1,2,j_ind)
    ))

@test(isequal(
    reorder(σ(1,2,k_ind)*σ(1,2,j_ind)*σ(1,2,i_ind),[(i_ind,j_ind)]),
    SpecialIndexedTerm(σ(1,2,k_ind)*σ(1,2,i_ind)*σ(1,2,j_ind),[(i_ind,j_ind)])
))
@test(isequal(σ(1,2,k_ind) * sum1, simplify(Sum(σ(1,2,k_ind)*σ(1,2,i_ind)*a',i_ind))
))
σ(1,2,k_ind) * sum1
qqq = simplify(Sum(σ(1,2,k_ind)*σ(1,2,i_ind)*a',i_ind))
# qqq = Sum(σ(1,2,k_ind)*σ(1,2,i_ind)*a',i_ind)
QuantumCumulants.get_indices(qqq)
SymbolicUtils._iszero(a'*σ(1,2,i_ind)*σ(1,2,k_ind))

@test(isequal(simplify(σ(2,1,k_ind) * sum1), simplify(Sum(σ(2,1,k_ind)*σ(1,2,i_ind)*a',i_ind,[k_ind]) + a'*σ(2,2,k_ind))
))
innerSum = Sum(σ(2,1,i_ind)*σ(1,2,j_ind),i_ind)
@test(isequal(
    DoubleSum(innerSum,j_ind), DoubleSum(Sum(σ(2,1,i_ind)*σ(1,2,j_ind),i_ind,[j_ind]),j_ind) + Sum(σ(2,2,j_ind),j_ind)
))
@test(isequal(SymbolicUtils.arguments(σ(1,2,indT(:i))*a'),SymbolicUtils.arguments(sum1)))

@test isequal(N*g(ind(:j)),Σ(g(ind(:j)),ind(:i)))
@test Σ(g(ind(:j)),ind(:j)) isa qc.Sum

@test isequal(N*Γij,Σ(Γij,ind(:k)))
@test Σ(Γij,ind(:i)) isa qc.Sum

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
#@test isequal(sum1+qmul,qmul+sum1)
#@test isequal((σ(1,2,i_ind)+qmul),(qmul + σ(1,2,i_ind)))
@test isequal((g(i_ind)+qmul),(qmul + g(i_ind)))
@test isequal(g(i_ind) + σ(1,2,j_ind),σ(1,2,j_ind) + g(i_ind))

specTerm = qc.SpecialIndexedTerm(σ(1,2,i_ind)*σ(1,2,j_ind),[(i_ind,j_ind)])
#@test isequal((sum1+specTerm),(specTerm + sum1))
#@test isequal((σ(1,2,i_ind)+specTerm),(specTerm + σ(1,2,i_ind)))
@test isequal((g(i_ind)+specTerm),(specTerm + g(i_ind)))
@test isequal((specTerm+qadd),(qadd + specTerm))
#@test isequal((specTerm+qmul),(qmul + specTerm))
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

@test isequal([i_ind,j_ind],qc.get_indices(σ(1,2,i_ind) + σ(2,1,j_ind)))
@test isequal([i_ind,j_ind],sort(qc.get_indices(average(σ(1,2,i_ind)) + 3 + average(σ(2,1,j_ind)))))

@test isequal(IndexedVariable(:Ω,1,2),qc.DoubleNumberedVariable(:Ω,1,2))
@test isequal(IndexedVariable(:Ω,2),qc.SingleNumberedVariable(:Ω,2))

# @test isequal(σ(1,2,1.0),σ(1,2,1))
# @test isequal(σ(1,2,1.9),σ(1,2,2))
#
# @test isequal(g(1.1),qc.SingleNumberedVariable(:g,1))
# @test isequal(IndexedVariable(:Ω,1.1,2.1),qc.DoubleNumberedVariable(:Ω,1,2))

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

# Test fock basis conversion
hfock = FockSpace(:fock)
bfock = FockBasis(3)

hnlevel = NLevelSpace(:nlevel,2)
bnlevel = NLevelBasis(2)

h_ = hnlevel ⊗ hfock

N_n = 4
N_f = 2
ranges = [N_n,N_f]

b_1 = ⊗([bnlevel for i = 1:N_n]...)
b_2 = ⊗([bfock for i =1:N_f]...)
b_ = b_1 ⊗ b_2

i1 = Index(h_,:i1,N_n,hnlevel)
j1 = Index(h_,:j1,N_f,hfock)

ai(i) = IndexedOperator(Destroy(h_,:a),i)
σi(i,j,k) = IndexedOperator(Transition(h_,:σ,i,j),k)

@test to_numeric(ai(1),b_;ranges=ranges) == LazyTensor(b_, [5], (destroy(bfock),))
@test to_numeric(ai(2),b_;ranges=ranges) == LazyTensor(b_, [6], (destroy(bfock),))
@test to_numeric(σi(1,2,4),b_;ranges=ranges) == LazyTensor(b_, [4], (QuantumOpticsBase.transition(bnlevel,1,2),))
@test_throws MethodError to_numeric(σi(1,2,5),b_;ranges=ranges)

ai2(i) = IndexedOperator(Destroy(hfock,:a),i)
@test to_numeric(ai2(1),b_2;ranges=[2]) isa LazyTensor
@test to_numeric(ai2(2),b_2;ranges=[2]) isa LazyTensor
@test_throws BoundsError to_numeric(ai2(3),b_2;ranges=[2])

# Indices and only one HilbertSpace
h = NLevelSpace(:atom,2)
i = Index(h,:i,N,h)
σ(x,y,k) = IndexedOperator(Transition(h,:σ,x,y),k)
i.hilb
@test isa(i.hilb,ProductSpace) == false
@test σ(1,2,i) == (σ(1,2,i)')'

# Multiplication IndexedOperator*Transition
hc = NLevelSpace(:cavity, 3)
ha = NLevelSpace(:atom,2)
h = hc ⊗ ha
@cnumbers N α
i = Index(h,:i,N,ha)
S(x,y) = Transition(h,:S,x,y,1)
σ(x,y,k) = IndexedOperator(Transition(h,:σ,x,y,2),k)
@test S(2,1)*σ(1,2,i) isa QuantumCumulants.QMul
@test σ(1,2,i)*S(2,1) isa QuantumCumulants.QMul
@test σ(1,2,2)*S(2,1) isa QuantumCumulants.QMul
@test S(2,1)*σ(1,2,3) isa QuantumCumulants.QMul
@test σ(1,2,2) isa NumberedOperator
@test isequal(S(2,1)*σ(1,2,i), σ(1,2,i)*S(2,1))

j = Index(h,:j,N,hc)
j2 = Index(h,:j,N,ha)
ranges1 = [1:10,1:5]

arr = qc.create_index_arrays([i,j],ranges1)
@test isequal(vec(collect(collect(Iterators.product(ranges1...)))),arr)
arr = qc.create_index_arrays([i],[1:10])
@test isequal(1:10,arr)

@test isequal(qc.inorder!(σ(2,1,1)*σ(2,2,2)*σ(1,2,1)),σ(2,2,1)*σ(2,2,2))
@test isequal(qc.inadjoint(σ(2,1,1)*σ(2,2,2)*σ(1,2,1)),σ(2,2,1)*σ(2,2,2))
@test isequal(qc._inconj(average(σ(2,1,1)*σ(2,2,2)*σ(1,2,1))),(average(σ(2,2,1)*σ(2,2,2))))
@test qc.ismergeable(σ(2,1,5),σ(1,2,5))

# issue 188
gi = IndexedVariable(:g, i)
@test isa(∑(5gi,i), Sum)
@test isa(∑(gi*α,i), Sum)
@test isequal(∑(α,i), N*α)
@test isequal(∑(5α,i), 5*N*α)

# end
