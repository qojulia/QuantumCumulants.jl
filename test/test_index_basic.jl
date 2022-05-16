using QuantumOptics
using OrdinaryDiffEq
using ModelingToolkit
using LinearAlgebra
using Symbolics
using SymbolicUtils
using DifferentialEquations
using Plots
include("../src/indexing.jl")
include("../src/doubleSums.jl")
include("../src/averageSums.jl")
include("../src/indexedMeanfield.jl")
using Test

@testset "index_basic" begin

N = 10
ha = NLevelSpace(Symbol(:atom),2)
hf = FockSpace(:cavity)
h = hf⊗ha

i_ind = Index(h,:i,N)
j_ind = Index(h,:j,N)

@test(!isequal(i_ind,j_ind))
@test(isequal(i_ind,Index(h,:i,10)))

g(k) = IndexedVariable(:g,k)
@test(!isequal(g(i_ind),g(j_ind)))
@test(isequal(g(i_ind),g(Index(h,:i,10))))

σ(i,j,k) = IndexedOperator(Transition(h,:σ,i,j),k)
σ12i = σ(1,2,i_ind)
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

k_ind = Index(h,:k,N)
Γij = DoubleIndexedVariable(:Γ,i_ind,j_ind,true)

@test(isequal(changeIndex(Γij,j_ind,k_ind), DoubleIndexedVariable(:Γ,i_ind,k_ind,true)))
@test(isequal(changeIndex(σ(1,2,j_ind)*σ(1,2,i_ind),j_ind,i_ind),0))
@test(isequal(changeIndex(g(k_ind),k_ind,j_ind),g(j_ind)))

@test(isequal(
    orderByIndex(σ(1,2,k_ind)*σ(1,2,j_ind)*σ(1,2,i_ind),[i_ind]), σ(1,2,i_ind)*σ(1,2,k_ind)*σ(1,2,j_ind)
    ))

@test(isequal(reorder(σ(1,2,k_ind)*σ(1,2,j_ind)*σ(1,2,i_ind),[(i_ind,j_ind)]), SpecialIndexedTerm(σ(1,2,k_ind)*σ(1,2,i_ind)*σ(1,2,j_ind),[(i_ind,j_ind)])))

end

