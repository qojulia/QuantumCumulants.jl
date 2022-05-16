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

@testset "double_sums" begin


N = 10
ha = NLevelSpace(Symbol(:atom),2)
hf = FockSpace(:cavity)
h = hf⊗ha

i_ind = Index(h,:i,N)
j_ind = Index(h,:j,N)

σ(i,j,k) = IndexedOperator(Transition(h,:σ,i,j),k)
g(k) = IndexedVariable(:g,k)

innerSum = IndexedSingleSum(σ(2,1,i_ind)*σ(1,2,j_ind),i_ind)
@test(isequal(
    IndexedDoubleSum(innerSum,j_ind), IndexedDoubleSum(IndexedSingleSum(σ(2,1,i_ind)*σ(1,2,j_ind),i_ind,[j_ind]),j_ind) + IndexedSingleSum(σ(2,2,j_ind),j_ind)
))

end