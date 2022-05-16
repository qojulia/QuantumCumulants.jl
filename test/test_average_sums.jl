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

@testset "average_sums" begin

N = 2
ha = NLevelSpace(Symbol(:atom),2)
hf = FockSpace(:cavity)
h = hf⊗ha

i_ind = Index(h,:i,N)
j_ind = Index(h,:j,N)
k_ind = Index(h,:k,N)

g(k) = IndexedVariable(:g,k)
σ(i,j,k) = IndexedOperator(Transition(h,:σ,i,j),k)

@test(isequal(average(2*σ(1,2,k_ind)),2*average(σ(1,2,k_ind))))
@test(isequal(average(g(k_ind)*σ(2,2,k_ind)),g(k_ind)*average(σ(2,2,k_ind))))
@test(isequal(average(g(k_ind)),g(k_ind)))

sum1 = IndexedSingleSum(σ(1,2,k_ind),k_ind)
σn(i,j,k) = NumberedOperator(Transition(h,:σ,i,j),k)
@test(isequal(evalTerm(average(sum1)),average(σn(1,2,1)) + average(σn(1,2,2))))

end

