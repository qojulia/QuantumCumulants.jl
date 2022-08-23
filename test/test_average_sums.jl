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

ind(i) = Index(h,i,N,true)

g(k) = IndexedVariable(:g,k)
σ(i,j,k) = IndexedOperator(Transition(h,:σ,i,j),k)

@test(isequal(average(2*σ(1,2,ind(:k))),2*average(σ(1,2,ind(:k)))))
@test(isequal(average(g(ind(:k))*σ(2,2,ind(:k))),g(ind(:k))*average(σ(2,2,ind(:k)))))
@test(isequal(average(g(ind(:k))),g(ind(:k))))

sum1 = IndexedSingleSum(σ(1,2,ind(:k)),ind(:k))
σn(i,j,k) = NumberedOperator(Transition(h,:σ,i,j),k)
@test(isequal(evalTerm(average(sum1)),average(σn(1,2,1)) + average(σn(1,2,2))))
@test(isequal(σn(1,2,1)+σn(2,1,1),NumberedOperator(Transition(h,:σ,1,2)+Transition(h,:σ,2,1),1)))

@test(isequal(sum1,undo_average(average(sum1))))

#test insertIndex
@test(isequal(σn(2,2,1),insertIndex(σ(2,2,ind(:j)),ind(:j),1)))
@test(isequal(σ(1,2,ind(:j)),insertIndex(σ(1,2,ind(:j)),ind(:k),2)))
@test(isequal(1,insertIndex(1,ind(:k),1)))

sum2 = average(sum1*σ(1,2,ind(:l)))

@test(!isequal(σn(2,2,1),insertIndex(sum2,ind(:j),1)))
@test(isequal(
    indexedAverageSum(average(σ(1,2,ind(:k))),ind(:k)
)))


end

