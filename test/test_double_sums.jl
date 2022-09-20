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

@test(isequal(hilbert(Dsum),hilbert(innerSum)))
@test(isequal(hilbert(Dsum),hilbert(ind(:j))))
@test(isequal(hilbert(Dsum),h))


end