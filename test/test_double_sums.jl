using Test

import SymbolicUtils
import SymbolicUtils: substitute

import Symbolics
import TermInterface

import SciMLBase

import ModelingToolkit
const MTK = ModelingToolkit

using Combinatorics: partitions, combinations
using LinearAlgebra

using QuantumOpticsBase
import QuantumOpticsBase: ⊗, tensor

const NO_METADATA = SymbolicUtils.NO_METADATA

source_metadata(source, name) = 
    Base.ImmutableDict{DataType, Any}(Symbolics.VariableSource, (source, name))

include("../src/hilbertspace.jl")
include("../src/qnumber.jl")
include("../src/cnumber.jl")
include("../src/fock.jl")
include("../src/nlevel.jl")
include("../src/equations.jl")
include("../src/meanfield.jl")
include("../src/average.jl")
include("../src/utils.jl")
include("../src/diffeq.jl")
include("../src/correlation.jl")
include("../src/cluster.jl")
include("../src/scale.jl")
include("../src/latexify_recipes.jl")
include("../src/printing.jl")
include("../src/indexing.jl")
include("../src/doubleSums.jl")
include("../src/averageSums.jl")
include("../src/indexedMeanfield.jl")
include("../src/indexedScale.jl")
include("../src/indexedCorrelation.jl")

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

@test H isa QAdd

for arg in H.arguments
    @test arg isa IndexedDoubleSum
end

H2 = H*a(l_ind)

@test Ssum1*a(l_ind) isa IndexedSingleSum

DSum1 = H.arguments[1]
Dsum2 = H.arguments[2]

@test isequal(DSum1*a(l_ind),Σ(Ssum1*a(l_ind),k_ind))
@test DSum1*a(l_ind) isa QAdd

@test isequal(Σ(Σ(σ(2,1,i_ind)*σ(1,2,j_ind),i_ind,[j_ind]),j_ind,[i_ind]),Σ(Σ(σ(2,1,i_ind)*σ(1,2,j_ind),i_ind,[j_ind]),j_ind))
@test isequal(Σ(Σ(σ(2,1,i_ind)*σ(1,2,j_ind),i_ind,[j_ind]),j_ind,[i_ind]),Σ(Σ(σ(2,1,i_ind),i_ind,[j_ind])*σ(1,2,j_ind),j_ind))

#for arg in H.arguments
#    @test arg isa IndexedDoubleSum
#end





end