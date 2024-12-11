module QuantumCumulants

import SymbolicUtils
import SymbolicUtils: substitute, BasicSymbolic

import Symbolics
import TermInterface

import SciMLBase

import ModelingToolkit
import ModelingToolkit: complete, complete!
const MTK = ModelingToolkit

using Combinatorics: partitions, combinations, levicivita
using LinearAlgebra

using QuantumOpticsBase
import QuantumOpticsBase: ⊗, tensor

export HilbertSpace, ProductSpace, ⊗, tensor,
        QSym, QTerm, @qnumbers,
        FockSpace, Destroy, Create, GroundStateProjection,
        NLevelSpace, Transition,
        PauliSpace, Pauli, SpinSpace, Spin,
        MeanfieldEquations,
        meanfield, commutator, acts_on,
        CNumber, Parameter, @cnumbers, cnumbers, cnumber, RNumber, RealParameter, @rnumbers, rnumbers, rnumber,
        Average, average, cumulant_expansion, get_order, cumulant,
        find_missing, complete, complete!, find_operators, fundamental_operators,
            unique_ops, unique_ops!, to_numeric, numeric_average, initial_values, get_solution, get_scale_solution,
        CorrelationFunction, Spectrum, correlation_u0, correlation_p0,
        ClusterSpace,
        scale,
        transition_superscript,
        Index, reorder, IndexedOperator, SingleSum, IndexedVariable, DoubleIndexedVariable,
        DoubleSum, indexed_complete, IndexedCorrelationFunction, scale_term,
        scaleME, evalME, indexed_complete!, indexed_meanfield, subst_reds, AvgSums, plotME,
        IndexedAverageSum, IndexedAverageDoubleSum, SpecialIndexedTerm, find_missing_sums, Σ, ∑,
        evaluate, value_map, NumberedOperator, change_index, order_by_index, split_sums, insert_index, eval_term,
        MeanfieldNoiseEquations, 
        IndexedMeanfieldNoiseEquations#, indexed_arithmetic, indexed_noise, simplified_indexed_complete!


const NO_METADATA = SymbolicUtils.NO_METADATA

source_metadata(source, name) =
    Base.ImmutableDict{DataType, Any}(Symbolics.VariableSource, (source, name))

include("hilbertspace.jl")
include("qnumber.jl")
include("cnumber.jl")
include("fock.jl")
include("nlevel.jl")
include("spin.jl")
include("equations.jl")
include("meanfield.jl")
include("average.jl")
include("utils.jl")
include("diffeq.jl")
include("correlation.jl")
include("cluster.jl")
include("scale.jl")
include("measurement_backaction.jl")
include("measurement_backaction_indices.jl")
include("latexify_recipes.jl")
include("printing.jl")
include("indexing.jl")
include("index_double_sums.jl")
include("index_average.jl")
include("index_meanfield.jl")
include("index_scale.jl")
include("index_correlation.jl")
include("index_utils.jl")


@deprecate heisenberg(args...; kwargs...) meanfield(args...; kwargs...)

end # module
