module QuantumCumulants

using Latexify: Latexify, @latexrecipe, latexify
using MacroTools: MacroTools
using LaTeXStrings: LaTeXStrings, latexstring, @L_str

using Combinatorics: partitions
using LinearAlgebra: LinearAlgebra
using OrderedCollections: OrderedCollections

using SciMLBase: SciMLBase
using QuantumOpticsBase: QuantumOpticsBase

using SymbolicUtils: SymbolicUtils, substitute, BasicSymbolic, operation, arguments, iscall

using Symbolics: Symbolics
using TermInterface: TermInterface

using ModelingToolkit: ModelingToolkit, complete, complete!
const MTK = ModelingToolkit

using SecondQuantizedAlgebra:
    SecondQuantizedAlgebra,
    QTerm,
    QSym,
    CNumber,
    Average,
    QNumber,
    QMul,
    QAdd,
    SpecialIndexedTerm,
    IndexedOperator,
    ClusterAon,
    CallableTransition,
    NumberedOperator,
    Parameter,
    ProductSpace,
    FockSpace,
    Destroy,
    Create,
    PhaseSpace,
    Position,
    Momentum,
    Index,
    SingleSum,
    DoubleSum,
    Transition,
    NLevelSpace,
    get_indices,
    commutator,
    _conj,
    _adjoint,
    get_i,
    hilbert,
    inorder!,
    levels,
    has_cluster,
    undo_average,
    _average,
    sym_average,
    IndexedVariable,
    DoubleIndexedVariable,
    SpecialIndexedAverage,
    _inconj,
    DoubleNumberedVariable,
    SingleNumberedVariable,
    create_index_arrays,
    _to_expression,
    AvgSums,
    get_numbers,
    @cnumbers,
    @qnumbers,
    @rnumbers,
    ClusterSpace,
    HilbertSpace,
    Pauli,
    PauliSpace,
    SpinSpace,
    RNumber,
    RealParameter,
    Spin,
    cnumber,
    cnumbers,
    find_operators,
    fundamental_operators,
    insert_index,
    order_by_index,
    reorder,
    rnumber,
    rnumbers,
    to_numeric,
    transition_superscript,
    unique_ops,
    unique_ops!,
    Σ,
    ∑,
    ⊗,
    tensor,
    acts_on,
    average,
    change_index,
    numeric_average,
    IndexedAverageSum,
    IndexedAverageDoubleSum

const SQA = SecondQuantizedAlgebra

export HilbertSpace,
    ProductSpace,
    ⊗,
    tensor,
    QSym,
    QTerm,
    @qnumbers,
    FockSpace,
    Destroy,
    Create,
    NLevelSpace,
    Transition,
    PauliSpace,
    Pauli,
    SpinSpace,
    Spin,
    PhaseSpace,
    Position,
    Momentum,
    MeanfieldEquations,
    meanfield,
    commutator,
    acts_on,
    CNumber,
    Parameter,
    @cnumbers,
    cnumbers,
    cnumber,
    RNumber,
    RealParameter,
    @rnumbers,
    rnumbers,
    rnumber,
    Average,
    average,
    cumulant_expansion,
    get_order,
    cumulant,
    find_missing,
    complete,
    complete!,
    find_operators,
    fundamental_operators,
    unique_ops,
    unique_ops!,
    to_numeric,
    numeric_average,
    initial_values,
    get_solution,
    get_scale_solution,
    CorrelationFunction,
    Spectrum,
    correlation_u0,
    correlation_p0,
    ClusterSpace,
    scale,
    transition_superscript,
    Index,
    reorder,
    IndexedOperator,
    SingleSum,
    IndexedVariable,
    DoubleIndexedVariable,
    DoubleSum,
    indexed_complete,
    IndexedCorrelationFunction,
    scale_term,
    scaleME,
    evalME,
    indexed_complete!,
    indexed_meanfield,
    subst_reds,
    AvgSums,
    plotME,
    IndexedAverageSum,
    IndexedAverageDoubleSum,
    SpecialIndexedTerm,
    find_missing_sums,
    Σ,
    ∑,
    evaluate,
    value_map,
    NumberedOperator,
    change_index,
    order_by_index,
    split_sums,
    insert_index,
    eval_term,
    MeanfieldNoiseEquations,
    IndexedMeanfieldNoiseEquations


include("equations.jl")
include("meanfield.jl")
include("utils.jl")
include("cumulant_expansion.jl")
include("diffeq.jl")
include("correlation.jl")
include("scale.jl")
include("measurement_backaction.jl")
include("measurement_backaction_indices.jl")
include("latexify_recipes.jl")
include("printing.jl")
include("index_average.jl")
include("index_meanfield.jl")
include("index_scale.jl")
include("index_correlation.jl")
include("index_utils.jl")


@deprecate heisenberg(args...; kwargs...) meanfield(args...; kwargs...)

end # module
