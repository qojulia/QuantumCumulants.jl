module QuantumCumulants

using SymbolicUtils: SymbolicUtils
import SymbolicUtils: substitute, BasicSymbolic, arguments

using Symbolics: Symbolics
using TermInterface: TermInterface

using SciMLBase: SciMLBase

using LinearAlgebra
using Combinatorics

using QuantumOpticsBase
import QuantumOpticsBase: ⊗, tensor

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
    find_missing,
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
    evaluate,
    value_map,
    NumberedOperator,
    change_index,
    order_by_index,
    split_sums,
    insert_index,
    eval_term

const NO_METADATA = SymbolicUtils.NO_METADATA

function source_metadata(source, name)
    Base.ImmutableDict{DataType,Any}(Symbolics.VariableSource, (source, name))
end

include("hilbertspace.jl")
include("qnumber.jl")
include("cnumber.jl")
include("fock.jl")
include("nlevel.jl")
include("spin.jl")
include("average.jl")
include("commutator.jl")
include("utils.jl")
include("cluster.jl")
include("latexify_recipes.jl")
include("printing.jl")
include("indexing.jl")
include("index_double_sums.jl")
include("index_average.jl")
include("index_utils.jl")
include("index_scale.jl")

end # module
