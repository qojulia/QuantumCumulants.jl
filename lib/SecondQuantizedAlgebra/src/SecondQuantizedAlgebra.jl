module SecondQuantizedAlgebra

using SymbolicUtils: SymbolicUtils
import SymbolicUtils: substitute, BasicSymbolic, arguments, iscall, operation

using Symbolics: Symbolics
using TermInterface: TermInterface

using Combinatorics: combinations, levicivita

import SciMLBase
using SciMLBase: SciMLBase
using QuantumOpticsBase: QuantumOpticsBase
import QuantumOpticsBase: ⊗, tensor

using Latexify
import MacroTools
using LaTeXStrings

const NO_METADATA = SymbolicUtils.NO_METADATA

source_metadata(source, name) =
    Base.ImmutableDict{DataType, Any}(Symbolics.VariableSource, (source, name))

import ModelingToolkit as MTK

include("hilbertspace.jl")
include("qnumber.jl")
include("cnumber.jl")
include("fock.jl")
include("nlevel.jl")
include("spin.jl")
include("commutator.jl")

include("average.jl")
include("utils.jl")
include("cluster.jl")

include("latexify_recipes.jl")
include("printing.jl")

include("indexing.jl")
include("index_numbered_operator.jl")
include("index_double_sums.jl")
include("index_average.jl")
include("index_utils.jl")

export HilbertSpace, ProductSpace, ⊗, tensor,
    QSym, QTerm, @qnumbers,
    FockSpace, Destroy, Create,
    NLevelSpace, Transition,
    PauliSpace, Pauli, SpinSpace, Spin,
    commutator, acts_on,
    CNumber, Parameter, @cnumbers, cnumbers, cnumber, RNumber, RealParameter, @rnumbers,
    rnumbers, rnumber,
    unique_ops, unique_ops!, to_numeric, numeric_average,
    ClusterSpace, find_operators, fundamental_operators,
    transition_superscript, Average, average,
    Index, reorder, IndexedOperator, SingleSum, IndexedVariable, DoubleIndexedVariable,
    DoubleSum, SpecialIndexedTerm, Σ, ∑,
    NumberedOperator, change_index, order_by_index, insert_index

end
