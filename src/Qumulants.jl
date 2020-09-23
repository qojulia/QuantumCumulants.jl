module Qumulants

import SymbolicUtils
using Combinatorics: partitions, combinations, permutations

export HilbertSpace, ProductSpace, ⊗,
        simplify_operators, substitute, expand,
        AbstractOperator, BasicOperator, OperatorTerm, embed,
        FockSpace, Destroy, Create,
        NLevelSpace, Transition,
        AbstractEquation, DifferentialEquation,
        heisenberg, commutator, acts_on,
        SymbolicNumber, NumberTerm, Parameter, @parameters, parameters,
                simplify_constants,
        Average, average, cumulant_expansion, get_order, cumulant,
        find_missing, complete, find_operators, fundamental_operators,
            unique_ops, get_symbolics, get_operators,
        build_ode, generate_ode,
        transition_superscript,
        IndexedDestroy, IndexedCreate, IndexedTransition, Index, Sum, find_index


"""
    SymbolicNumber <: Number

Abstract type for all symbolic numbers, i.e. [`Parameter`](@ref), [`Average`](@ref)
and corresponding expression trees.
"""
abstract type SymbolicNumber{T} <: Number end

"""
    NumberTerm <: SymbolicNumber

Expression tree consisting of [`SymbolicNumber`](@ref) variables.
"""
struct NumberTerm{T<:Number} <: SymbolicNumber{T}
    f::Function
    arguments::Vector
end

include("hilbertspace.jl")
include("operator.jl")
include("fock.jl")
include("nlevel.jl")
include("indexing.jl")
include("simplify.jl")
include("rules.jl")
include("equations.jl")
include("heisenberg.jl")
include("parameters.jl")
include("average.jl")
include("utils.jl")
include("diffeq.jl")
include("latexify_recipes.jl")
include("printing.jl")


end # module
