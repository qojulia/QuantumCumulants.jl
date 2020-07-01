module Qumulants

import SymbolicUtils
using Combinatorics: partitions, combinations, permutations

"""
    SymbolicNumber <: Number

Abstract type for all symbolic numbers, i.e. [`Parameter`](@ref), [`Average`](@ref)
and corresponding expression trees.
"""
abstract type SymbolicNumber <: Number end


export HilbertSpace, ProductSpace,
        simplify_operators, substitute, expand,
        AbstractOperator, BasicOperator, Identity, Zero, OperatorTerm, ⊗, embed,
        FockSpace, Destroy, Create,
        NLevelSpace, Transition, levels, ground_state,
        AbstractEquation, DifferentialEquation,
        heisenberg, commutator, acts_on,
        SymbolicNumber, NumberTerm, Parameter, @parameters, parameters,
                simplify_constants,
        Index, KroneckerDelta,
        Average, average, cumulant_expansion, get_order,
        find_missing, complete, find_operators, fundamental_operators,
            unique_ops, get_symbolics, get_operators,
        build_ode, generate_ode

include("ancrule.jl")
include("indexing.jl")
include("hilbertspace.jl")
include("operator.jl")
include("simplify.jl")
include("rules.jl")
include("fock.jl")
include("nlevel.jl")
include("equations.jl")
include("heisenberg.jl")
include("parameters.jl")
include("average.jl")
include("utils.jl")
include("diffeq.jl")
include("latexify_recipes.jl")
include("printing.jl")


end # module
