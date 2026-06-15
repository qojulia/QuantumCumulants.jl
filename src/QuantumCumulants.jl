module QuantumCumulants

using Reexport: @reexport
@reexport using SecondQuantizedAlgebra

using SecondQuantizedAlgebra: SecondQuantizedAlgebra, QField, QAdd,
    average, commutator, undo_average, operators
import SecondQuantizedAlgebra as SQA
using SymbolicUtils: SymbolicUtils
using Symbolics: Symbolics, @variables
using Latexify: Latexify, latexify, @latexrecipe
using LaTeXStrings: latexstring
using ModelingToolkitBase: ModelingToolkitBase, complete, System
using OrderedCollections: OrderedCollections
using Combinatorics: Combinatorics, partitions
using TermInterface: TermInterface
using LinearAlgebra: I
const MTK = ModelingToolkitBase

export AbstractMeanfieldEquations, MeanfieldEquations, NoiseMeanfieldEquations
export EvolutionDirection, Forward, Backward
export meanfield, cumulant_expansion, cumulant, get_order
export complete, complete!, find_missing
export scale, scale!
export evaluate
export System, initial_values, get_solution, parameter_map
export CorrelationFunction, Spectrum, correlation_u0, correlation_p0
export translate_W_to_Y, modify_equations, modify_equations!
export simplify!

#  early types (abstract eqs supertype, treatment enum, direction tags), identity
include("equations.jl")
include("tree.jl")
include("canonical.jl")

#  operator algebra to moments
include("operator_drift.jl")
include("moments.jl")

#  the cumulant hierarchy (coupled moment equations) and its graph-backed wrappers
include("graph.jl")
include("equations_concrete.jl")
include("cumulant.jl")

# public orchestration (passes over the graph)
include("meanfield.jl")
include("completion.jl")
include("scaling.jl")
include("evaluate.jl")
include("mtk.jl")
include("correlation.jl")

# plain-text and LaTeX display (after all displayed types are defined)
include("printing.jl")

end # module
