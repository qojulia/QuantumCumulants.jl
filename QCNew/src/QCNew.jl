module QCNew

using Reexport: @reexport
@reexport using SecondQuantizedAlgebra

using SecondQuantizedAlgebra: SecondQuantizedAlgebra, QField, QAdd,
    average, commutator, undo_average, operators
import SecondQuantizedAlgebra as SQA
using SymbolicUtils: SymbolicUtils
using Symbolics: Symbolics, @variables
using ModelingToolkitBase: ModelingToolkitBase, complete, System
using OrderedCollections: OrderedCollections
using Combinatorics: Combinatorics
using TermInterface: TermInterface
using LinearAlgebra: I
const MTK = ModelingToolkitBase

export AbstractMeanFieldEquations, MeanFieldEquations, NoiseMeanFieldEquations
export EvolutionDirection, Forward, Backward
export meanfield, cumulant_expansion, cumulant, get_order
export complete, complete!, find_missing
export scale, scale!
export evaluate
export System, initial_values, get_solution, parameter_map
export CorrelationFunction, Spectrum, correlation_u0, correlation_p0
export translate_W_to_Y, modify_equations, modify_equations!
export simplify!

# The clean rebuild grows part by part. The tested moment-graph KERNEL is ported
# in first (Layers 1 to 3); public orchestration (meanfield, complete, scale,
# evaluate, correlation, mtk) is built on top of it incrementally.

# Layer 1: primitives, equation containers, identity
include("channels.jl")        # jump classification, index/average predicates
include("equations.jl")       # MeanFieldEquations, NoiseMeanFieldEquations, directions
include("tree.jl")            # Tree traversal (mapleaves/foldleaves/eachleaf)
include("canonical.jl")       # CanonCtx, build_ctx, canon_key, orbit_key

# Layer 2: operator algebra to moments
include("operator_drift.jl")  # coherent + Lindblad + backward recycling + noise dispatch
include("cumulant.jl")        # cumulant_expansion, cumulant, get_order
include("noise.jl")           # measurement-backaction noise drift builders
include("moments.jl")         # average_and_truncate, NodeData, derive

# Layer 3: the moment dependency graph IR
include("graph.jl")           # SystemSpec, MomentGraph, seed, closure!, frontier, lower_to_eqs

# Layer 5: public orchestration (thin passes over the graph)
include("meanfield.jl")       # meanfield (seed + lower; det and noise unified)
include("completion.jl")      # complete, complete!, find_missing (closure!/frontier)
include("scaling.jl")         # scale, scale! (quotient! pass over orbit_key)
include("evaluate.jl")        # evaluate (specialize pass: unroll indexed sums)
include("mtk.jl")             # System / initial_values / get_solution / parameter_map
include("correlation.jl")     # CorrelationFunction, Spectrum, System(c), correlation_u0/p0
include("measurement_backaction.jl")  # translate_W_to_Y (dW -> dY record)
include("modify.jl")          # modify_equations(!) RHS post-processing hook

end # module
