module CumulantHomotopy

using QuantumCumulants: QuantumCumulants, MeanfieldEquations, complete
using SecondQuantizedAlgebra: SecondQuantizedAlgebra
using Symbolics: Symbolics
using SymbolicUtils: SymbolicUtils
using HomotopyContinuation: HomotopyContinuation
using OrderedCollections: OrderedDict

const SQA = SecondQuantizedAlgebra
const HC = HomotopyContinuation

include("realification.jl")
include("homotopy.jl")

export StationaryPolynomialSystem, realify
export StationaryResult, stationary_states

"""
    stationary_state(args...; kwargs...)

Placeholder for the proposed physical stationary branch solver (ROADMAP Phase 2).
"""
function stationary_state(args...; kwargs...)
    error("`stationary_state` is not implemented yet; use `stationary_states` to enumerate all steady states (Phase 0).")
end

"""
    stationary_sequence(args...; kwargs...)

Placeholder for continuation across cumulant hierarchy orders (ROADMAP Phase 3).
"""
function stationary_sequence(args...; kwargs...)
    error("`stationary_sequence` is not implemented yet (ROADMAP Phase 3).")
end

export stationary_state, stationary_sequence

end
