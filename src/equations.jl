"""
    AbstractMeanfieldEquations

Supertype of [`MeanfieldEquations`](@ref) and [`NoiseMeanfieldEquations`](@ref).
"""
abstract type AbstractMeanfieldEquations end

"""
How one Hilbert-space factor (an SQA *subspace*, the tensor factor of the system's
`ProductSpace` identified by its `space_index`) is treated when building a moment's
canonical label:
- `Free`: its atom index stays a free symbolic index.
- `Scaled`: reduced under the permutation symmetry ``S_n`` of the identical atoms (the `scale` reduction).
- `Concrete`: pinned to fixed sites `1..M` (after `evaluate`).

The keying logic that consumes this lives in [canonical.jl](canonical.jl); it is defined
here only because the cumulant hierarchy stores a `Dict{Int, SubspaceTreatment}`.
"""
@enum SubspaceTreatment Free Scaled Concrete

"""
    EvolutionDirection

Supertype of the singleton tags [`Forward`](@ref) and [`Backward`](@ref) used to
dispatch the noise/retrodiction code path at compile time.
"""
abstract type EvolutionDirection end

"""Forward Heisenberg evolution (positive sign on `i[H, ·]`)."""
struct Forward <: EvolutionDirection end

"""Backward Heisenberg evolution (negative sign), used for retrodiction."""
struct Backward <: EvolutionDirection end

Base.length(eqs::AbstractMeanfieldEquations) = length(eqs.equations)
Base.getindex(eqs::AbstractMeanfieldEquations, i) = eqs.equations[i]
Base.lastindex(eqs::AbstractMeanfieldEquations) = lastindex(eqs.equations)
Base.iterate(eqs::AbstractMeanfieldEquations, st = 1) =
    st > length(eqs) ? nothing : (eqs.equations[st], st + 1)
Base.eltype(::Type{<:AbstractMeanfieldEquations}) = Symbolics.Equation
