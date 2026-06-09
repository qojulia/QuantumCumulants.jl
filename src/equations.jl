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
here only because `MeanfieldEquations` stores a `Dict{Int, SubspaceTreatment}` field.
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

"""
    MeanfieldEquations

Concrete equation set produced by [`meanfield`](@ref) when no measurement backaction is
requested. All type parameters are bound concretely.

# Fields
* `equations`: the averaged differential equations (left-hand side average, right-hand
  side drift).
* `operator_equations`: the same equations at the operator level.
* `states`: the averages on the left-hand sides.
* `operators`: the operators on the left-hand sides.
* `hamiltonian`: the system Hamiltonian.
* `jumps`, `jumps_dagger`: the collapse operators and their adjoints.
* `rates`: the decay rates corresponding to `jumps`.
* `iv`: the independent (time) variable.
* `order`: the cumulant-expansion order, or `nothing`.
* `direction`: [`Forward`](@ref) or [`Backward`](@ref) evolution.
"""
struct MeanfieldEquations{
        O <: Union{Nothing, Vector{Int}},
        H <: QField,
        Op <: QField,
        Jt,
        Jdt,
        R,
        S <: SymbolicUtils.BasicSymbolic,
        D <: EvolutionDirection,
    } <: AbstractMeanfieldEquations
    equations::Vector{Symbolics.Equation}
    operator_equations::Vector{Symbolics.Equation}
    states::Vector{S}
    operators::Vector{Op}
    hamiltonian::H
    jumps::Vector{Jt}
    jumps_dagger::Vector{Jdt}
    rates::Vector{R}
    iv::Symbolics.Num
    order::O
    direction::D
    treatments::Dict{Int, SubspaceTreatment}

    function MeanfieldEquations(
            equations::Vector{Symbolics.Equation},
            operator_equations::Vector{Symbolics.Equation},
            states::Vector{S},
            operators::Vector{Op},
            hamiltonian::H,
            jumps::Vector{Jt},
            jumps_dagger::Vector{Jdt},
            rates::Vector{R},
            iv::Symbolics.Num,
            order::O,
            direction::D = Forward();
            treatments::Dict{Int, SubspaceTreatment} = Dict{Int, SubspaceTreatment}(),
        ) where {
            O, H <: QField, Op <: QField, Jt, Jdt, R,
            S <: SymbolicUtils.BasicSymbolic, D <: EvolutionDirection,
        }
        n = length(equations)
        @assert n == length(operator_equations) == length(states) == length(operators) (
            "equations/states/operators must have matching lengths"
        )
        @assert length(jumps) == length(jumps_dagger) == length(rates) (
            "jumps/jumps_dagger/rates must have matching lengths"
        )
        return new{O, H, Op, Jt, Jdt, R, S, D}(
            equations, operator_equations, states, operators,
            hamiltonian, jumps, jumps_dagger, rates, iv, order, direction, treatments,
        )
    end
end

"""
    NoiseMeanfieldEquations

Concrete equation set produced by [`meanfield`](@ref) when `efficiencies` is supplied,
adding a measurement-backaction noise drift to the deterministic equations. The
`direction` tag ([`Forward`](@ref) or [`Backward`](@ref)) drives compile-time dispatch of
the noise assembly.

# Fields
* `equations`: the deterministic averaged differential equations.
* `noise_equations`: the averaged noise (Brownian) drift per state.
* `operator_equations`, `operator_noise_equations`: the operator-level counterparts.
* `states`: the averages on the left-hand sides.
* `operators`: the operators on the left-hand sides.
* `hamiltonian`: the system Hamiltonian.
* `jumps`, `jumps_dagger`: the collapse operators and their adjoints.
* `rates`: the decay rates corresponding to `jumps`.
* `efficiencies`: the detector efficiencies per jump.
* `iv`: the independent (time) variable.
* `order`: the cumulant-expansion order, or `nothing`.
* `direction`: [`Forward`](@ref) or [`Backward`](@ref) evolution.
"""
struct NoiseMeanfieldEquations{
        O <: Union{Nothing, Vector{Int}},
        H <: QField,
        Op <: QField,
        Jt,
        Jdt,
        R, E,
        S <: SymbolicUtils.BasicSymbolic,
        D <: EvolutionDirection,
    } <: AbstractMeanfieldEquations
    equations::Vector{Symbolics.Equation}
    noise_equations::Vector{Symbolics.Equation}
    operator_equations::Vector{Symbolics.Equation}
    operator_noise_equations::Vector{Symbolics.Equation}
    states::Vector{S}
    operators::Vector{Op}
    hamiltonian::H
    jumps::Vector{Jt}
    jumps_dagger::Vector{Jdt}
    rates::Vector{R}
    efficiencies::Vector{E}
    iv::Symbolics.Num
    order::O
    direction::D
    # Per-subspace treatment state; see MeanfieldEquations.
    treatments::Dict{Int, SubspaceTreatment}

    function NoiseMeanfieldEquations(
            equations::Vector{Symbolics.Equation},
            noise_equations::Vector{Symbolics.Equation},
            operator_equations::Vector{Symbolics.Equation},
            operator_noise_equations::Vector{Symbolics.Equation},
            states::Vector{S},
            operators::Vector{Op},
            hamiltonian::H,
            jumps::Vector{Jt},
            jumps_dagger::Vector{Jdt},
            rates::Vector{R},
            efficiencies::Vector{E},
            iv::Symbolics.Num,
            order::O,
            direction::D;
            treatments::Dict{Int, SubspaceTreatment} = Dict{Int, SubspaceTreatment}(),
        ) where {O, H, Op, Jt, Jdt, R, E, S, D}
        n = length(equations)
        @assert n == length(noise_equations) == length(operator_equations) ==
            length(operator_noise_equations) == length(states) == length(operators) (
            "equations/states/operators/noise must have matching lengths"
        )
        @assert length(jumps) == length(jumps_dagger) == length(rates) == length(efficiencies) (
            "jumps/jumps_dagger/rates/efficiencies must have matching lengths"
        )
        return new{O, H, Op, Jt, Jdt, R, E, S, D}(
            equations, noise_equations, operator_equations,
            operator_noise_equations, states, operators,
            hamiltonian, jumps, jumps_dagger, rates,
            efficiencies, iv, order, direction, treatments,
        )
    end
end

"""
    MeanfieldEquations(eqs::NoiseMeanfieldEquations)

Strip the measurement-backaction columns (`noise_equations`,
`operator_noise_equations`, `efficiencies`) and return a plain
[`MeanfieldEquations`](@ref) carrying only the deterministic drift. The `direction`
tag is preserved.

Use this to compare a stochastic simulation against the no-measurement
evolution: build the noisy system once, then pass `MeanfieldEquations(eqs)`
to `System` / `ODEProblem` for the deterministic reference.
"""
function MeanfieldEquations(eqs::NoiseMeanfieldEquations)
    return MeanfieldEquations(
        eqs.equations, eqs.operator_equations, eqs.states, eqs.operators,
        eqs.hamiltonian, eqs.jumps, eqs.jumps_dagger, eqs.rates,
        eqs.iv, eqs.order, eqs.direction;
        treatments = eqs.treatments,
    )
end

Base.length(eqs::AbstractMeanfieldEquations) = length(eqs.equations)
Base.getindex(eqs::AbstractMeanfieldEquations, i) = eqs.equations[i]
Base.lastindex(eqs::AbstractMeanfieldEquations) = lastindex(eqs.equations)
Base.iterate(eqs::AbstractMeanfieldEquations, st = 1) =
    st > length(eqs) ? nothing : (eqs.equations[st], st + 1)
Base.eltype(::Type{<:AbstractMeanfieldEquations}) = Symbolics.Equation
