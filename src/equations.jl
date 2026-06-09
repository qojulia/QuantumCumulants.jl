"""
    AbstractMeanFieldEquations

Supertype of [`MeanFieldEquations`](@ref) and [`NoiseMeanFieldEquations`](@ref).
"""
abstract type AbstractMeanFieldEquations end

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
    MeanFieldEquations

Concrete equation set produced by [`meanfield`](@ref) when no measurement backaction is
requested. All type parameters are bound concretely.

# Fields
* `equations`: the averaged differential equations (left-hand side average, right-hand
  side drift).
* `operator_equations`: the same equations at the operator level.
* `states`: the averages on the left-hand sides.
* `operators`: the operators on the left-hand sides.
* `initial_operators`: the operators originally requested, before completion added any.
* `hamiltonian`: the system Hamiltonian.
* `jumps`, `jumps_dagger`: the collapse operators and their adjoints.
* `rates`: the decay rates corresponding to `jumps`.
* `iv`: the independent (time) variable.
* `order`: the cumulant-expansion order, or `nothing`.
* `direction`: [`Forward`](@ref) or [`Backward`](@ref) evolution.
"""
struct MeanFieldEquations{
        O <: Union{Nothing, Vector{Int}},
        H <: QField,
        Op <: QField,
        Jt,
        Jdt,
        R,
        S <: SymbolicUtils.BasicSymbolic,
        D <: EvolutionDirection,
    } <: AbstractMeanFieldEquations
    equations::Vector{Symbolics.Equation}
    operator_equations::Vector{Symbolics.Equation}
    states::Vector{S}
    operators::Vector{Op}
    initial_operators::Vector{Op}
    hamiltonian::H
    jumps::Vector{Jt}
    jumps_dagger::Vector{Jdt}
    rates::Vector{R}
    iv::Symbolics.Num
    order::O
    direction::D
    # Per-subspace treatment state, keyed by space index and valued by the `Int`
    # of the `SubspaceTreatment` enum (defined in canonical.jl, included after this file).
    # An empty map means every subspace is Free.
    treatments::Dict{Int, Int}

    function MeanFieldEquations(
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
            initial_operators::Vector{Op} = copy(operators),
            treatments::Dict{Int, Int} = Dict{Int, Int}(),
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
            equations, operator_equations, states, operators, initial_operators,
            hamiltonian, jumps, jumps_dagger, rates, iv, order, direction, treatments,
        )
    end
end

"""
    NoiseMeanFieldEquations

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
* `initial_operators`: the operators originally requested, before completion added any.
* `hamiltonian`: the system Hamiltonian.
* `jumps`, `jumps_dagger`: the collapse operators and their adjoints.
* `rates`: the decay rates corresponding to `jumps`.
* `efficiencies`: the detector efficiencies per jump.
* `iv`: the independent (time) variable.
* `order`: the cumulant-expansion order, or `nothing`.
* `direction`: [`Forward`](@ref) or [`Backward`](@ref) evolution.
"""
struct NoiseMeanFieldEquations{
        O <: Union{Nothing, Vector{Int}},
        H <: QField,
        Op <: QField,
        Jt,
        Jdt,
        R, E,
        S <: SymbolicUtils.BasicSymbolic,
        D <: EvolutionDirection,
    } <: AbstractMeanFieldEquations
    equations::Vector{Symbolics.Equation}
    noise_equations::Vector{Symbolics.Equation}
    operator_equations::Vector{Symbolics.Equation}
    operator_noise_equations::Vector{Symbolics.Equation}
    states::Vector{S}
    operators::Vector{Op}
    initial_operators::Vector{Op}
    hamiltonian::H
    jumps::Vector{Jt}
    jumps_dagger::Vector{Jdt}
    rates::Vector{R}
    efficiencies::Vector{E}
    iv::Symbolics.Num
    order::O
    direction::D
    # Per-subspace treatment state; see MeanFieldEquations.
    treatments::Dict{Int, Int}

    function NoiseMeanFieldEquations(
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
            initial_operators::Vector{Op} = copy(operators),
            treatments::Dict{Int, Int} = Dict{Int, Int}(),
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
            operator_noise_equations, states, operators, initial_operators,
            hamiltonian, jumps, jumps_dagger, rates,
            efficiencies, iv, order, direction, treatments,
        )
    end
end

"""
    MeanFieldEquations(eqs::NoiseMeanFieldEquations)

Strip the measurement-backaction columns (`noise_equations`,
`operator_noise_equations`, `efficiencies`) and return a plain
[`MeanFieldEquations`](@ref) carrying only the deterministic drift. The `direction`
tag is preserved.

Use this to compare a stochastic simulation against the no-measurement
evolution: build the noisy system once, then pass `MeanFieldEquations(eqs)`
to `System` / `ODEProblem` for the deterministic reference.
"""
function MeanFieldEquations(eqs::NoiseMeanFieldEquations)
    return MeanFieldEquations(
        eqs.equations, eqs.operator_equations, eqs.states, eqs.operators,
        eqs.hamiltonian, eqs.jumps, eqs.jumps_dagger, eqs.rates,
        eqs.iv, eqs.order, eqs.direction;
        initial_operators = copy(eqs.initial_operators),
        treatments = copy(eqs.treatments),
    )
end

Base.length(eqs::AbstractMeanFieldEquations) = length(eqs.equations)
Base.getindex(eqs::AbstractMeanFieldEquations, i) = eqs.equations[i]
Base.lastindex(eqs::AbstractMeanFieldEquations) = lastindex(eqs.equations)
Base.iterate(eqs::AbstractMeanFieldEquations, st = 1) =
    st > length(eqs) ? nothing : (eqs.equations[st], st + 1)
Base.eltype(::Type{<:AbstractMeanFieldEquations}) = Symbolics.Equation
