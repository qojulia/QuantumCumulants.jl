"""
Abstract type for equations.
"""
abstract type AbstractMeanfieldEquations end

"""
    MeanfieldEquations <: AbstractMeanfieldEquations

Type defining a system of differential equations, where `lhs` is a vector of
derivatives and `rhs` is a vector of expressions. In addition, it keeps track
of the Hamiltonian, the collapse operators and the corresponding decay rates of
the system.

# Fields
*`equations`: Vector of the differential equations of averages.
*`operator_equations`: Vector of the operator differential equations.
*`states`: Vector containing the averages on the left-hand-side of the equations.
*`operators`: Vector containing the operators on the left-hand-side of the equations.
*`hamiltonian`: Operator defining the system Hamiltonian.
*`jumps`: Vector of operators specifying the decay processes.
*`jumps`: Vector of operators specifying the adjoint of the decay processes.
*`rates`: Decay rates corresponding to the `jumps`.
*`iv`: The independent variable (time parameter) of the system.
*`varmap`: Vector of pairs that map the averages to time-dependent variables.
    That format is necessary for ModelingToolkit functionality.
*`order`: The order at which the [`cumulant_expansion`](@ref) has been performed.

"""
struct MeanfieldEquations <: AbstractMeanfieldEquations
    equations::Vector{Symbolics.Equation}
    operator_equations::Vector{Symbolics.Equation}
    states::Vector
    operators::Vector{QNumber}
    hamiltonian::QNumber
    jumps::Vector
    jumps_dagger::Any
    rates::Vector
    iv::MTK.Num
    varmap::Vector{Pair}
    order::Union{Int,Vector{<:Int},Nothing}
end

"""
    IndexedMeanfieldEquations <: AbstractMeanfieldEquations

Type defining a system of differential equations, where `lhs` is a vector of
derivatives and `rhs` is a vector of expressions. In addition, it keeps track
of the Hamiltonian, the collapse operators and the corresponding decay rates of
the system. Similar to [`MeanfieldEquations`](@ref), specialized for equations,
that are using [`Index`](@ref) entities.

# Fields
*`equations`: Vector of the differential equations of averages.
*`operator_equations`: Vector of the operator differential equations.
*`states`: Vector containing the averages on the left-hand-side of the equations.
*`operators`: Vector containing the operators on the left-hand-side of the equations.
*`hamiltonian`: Operator defining the system Hamiltonian.
*`jumps`: Vector of operators specifying the decay processes.
*`jumps_dagger`: Vector of operators specifying the adjoint of the decay processes.
*`rates`: Decay rates corresponding to the `jumps`.
*`iv`: The independent variable (time parameter) of the system.
*`varmap`: Vector of pairs that map the averages to time-dependent variables.
    That format is necessary for ModelingToolkit functionality.
*`order`: The order at which the [`cumulant_expansion`](@ref) has been performed.

"""
struct IndexedMeanfieldEquations <: AbstractMeanfieldEquations #these are for easier dispatching of meanfield, complete,... functions
    equations::Vector{Symbolics.Equation}
    operator_equations::Vector{Symbolics.Equation}
    states::Vector
    operators::Vector{QNumber}
    hamiltonian::QNumber
    jumps::Vector
    jumps_dagger::Any
    rates::Vector
    iv::MTK.Num
    varmap::Vector{Pair}
    order::Union{Int,Vector{<:Int},Nothing}
end

"""
    EvaledMeanfieldEquations <: AbstractMeanfieldEquations

Type defining a system of differential equations, where `lhs` is a vector of
derivatives and `rhs` is a vector of expressions. In addition, it keeps track
of the Hamiltonian, the collapse operators and the corresponding decay rates of
the system. Similar to [`MeanfieldEquations`](@ref), specialized for equations,
that are evaluated using the [`evaluate`](@ref) function.

# Fields
*`equations`: Vector of the differential equations of averages.
*`operator_equations`: Vector of the operator differential equations.
*`states`: Vector containing the averages on the left-hand-side of the equations.
*`operators`: Vector containing the operators on the left-hand-side of the equations.
*`hamiltonian`: Operator defining the system Hamiltonian.
*`jumps`: Vector of operators specifying the decay processes.
*`jumps_dagger`: Vector of operators specifying the adjoint of the decay processes.
*`rates`: Decay rates corresponding to the `jumps`.
*`iv`: The independent variable (time parameter) of the system.
*`varmap`: Vector of pairs that map the averages to time-dependent variables.
    That format is necessary for ModelingToolkit functionality.
*`order`: The order at which the [`cumulant_expansion`](@ref) has been performed.

"""
struct EvaledMeanfieldEquations <: AbstractMeanfieldEquations
    equations::Vector{Symbolics.Equation}
    operator_equations::Vector{Symbolics.Equation}
    states::Vector
    operators::Vector{QNumber}
    hamiltonian::QNumber
    jumps::Vector
    jumps_dagger::Any
    rates::Vector
    iv::MTK.Num
    varmap::Vector{Pair}
    order::Union{Int,Vector{<:Int},Nothing}
end

Base.getindex(de::AbstractMeanfieldEquations, i::Int) = de.equations[i]
Base.getindex(de::AbstractMeanfieldEquations, i) = de.equations[i]
Base.lastindex(de::AbstractMeanfieldEquations) = lastindex(de.equations)
Base.length(de::AbstractMeanfieldEquations) = length(de.equations)

function _append!(de::T, me::T) where {T<:AbstractMeanfieldEquations}
    append!(de.equations, me.equations)
    append!(de.operator_equations, me.operator_equations)
    append!(de.states, me.states)
    append!(de.operators, me.operators)
    append!(de.varmap, me.varmap)
    return de
end

# Substitution
function SymbolicUtils.substitute(de::T, dict) where {T<:AbstractMeanfieldEquations}
    eqs = [substitute(eq, dict) for eq ∈ de.equations]
    states = getfield.(eqs, :lhs)
    fields = [getfield(de, s) for s ∈ fieldnames(T)[4:end]]
    return T(eqs, de.operator_equations, states, fields...)
end

function SymbolicUtils.substitute(de::IndexedMeanfieldEquations, dict)
    eqs = [
        Symbolics.Equation(
            inorder!(substitute(eq.lhs, dict)),
            inorder!(substitute(eq.rhs, dict)),
        ) for eq ∈ de.equations
    ]
    states = getfield.(eqs, :lhs)
    fields = [getfield(de, s) for s ∈ fieldnames(IndexedMeanfieldEquations)[4:end]]
    return IndexedMeanfieldEquations(eqs, de.operator_equations, states, fields...)
end

# Simplification
function SymbolicUtils.simplify(de::T; kwargs...) where {T<:AbstractMeanfieldEquations}
    eqs = [SymbolicUtils.simplify(eq; kwargs...) for eq ∈ de.equations]
    eqs_op = [SymbolicUtils.simplify(eq; kwargs...) for eq ∈ de.operator_equations]
    fields = [getfield(de, s) for s ∈ fieldnames(T)[3:end]]
    return T(eqs, eqs_op, fields...)
end

# Adding MTK variables
function add_vars!(varmap, vs, t)
    keys = getindex.(varmap, 1)
    vals = getindex.(varmap, 2)
    hashkeys = map(hash, keys)
    hashvals = map(hash, vals)
    hashvs = map(hash, vs)
    for i = 1:length(vs)
        if !(hashvs[i] ∈ hashkeys)
            var = make_var(vs[i], t)
            !(hash(var) ∈ hashvals) || @warn string(
                "Two different averages have the exact same name. ",
                "This may lead to unexpected behavior when trying to access the solution for $(vals[i])",
            )
            push!(keys, vs[i])
            push!(vals, var)
            push!(hashkeys, hashvs[i])
        end
    end
    for i = (length(varmap)+1):length(keys)
        push!(varmap, keys[i]=>vals[i])
    end
    return varmap
end

function make_var(v, t)
    sym = Symbol(string(v))
    d = source_metadata(:make_var, sym)
    var_f = SymbolicUtils.Sym{SymbolicUtils.FnType{Tuple{Any},Complex}}(sym; metadata = d)
    return SymbolicUtils.Term{Complex}(var_f, [t]; metadata = d)
end

source_metadata(source, name) =
    Base.ImmutableDict{DataType,Any}(Symbolics.VariableSource, (source, name))

function make_varmap(vs, t)
    varmap = Pair{Any,Any}[]
    add_vars!(varmap, vs, t)
    return varmap
end

struct ScaledMeanfieldEquations <: AbstractMeanfieldEquations
    equations::Vector{Symbolics.Equation}
    operator_equations::Vector{Symbolics.Equation}
    states::Vector
    operators::Vector{QNumber}
    hamiltonian::QNumber
    jumps::Vector
    jumps_dagger::Any
    rates::Vector
    iv::MTK.Num
    varmap::Vector{Pair}
    order::Union{Int,Vector{<:Int},Nothing}
    scale_aons::Any
    names::Vector
    was_scaled::Vector{Bool}
end


# noise equations
"""
    MeanfieldNoiseEquations

Mean field equations including a separate set of equations describing the
noise generated by measurement backactions.
"""
struct MeanfieldNoiseEquations <: AbstractMeanfieldEquations
    equations::Vector{Symbolics.Equation}
    operator_equations::Vector{Symbolics.Equation}
    noise_equations::Vector{Symbolics.Equation}
    operator_noise_equations::Vector{Symbolics.Equation} # useless but needed to create eqs
    states::Vector
    operators::Vector{QNumber}
    hamiltonian::QNumber
    jumps::Vector
    jumps_dagger::Any
    rates::Vector
    efficiencies::Vector
    iv::MTK.Num
    varmap::Vector{Pair}
    order::Union{Int,Vector{<:Int},Nothing}
end

"""
    IndexedMeanfieldNoiseEquations

Like a [`MeanfieldNoiseEquations`](@ref), but with symbolic indices.
"""
struct IndexedMeanfieldNoiseEquations <: AbstractMeanfieldEquations
    equations::Vector{Symbolics.Equation}
    operator_equations::Vector{Symbolics.Equation}
    noise_equations::Vector{Symbolics.Equation}
    operator_noise_equations::Vector{Symbolics.Equation}
    states::Vector
    operators::Vector{QNumber}
    hamiltonian::QNumber
    jumps::Vector
    jumps_dagger::Any
    rates::Vector
    efficiencies::Vector
    iv::MTK.Num
    varmap::Vector{Pair}
    order::Union{Int,Vector{<:Int},Nothing}
end

# TODO: ScaledMeanfieldNoiseEquations

"""
    BackwardMeanfieldNoiseEquations

Mean field equations for the backward propagation with 
noise generated by measurement backactions. 

See also: [`MeanfieldNoiseEquations`](@ref)
"""
struct BackwardMeanfieldNoiseEquations <: AbstractMeanfieldEquations
    equations::Vector{Symbolics.Equation}
    operator_equations::Vector{Symbolics.Equation}
    noise_equations::Vector{Symbolics.Equation}
    operator_noise_equations::Vector{Symbolics.Equation} # useless but needed to create eqs # TODO: delete?
    states::Vector
    operators::Vector{QNumber}
    hamiltonian::QNumber
    jumps::Vector
    jumps_dagger::Any
    rates::Vector
    efficiencies::Vector
    iv::MTK.Num
    varmap::Vector{Pair}
    order::Union{Int,Vector{<:Int},Nothing}
end

const NoiseEquations = Union{MeanfieldNoiseEquations,BackwardMeanfieldNoiseEquations}

# TODO: ScaledBackwardMeanfieldNoiseEquations, EvaledBackwardMeanfieldEquations