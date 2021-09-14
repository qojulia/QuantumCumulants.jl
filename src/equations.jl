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
    jumps_dagger
    rates::Vector
    iv::SymbolicUtils.Sym
    varmap::Vector{Pair}
    order::Union{Int,Vector{<:Int},Nothing}
end

Base.getindex(de::AbstractMeanfieldEquations, i::Int) = de.equations[i]
Base.getindex(de::AbstractMeanfieldEquations, i) = de.equations[i]
Base.lastindex(de::AbstractMeanfieldEquations) = lastindex(de.equations)
Base.length(de::AbstractMeanfieldEquations) = length(de.equations)

function _append!(de::T, me::T) where T<:AbstractMeanfieldEquations
    append!(de.equations, me.equations)
    append!(de.operator_equations, me.operator_equations)
    append!(de.states, me.states)
    append!(de.operators, me.operators)
    append!(de.varmap, me.varmap)
    return de
end

# Substitution
function substitute(de::T,dict) where T<:AbstractMeanfieldEquations
    eqs = [substitute(eq, dict) for eq∈de.equations]
    states = getfield.(eqs, :lhs)
    fields = [getfield(de, s) for s∈fieldnames(T)[4:end]]
    return T(eqs, de.operator_equations, states, fields...)
end

# Simplification
function SymbolicUtils.simplify(de::T;kwargs...) where T<:AbstractMeanfieldEquations
    eqs = [SymbolicUtils.simplify(eq;kwargs...) for eq∈de.equations]
    eqs_op = [SymbolicUtils.simplify(eq;kwargs...) for eq∈de.operator_equations]
    fields = [getfield(de, s) for s∈fieldnames(T)[3:end]]
    return T(eqs,eqs_op,fields...)
end

# Adding MTK variables
function add_vars!(varmap, vs, t)
    keys = getindex.(varmap, 1)
    vals = getindex.(varmap, 2)
    hashkeys = map(hash, keys)
    hashvals = map(hash, vals)
    hashvs = map(hash, vs)
    for i=1:length(vs)
        if !(hashvs[i] ∈ hashkeys)
            var = make_var(vs[i], t)
            !(hash(var) ∈ hashvals) || @warn string("Two different averages have the exact same name. ",
                    "This may lead to unexpected behavior when trying to access the solution for $(vals[i])")
            push!(keys, vs[i])
            push!(vals, var)
            push!(hashkeys, hashvs[i])
        end
    end
    for i=length(varmap)+1:length(keys)
        push!(varmap, keys[i]=>vals[i])
    end
    return varmap
end

function make_var(v, t)
    sym = Symbol(string(v))
    var_f = SymbolicUtils.Sym{SymbolicUtils.FnType{Tuple{Any}, Complex}}(sym)
    return var_f(t)
end

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
    jumps_dagger
    rates::Vector
    iv::SymbolicUtils.Sym
    varmap::Vector{Pair}
    order::Union{Int,Vector{<:Int},Nothing}
    scale_aons
    names::Vector
    was_scaled::Vector{Bool}
end
