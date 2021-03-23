"""
Abstract type for equations.
"""
abstract type AbstractHeisenbergEquation end

"""
    HeisenbergEquation <: AbstractHeisenbergEquation

Type defining a system of differential equations, where `lhs` is a vector of
derivatives and `rhs` is a vector of expressions. In addition, it keeps track
of the Hamiltonian, the collapse operators and the corresponding decay rates of
the system.

# Fields
*`lhs`: Vector of operators or averages of which the derivatives are taken.
*`rhs`: Vector of expressions to which the derivatives are equal.
*`hamiltonian`: Operator defining the system Hamiltonian.
*`jumps`: Vector of operators specifying the decay processes.
*`rates`: Decay rates corresponding to the `jumps`.

"""
struct HeisenbergEquation <: AbstractHeisenbergEquation
    equations::Vector{Symbolics.Equation}
    operator_equations::Vector{Symbolics.Equation}
    states::Vector
    operators::Vector{QNumber}
    hamiltonian::QNumber
    jumps::Vector{QNumber}
    rates::Vector
    iv::SymbolicUtils.Sym
    varmap::Vector{Pair}
    order::Union{Int,Vector{<:Int},Nothing}
end

Base.getindex(de::HeisenbergEquation, i::Int) = de.equations[i]
Base.getindex(de::HeisenbergEquation, i) = de.equations[i]
Base.lastindex(de::HeisenbergEquation) = lastindex(de.equations)
Base.length(de::HeisenbergEquation) = length(de.equations)

function _append!(de::HeisenbergEquation, he::HeisenbergEquation)
    append!(de.equations, he.equations)
    append!(de.operator_equations, he.operator_equations)
    append!(de.states, he.states)
    append!(de.operators, he.operators)
    append!(de.varmap, he.varmap)
    return de
end

# Substitution
function substitute(de::HeisenbergEquation,dict)
    eqs = [substitute(eq, dict) for eq∈de.equations]
    states = getfield.(eqs, :lhs)
    return HeisenbergEquation(eqs, de.operator_equations, states, de.operators, de.hamiltonian, de.jumps, de.rates, de.iv, de.varmap, de.order)
end

# Simplification
function SymbolicUtils.simplify(de::HeisenbergEquation;kwargs...)
    eqs = [SymbolicUtils.simplify(eq;kwargs...) for eq∈de.equations]
    eqs_op = [SymbolicUtils.simplify(eq;kwargs...) for eq∈de.operator_equations]
    return HeisenbergEquation(eqs,eqs_op,de.states,de.operators,de.hamiltonian,de.jumps,de.rates,de.iv,de.varmap,de.order)
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
    # Vars need to be of symtype Real; otherwise they get converted internally
    # here: https://github.com/JuliaSymbolics/Symbolics.jl/blob/1f450def014628f7fe5321282f33972527fb0318/src/utils.jl#L96
    # which causes the missingvars check to fail due to a setdiff on lists with different types
    var_f = SymbolicUtils.Sym{SymbolicUtils.FnType{Tuple{Any}, Real}}(sym)
    return var_f(t)
end

function make_varmap(vs, t)
    varmap = Pair{Any,Any}[]
    add_vars!(varmap, vs, t)
    return varmap
end
