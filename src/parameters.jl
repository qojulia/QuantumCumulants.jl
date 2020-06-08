"""
    SymbolicNumber <: Number

Abstract type for all symbolic numbers, i.e. [`Parameter`](@ref), [`Average`](@ref)
and corresponding expression trees.
"""
abstract type SymbolicNumber <: Number end

Base.:(==)(s::SymbolicNumber,x::Number) = false
Base.:(==)(x::Number,s::SymbolicNumber) = false
Base.:(==)(s1::SymbolicNumber,s2::SymbolicNumber) = false

"""
    NumberTerm <: SymbolicNumber

Expression tree consisting of [`SymbolicNumber`](@ref) variables.
"""
struct NumberTerm{T<:Number} <: SymbolicNumber
    f::Function
    arguments::Vector
end
function NumberTerm(f,args;type=promote_type(typeof.(args)...))
    return NumberTerm{type}(f,args)
end
Base.:(==)(t1::NumberTerm,t2::NumberTerm) = (t1.f===t2.f && t1.arguments==t2.arguments)

Base.zero(::Type{<:SymbolicNumber}) = 0
Base.one(::Type{<:SymbolicNumber}) = 1
Base.conj(t::NumberTerm) = NumberTerm(t.f, conj.(t.arguments))

"""
    Parameter <: SymbolicNumber

A parameter represented as a symbolic.
See also: [`parameters`](@ref), [`@parameters`](@ref)
"""
struct Parameter{T<:Number} <: SymbolicNumber
    name::Symbol
    index
end
Parameter{T}(name::Symbol) where T = Parameter{T}(name, default_index())
Parameter(name::Symbol) = Parameter{Number}(name)
Base.getindex(p::Parameter{T},index::Index) where T = Parameter{T}(p.name, index)
Base.getindex(p::Parameter{T},index::Index...) where T = Parameter{T}(p.name, index)
Base.:(==)(p::T, q::T) where T<:Parameter = (p.name==q.name && p.index==q.index)

get_index(p::Parameter) = p.index
function get_index(t::NumberTerm)
    idx = Index[]
    for arg in t.arguments
        append!(idx, get_index(arg))
    end
    return idx
end

# Methods
Base.conj(p::Parameter{<:Real}) = p

for f = [:+,:-,:*,:/,:^]
    @eval Base.$f(a::SymbolicNumber,b::Number) = NumberTerm($f, [a,b])
    @eval Base.$f(a::Number,b::SymbolicNumber) = NumberTerm($f, [a,b])
    @eval Base.$f(a::SymbolicNumber,b::SymbolicNumber) = NumberTerm($f, [a,b])
end
for f = [:cos,:sin,:tan,:sqrt,:conj]
    @eval Base.$f(a::SymbolicNumber) = NumberTerm($f, [a])
end

# Variadic methods
Base.:-(x::SymbolicNumber) = -1*x
for f in [:+,:*]
    @eval Base.$f(x::SymbolicNumber) = x
    @eval Base.$f(x::SymbolicNumber, w::SymbolicNumber...) = NumberTerm($f, [x;w...])
    @eval Base.$f(x::Number, y::SymbolicNumber, w::Number...) = NumberTerm($f, [x;y;w...])
    @eval Base.$f(x::SymbolicNumber, y::SymbolicNumber, w::Number...) = NumberTerm($f, [x;y;w...])
end

# Substitution
function substitute(t::NumberTerm, dict)
    if haskey(dict, t)
        return dict[t]
    else
        return NumberTerm(t.f, [substitute(arg, dict) for arg in t.arguments])
    end
end
substitute(x::SymbolicNumber, dict) = haskey(dict, x) ? dict[x] : x

# Conversion to SymbolicUtils
function parameter end
function average end
function _to_symbolic(p::Parameter{T}) where T
    SymbolicUtils.term(parameter, p.name, p.index; type=T)
end
_to_symbolic(n::NumberTerm{T}) where T = SymbolicUtils.Term{T}(n.f, _to_symbolic.(n.arguments))

# Convert back
function _to_qumulants(t::SymbolicUtils.Term{T}) where T<:Number
    if t.f===parameter
        return Parameter{T}(t.arguments...)
    elseif t.f===average
        return average(_to_qumulants(t.arguments[1]))
    else
        return NumberTerm(t.f, _to_qumulants.(t.arguments))
    end
end

"""
    @parameters(ps...)

Convenience macro to quickly define symbolic parameters.

Examples
========
```
julia> @parameters ω κ
(ω, κ)
```
"""
macro parameters(ps...)
    ex = Expr(:block)
    pnames = []
    for p in ps
        @assert p isa Symbol
        push!(pnames, p)
        ex_ = Expr(:(=), esc(p), Expr(:call, :Parameter, Expr(:quote, p)))
        push!(ex.args, ex_)
    end
    push!(ex.args, Expr(:tuple, map(esc, pnames)...))
    return ex
end

"""
    parameters(symbols::Symbol...)
    paramters(s::String)

Create symbolic parameters.

Expamples
=========
```
julia> ps = parameters(:a, :b)
(a, b)

julia> parameters("a b") == ps
true
```
"""
function parameters(syms::Symbol...)
    ps = Tuple(Parameter{Number}(s) for s in syms)
    return ps
end
function parameters(s::String)
    syms = [Symbol(p) for p in split(s, " ")]
    return parameters(syms...)
end

simplify_constants(s::SymbolicNumber;kwargs...) = s
function simplify_constants(t::NumberTerm;kwargs...)
    s = _to_symbolic(t)
    s_ = SymbolicUtils.simplify(s;kwargs...)
    return _to_qumulants(s_)
end
