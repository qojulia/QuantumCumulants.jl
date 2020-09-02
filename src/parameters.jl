"""
    SymbolicNumber <: Number

Abstract type for all symbolic numbers, i.e. [`Parameter`](@ref), [`Average`](@ref)
and corresponding expression trees.
"""
abstract type SymbolicNumber <: Number end

Base.:(==)(s::SymbolicNumber,x::Number) = false
Base.:(==)(x::Number,s::SymbolicNumber) = false
Base.:(==)(s1::SymbolicNumber,s2::SymbolicNumber) = false
Base.iszero(::SymbolicNumber) = false
Base.isone(::SymbolicNumber) = false

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
Base.hash(t::NumberTerm{T}, h::UInt) where T = hash(t.arguments, hash(t.f, hash(T, h)))

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
end
Parameter(name::Symbol) = Parameter{Number}(name)
Base.:(==)(p::T, q::T) where T<:Parameter = (p.name==q.name)
Base.hash(p::Parameter{T}, h::UInt) where T = hash(p.name, hash(T, h))

# Methods
Base.conj(p::Parameter{<:Real}) = p

for f = [:+,:-,:*,:/,:^]
    @eval Base.$f(a::SymbolicNumber,b::Number) = NumberTerm($f, [a,b])
    @eval Base.$f(a::Number,b::SymbolicNumber) = NumberTerm($f, [a,b])
    @eval Base.$f(a::SymbolicNumber,b::SymbolicNumber) = NumberTerm($f, [a,b])
end
Base.:^(a::SymbolicNumber, b::Int) = NumberTerm(^, [a,b])

const ftrig_ps = [:cos,:sin,:tan,:sqrt,:conj,:exp]
for f = ftrig_ps
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
function substitute(t::NumberTerm, dict; simplify=true, kwargs...)
    if haskey(dict, t)
        return dict[t]
    elseif haskey(dict, t')
        return dict[t]'
    else
        if simplify
            return simplify_constants(t.f([substitute(arg, dict; simplify=simplify) for arg in t.arguments]...), kwargs...)
        else
            return t.f([substitute(arg, dict; simplify=simplify) for arg in t.arguments]...)
        end
    end
end
function substitute(x::SymbolicNumber, dict; kwargs...)
    if haskey(dict, x)
        return dict[x]
    elseif haskey(dict, x')
        return dict[x']'
    else
        return x
    end
end

# Conversion to SymbolicUtils
_to_symbolic(p::Parameter{T}) where T = SymbolicUtils.Sym{T}(p.name)
_to_symbolic(n::NumberTerm{T}) where T = SymbolicUtils.Term{T}(n.f, _to_symbolic.(n.arguments))
function _to_qumulants(s::SymbolicUtils.Sym{T}) where T<:Number
    return Parameter{T}(s.name)
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

"""
    simplify_constants(t::NumberTerm;kwargs...)

Standard simplification for [`SymbolicNumber`](@ref) types. Converts to a
`SymbolicUtils` expression and uses its standard simplification routines for
symbolic number variables.
"""
simplify_constants(s::Number;kwargs...) = s
function simplify_constants(t::NumberTerm;kwargs...)
    s = _to_symbolic(t)
    s_ = SymbolicUtils.simplify(s;kwargs...)
    return _to_qumulants(s_)
end
