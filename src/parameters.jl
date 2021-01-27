"""
    SymbolicNumber <: Number

Abstract type for all symbolic numbers, i.e. [`Parameter`](@ref), [`Average`](@ref)
and corresponding expression trees.
"""
abstract type SymbolicNumber{T} <: Number end
Base.isequal(::SymbolicNumber,::Number) = false
Base.isequal(::Number,::SymbolicNumber) = false
Base.isequal(::SymbolicNumber,::SymbolicNumber) = false

struct Parameter{T} <: SymbolicNumber{T}
    function Parameter{T}(name) where T<:Number
        return SymbolicUtils.Sym{Parameter{T}}(name)
    end
end
Parameter(name) = Parameter{Number}(name)

# Promoting to SymbolicNumber ensures we own the symtype; could be used to dispatch
# on Base methods (e.g. latex printing); not sure in how far this is type piracy
# Base.promote_rule(::Type{<:SymbolicNumber},::Type{<:Number}) = SymbolicNumber

"""
    @params(ps...)

Convenience macro to quickly define symbolic params.

Examples
========
```
julia> @params ω κ
(ω, κ)
```
"""
macro params(ps...)
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
    params(symbols::Symbol...)
    params(s::String)

Create symbolic params.

Expamples
=========
```
julia> ps = params(:a, :b)
(a, b)

julia> params("a b") == ps
true
```
"""
function params(syms::Symbol...)
    ps = Tuple(Parameter{Number}(s) for s in syms)
    return ps
end
function params(s::String)
    syms = [Symbol(p) for p in split(s, " ")]
    return params(syms...)
end

# """
#     simplify_constants(t::NumberTerm;kwargs...)
#
# Standard simplification for [`SymbolicNumber`](@ref) types. Converts to a
# `SymbolicUtils` expression and uses its standard simplification routines for
# symbolic number variables.
# """
# simplify_constants(s::Number;kwargs...) = s
# function simplify_constants(t::NumberTerm;kwargs...)
#     s = _to_symbolic(t)
#     s_ = SymbolicUtils.simplify(s;kwargs...)
#     return _to_qumulants(s_)
# end
