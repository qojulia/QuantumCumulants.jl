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
