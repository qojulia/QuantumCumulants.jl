"""
    CNumber <: Number

Abstract type for all symbolic numbers, i.e. [`Parameter`](@ref), [`Average`](@ref)
and corresponding expression trees.
"""
abstract type CNumber <: Number end

"""
    Parameter <: CNumber

Type used as symbolic type in a `SymbolicUtils.Sym` variable to represent
a parameter.
"""
struct Parameter <: CNumber
    function Parameter(name)
        return SymbolicUtils.Sym{Parameter}(name)
    end
end

# Promoting to CNumber ensures we own the symtype; could be used to dispatch
# on Base methods (e.g. latex printing); not sure in how far this is type piracy
Base.promote_rule(::Type{<:CNumber},::Type{<:Number}) = CNumber

"""
    @cnumbers(ps...)

Convenience macro to quickly define symbolic cnumbers.

Examples
========
```
julia> @cnumbers ω κ
(ω, κ)
```
"""
macro cnumbers(ps...)
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
    cnumbers(symbols::Symbol...)
    cnumbers(s::String)

Create symbolic cnumbers.

Expamples
=========
```
julia> ps = cnumbers(:a, :b)
(a, b)

julia> cnumbers("a b") == ps
true
```
"""
function cnumbers(syms::Symbol...)
    ps = Tuple(Parameter(s) for s in syms)
    return ps
end
function cnumbers(s::String)
    syms = [Symbol(p) for p in split(s, " ")]
    return cnumbers(syms...)
end
