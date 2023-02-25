"""
    CNumber <: Number

Abstract type for all symbolic numbers, i.e. [`Parameter`](@ref), [`average`](@ref).
"""
abstract type CNumber <: Number end

"""
    Parameter <: CNumber

Type used as symbolic type in a `SymbolicUtils.Sym` variable to represent
a parameter.
"""
struct Parameter <: CNumber
    function Parameter(name; metadata=source_metadata(:Parameter, name))
        s = SymbolicUtils.Sym{Parameter}(name)
        s = SymbolicUtils.setmetadata(s, MTK.VariableSource, (:Parameter, name))
        # return SymbolicUtils.setmetadata(s, MTK.MTKParameterCtx, true) # this is not the case anymore for the current SymbolicUtils + MTK versions
        return SymbolicUtils.setmetadata(s,MTK.MTKVariableTypeCtx,MTK.PARAMETER) #this works with current version of MTK
    end
end

# Promoting to CNumber ensures we own the symtype; could be used to dispatch
# on Base methods (e.g. latex printing)
Base.promote_rule(::Type{<:CNumber},::Type{<:Number}) = CNumber

Base.one(::Type{Parameter}) = 1
Base.zero(::Type{Parameter}) = 0
Base.adjoint(x::SymbolicUtils.Symbolic{<:CNumber}) = conj(x)

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
        d = source_metadata(:cnumbers, p)
        ex_ = Expr(:(=), esc(p), Expr(:call, :Parameter, Expr(:quote, p), Expr(:kw, :metadata, Expr(:quote, d))))
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
    ps = Tuple(Parameter(s; metadata=source_metadata(:cnumbers, s)) for s in syms)
    return ps
end
function cnumbers(s::String)
    syms = [Symbol(p) for p in split(s, " ")]
    return cnumbers(syms...)
end
