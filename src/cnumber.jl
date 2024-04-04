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
        return SymbolicUtils.setmetadata(s,MTK.MTKVariableTypeCtx,MTK.PARAMETER)
    end
end

# Promoting to CNumber ensures we own the symtype; could be used to dispatch
# on Base methods (e.g. latex printing)
Base.promote_rule(::Type{<:CNumber},::Type{<:Number}) = CNumber

Base.one(::Type{Parameter}) = 1
Base.zero(::Type{Parameter}) = 0
Base.adjoint(x::SymbolicUtils.Symbolic{<:CNumber}) = conj(x)

# TODO: this doesn't work with just setting Complex for some reason; am I doing this right?
MTK.concrete_symtype(::Symbolics.BasicSymbolic{T}) where T <: CNumber = ComplexF64

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

"""
    cnumber(symbols::Symbol)
    cnumber(s::String)

Create symbolic cnumber.

Expamples
=========
```
julia> ps = cnumber(:a)
a

julia> cnumber("a") == ps
true
```
"""
cnumber(s::Symbol) = Parameter(s; metadata=source_metadata(:cnumbers, s))
cnumber(s::String) = cnumber(Symbol(s))


### real paramters ###
"""
    RNumber <: Real

Abstract type for real symbolic numbers [`RealParameter`](@ref).
"""
abstract type RNumber <: Real end

"""
    RealParameter <: RNumber

Type used as symbolic type in a `SymbolicUtils.Sym` variable to represent
a real parameter.
"""
struct RealParameter <: RNumber
    function RealParameter(name; metadata=source_metadata(:RealParameter, name))
        s = SymbolicUtils.Sym{RealParameter}(name)
        s = SymbolicUtils.setmetadata(s, MTK.VariableSource, (:RealParameter, name))
        return SymbolicUtils.setmetadata(s,MTK.MTKVariableTypeCtx,MTK.PARAMETER)
    end
end

# Promoting to RNumber ensures we own the symtype; could be used to dispatch
# on Base methods (e.g. latex printing)
Base.promote_rule(::Type{<:RNumber},::Type{<:Real}) = RNumber

Base.one(::Type{RealParameter}) = 1
Base.zero(::Type{RealParameter}) = 0
Base.adjoint(x::SymbolicUtils.Symbolic{<:RNumber}) = x
Base.adjoint(x::RNumber) = x
Base.conj(x::RNumber) = x

MTK.concrete_symtype(::Symbolics.BasicSymbolic{T}) where T<:RNumber = Real

"""
    @rnumbers(ps...)

Convenience macro to quickly define symbolic rnumbers.

Examples
========
```
julia> @rnumbers ω κ
(ω, κ)
```
"""
macro rnumbers(ps...)
    ex = Expr(:block)
    pnames = []
    for p in ps
        @assert p isa Symbol
        push!(pnames, p)
        d = source_metadata(:rnumbers, p)
        ex_ = Expr(:(=), esc(p), Expr(:call, :RealParameter, Expr(:quote, p), Expr(:kw, :metadata, Expr(:quote, d))))
        push!(ex.args, ex_)
    end
    push!(ex.args, Expr(:tuple, map(esc, pnames)...))
    return ex
end

"""
    rnumbers(symbols::Symbol...)
    rnumbers(s::String)

Create symbolic rnumbers.

Expamples
=========
```
julia> ps = rnumbers(:a, :b)
(a, b)

julia> rnumbers("a b") == ps
true
```
"""
function rnumbers(syms::Symbol...)
    ps = Tuple(RealParameter(s; metadata=source_metadata(:rnumbers, s)) for s in syms)
    return ps
end
function rnumbers(s::String)
    syms = [Symbol(p) for p in split(s, " ")]
    return rnumbers(syms...)
end

"""
    rnumber(symbols::Symbol)
    rnumber(s::String)

Create symbolic rnumber.

Expamples
=========
```
julia> ps = rnumber(:a)
a

julia> rnumber("a") == ps
true
```
"""
rnumber(s::Symbol) = RealParameter(s; metadata=source_metadata(:rnumbers, s))
rnumber(s::String) = rnumber(Symbol(s))

# this should be true for all analytic functions (write as Taylor-series)
function Base.adjoint(x::SymbolicUtils.BasicSymbolic{Complex{RNumber}})
    f = SymbolicUtils.operation(x)
    args = SymbolicUtils.arguments(x)
    return f(conj.(args)...)
end
function Base.conj(x::SymbolicUtils.BasicSymbolic{Complex{RNumber}})
    f = SymbolicUtils.operation(x)
    args = SymbolicUtils.arguments(x)
    return f(conj.(args)...)
end
function Base.adjoint(x::SymbolicUtils.BasicSymbolic{CNumber})
    f = SymbolicUtils.operation(x)
    args = SymbolicUtils.arguments(x)
    return f(conj.(args)...)
end
function Base.conj(x::SymbolicUtils.BasicSymbolic{CNumber})
    f = SymbolicUtils.operation(x)
    args = SymbolicUtils.arguments(x)
    return f(conj.(args)...)
end


const AbstractQCParameter = Union{CNumber, RNumber}

# TODO: real IndexedVariables