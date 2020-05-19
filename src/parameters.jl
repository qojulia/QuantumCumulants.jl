abstract type SymbolicNumber <: Number end

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

struct Parameter{T<:Number} <: SymbolicNumber
    name::Symbol
end
Parameter(name::Symbol) = Parameter{Number}(name)

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

# Base.promote_rule(::Type{T},::Type{S}) where {T<:Number,S<:SymbolicNumber} = S
# Base.promote_rule(::Type{T},::Type{S}) where {T<:AbstractOperator,S<:SymbolicNumber} = T

# Conversion to SymbolicUtils
_to_symbolic(p::Parameter{T}) where T = SymbolicUtils.Sym{T}(p.name)
_to_symbolic(n::NumberTerm{T}) where T = SymbolicUtils.Term{T}(n.f, _to_symbolic.(n.arguments))
function _to_qumulants(s::SymbolicUtils.Sym{T}) where T<:Number
    return Parameter{T}(s.name)
end
_to_qumulants(t::SymbolicUtils.Term{T}) where T<:Number = NumberTerm{T}(t.f, _to_qumulants.(t.arguments))


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

function parameters(syms::Symbol...)
    ps = Tuple(Parameter{Number}(s) for s in syms)
    return ps
end
function parameters(s::String)
    syms = [Symbol(p) for p in split(s, " ")]
    return parameters(syms...)
end
