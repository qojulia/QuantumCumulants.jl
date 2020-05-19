abstract type SymbolicNumber <: Number end

Base.:(==)(s::SymbolicNumber,x::Number) = false
Base.:(==)(x::Number,s::SymbolicNumber) = false
Base.:(==)(s1::SymbolicNumber,s2::SymbolicNumber) = false

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
Base.:(==)(p::T, q::T) where T<:Parameter = (p.name==q.name)

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
_to_symbolic(p::Parameter{T}) where T = SymbolicUtils.Sym{T}(p.name)
_to_symbolic(n::NumberTerm{T}) where T = SymbolicUtils.Term{T}(n.f, _to_symbolic.(n.arguments))
function _to_qumulants(s::SymbolicUtils.Sym{T}) where T<:Number
    return Parameter{T}(s.name)
end

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

simplify_constants(s::SymbolicNumber;kwargs...) = s
function simplify_constants(t::NumberTerm;kwargs...)
    s = _to_symbolic(t)
    s_ = SymbolicUtils.simplify(s;kwargs...)
    return _to_qumulants(s_)
end