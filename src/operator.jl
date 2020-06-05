import SymbolicUtils

# Abstract types
abstract type AbstractOperator end

"""
    BasicOperator <: AbstractOperator

Abstract type representing fundamental operator types.
"""
abstract type BasicOperator <: AbstractOperator end

isoperator(x) = false
isoperator(x::Union{T,SymbolicUtils.Symbolic{T}}) where {T<:AbstractOperator} = true

"""
    OperatorTerm <: AbstractOperator

Symbolic expression tree consisting of [`AbstractOperator`](@ref) and `Number`
arguments.
"""
struct OperatorTerm{F,ARGS} <: AbstractOperator
    f::F
    arguments::ARGS
end
Base.:(==)(t1::OperatorTerm,t2::OperatorTerm) = (t1.f===t2.f && t1.arguments==t2.arguments)

for f = [:+,:-,:*]
    @eval Base.$f(a::AbstractOperator,b::AbstractOperator) = (check_hilbert(a,b); OperatorTerm($f, [a,b]))
    @eval Base.$f(a::AbstractOperator,b::Number) = OperatorTerm($f, [a,b])
    @eval Base.$f(a::Number,b::AbstractOperator) = OperatorTerm($f, [a,b])
end
# Base.:^(a::AbstractOperator,b) = OperatorTerm(^, [a,b]) TODO

# Variadic methods
Base.:-(x::AbstractOperator) = -1*x
for f in [:+,:*]
    @eval Base.$f(x::AbstractOperator) = x
    @eval Base.$f(x::AbstractOperator, w::AbstractOperator...) = (check_hilbert(x,w...); OperatorTerm($f, [x;w...]))
    @eval Base.$f(x, y::AbstractOperator, w...) = (check_hilbert(x,y,w...); OperatorTerm($f, [x;y;w...]))
    @eval Base.$f(x::AbstractOperator, y::AbstractOperator, w...) = (check_hilbert(x,y,w...); OperatorTerm($f, [x;y;w...]))
end

Base.adjoint(t::OperatorTerm) = OperatorTerm(t.f, adjoint.(t.arguments))
function Base.adjoint(t::OperatorTerm{<:typeof(*)})
    args = reverse(adjoint.(t.arguments))
    is_c = iscommutative.(*,args)
    args_c = args[is_c]
    args_nc = sort(args[.!is_c], by=acts_on)
    return OperatorTerm(t.f, [args_c;args_nc])
end

# Hilbert space checks
check_hilbert(a::BasicOperator,b::BasicOperator) = (a.hilbert == b.hilbert) || error("Incompatible Hilbert spaces $(a.hilbert) and $(b.hilbert)!")
function check_hilbert(a::OperatorTerm,b::AbstractOperator)
    a_ = findfirst(x->isa(x,AbstractOperator), a.arguments)
    return check_hilbert(a_,b)
end
function check_hilbert(a::AbstractOperator,b::OperatorTerm)
    b_ = findfirst(x->isa(x,AbstractOperator), b.arguments)
    return check_hilbert(a,b_)
end
function check_hilbert(a::OperatorTerm,b::OperatorTerm)
    a_ = findfirst(x->isa(x,AbstractOperator), a.arguments)
    b_ = findfirst(x->isa(x,AbstractOperator), b.arguments)
    return check_hilbert(a_,b_)
end
check_hilbert(args...) = nothing#reduce(check_hilbert, args) # TODO

acts_on(op::BasicOperator) = op.aon
function acts_on(t::OperatorTerm)
    ops = filter(isoperator, t.arguments)
    aon = Int[]
    for op in ops
        append!(aon, acts_on(op))
    end
    unique!(aon)
    sort!(aon)
    return aon
end
acts_on(::Number) = Int[]

get_index(t::BasicOperator) = t.index
function get_index(op::OperatorTerm)
    idx = Int[]
    for arg in op.arguments
        append!(idx, get_index(arg))
    end
    return idx
end
get_index(::Number) = Int[]

Base.one(::T) where T<:AbstractOperator = one(T)
Base.one(::Type{<:AbstractOperator}) = 1
Base.isone(::AbstractOperator) = false
Base.zero(::T) where T<:AbstractOperator = zero(T)
Base.zero(::Type{<:AbstractOperator}) = 0
Base.iszero(::AbstractOperator) = false

function Base.copy(op::T) where T<:BasicOperator
    fields = [getfield(op, n) for n in fieldnames(T)]
    return T(fields...)
end
function Base.copy(t::OperatorTerm)
    return OperatorTerm(t.f, copy.(t.arguments))
end
