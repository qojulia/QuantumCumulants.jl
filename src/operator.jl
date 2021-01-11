import SymbolicUtils

"""
    AbstractOperator

Abstract type representing any expression involving operators.
"""
abstract type AbstractOperator end

"""
    BasicOperator <: AbstractOperator

Abstract type representing fundamental operator types.
"""
abstract type BasicOperator <: AbstractOperator end

isoperator(x) = false
isoperator(x::Union{T,SymbolicUtils.Symbolic{T}}) where {A,T<:AbstractOperator} = true
Base.:(==)(a::T,b::T) where T<:BasicOperator = (a.hilbert==b.hilbert && a.name==b.name && a.aon==b.aon)

# Dicts for conversion
const OPERATORS_TO_SYMS = Dict{BasicOperator,SymbolicUtils.Sym}()
const SYMS_TO_OPERATORS = Dict{SymbolicUtils.Sym,BasicOperator}()

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
Base.hash(t::OperatorTerm, h::UInt) = hash(t.arguments, hash(t.f, h))

for f = [:+,:-,:*]
    @eval Base.$f(a::AbstractOperator,b::AbstractOperator) = (check_hilbert(a,b); OperatorTerm($f, [a,b]))
    @eval Base.$f(a::AbstractOperator,b::Number) = OperatorTerm($f, [a,b])
    @eval Base.$f(a::Number,b::AbstractOperator) = OperatorTerm($f, [a,b])
end
Base.:^(a::AbstractOperator,b::Integer) = OperatorTerm(^, [a,b])
Base.:/(a::AbstractOperator,b::Number) = OperatorTerm(/, [a,b])

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
    is_c = iscommutative.(args)
    args_c = args[is_c]
    args_nc = sort(args[.!is_c], lt=lt_aon)
    return OperatorTerm(t.f, [args_c;args_nc])
end

# Hilbert space checks
check_hilbert(a::BasicOperator,b::BasicOperator) = (a.hilbert == b.hilbert) || error("Incompatible Hilbert spaces $(a.hilbert) and $(b.hilbert)!")
function check_hilbert(a::OperatorTerm,b::BasicOperator)
    a_ = findfirst(x->isa(x,AbstractOperator), a.arguments)
    return check_hilbert(a_,b)
end
function check_hilbert(a::BasicOperator,b::OperatorTerm)
    b_ = findfirst(x->isa(x,AbstractOperator), b.arguments)
    return check_hilbert(a,b_)
end
function check_hilbert(a::OperatorTerm,b::OperatorTerm)
    a_ = findfirst(x->isa(x,AbstractOperator), a.arguments)
    b_ = findfirst(x->isa(x,AbstractOperator), b.arguments)
    return check_hilbert(a_,b_)
end
function check_hilbert(args...)
    for i=1:length(args)-1
        check_hilbert(args[i], args[i+1])
    end
end
check_hilbert(x,y) = true

"""
    acts_on(op::AbstractOperator)

Shows on which Hilbert space `op` acts. For [`BasicOperator`](@ref) types, this
returns an Integer, whereas for a [`OperatorTerm`](@ref) it returns a `Vector{Int}`
whose entries specify all subspaces on which the expression acts.
"""
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
acts_on(x) = Int[]

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
