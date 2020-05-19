import SymbolicUtils

# Abstract types
abstract type AbstractOperator end
abstract type BasicOperator <: AbstractOperator end

isoperator(x) = false
isoperator(x::Union{T,SymbolicUtils.Symbolic{T}}) where {A,T<:AbstractOperator} = true

# Keep track of operators when converting to SymbolicUtils
const OPERATORS_TO_SYMS = Dict{BasicOperator,SymbolicUtils.Sym}()
const SYMS_TO_OPERATORS = Dict{SymbolicUtils.Sym,BasicOperator}()

# Operator expressions
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
Base.adjoint(t::OperatorTerm{<:typeof(*)}) = OperatorTerm(t.f, reverse(adjoint.(t.arguments)))

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


# Tensor product via embedding operators
struct EmbeddedOperator{H<:ProductSpace,OP<:BasicOperator,A} <: BasicOperator
    hilbert::H
    operator::OP
    aon::A
    function EmbeddedOperator{H,OP,A}(hilbert::H,operator::OP,aon::A) where {H,OP,A}
        op = new(hilbert,operator,aon)
        if !haskey(OPERATORS_TO_SYMS, op)
            sym = generate_symbolic(operator)
            OPERATORS_TO_SYMS[op] = sym
            SYMS_TO_OPERATORS[sym] = op
        end
        return op
    end
end
EmbeddedOperator(hilbert::H,operator::OP,aon::A) where {H,OP,A} = EmbeddedOperator{H,OP,A}(hilbert,operator,aon)
Base.:(==)(a::EmbeddedOperator,b::EmbeddedOperator) = (a.hilbert==b.hilbert && a.aon==b.aon && a.operator==b.operator)
Base.adjoint(op::EmbeddedOperator) = embed(op.hilbert,adjoint(op.operator),op.aon)

acts_on(::BasicOperator) = 1
acts_on(op::EmbeddedOperator) = op.aon
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

# embed into product spaces
function embed(h::ProductSpace,op::BasicOperator,aon::Int)
    check_hilbert(h.spaces[aon],op.hilbert)
    EmbeddedOperator(h,op,aon)
end
function ⊗(a::BasicOperator,b::BasicOperator)
    h = a.hilbert⊗b.hilbert
    a_ = embed(h,a,1)
    b_ = embed(h,b,2)
    return a_*b_
end
function ⊗(a::EmbeddedOperator,b::BasicOperator)
    h = a.hilbert⊗b.hilbert
    a_ = embed(h,a.operator,a.aon)
    b_ = embed(h,b,a.aon+1)
    return a_*b_
end
function ⊗(a::EmbeddedOperator,b::EmbeddedOperator)
    h = a.hilbert⊗b.hilbert
    a_ = embed(h,a.operator,a.aon)
    b_ = embed(h,b.operator,a.aon+b.aon+1)
    return a_*b_
end

Base.one(::T) where T<:AbstractOperator = one(T)
Base.one(::Type{<:AbstractOperator}) = 1
Base.isone(::AbstractOperator) = false
Base.zero(::T) where T<:AbstractOperator = zero(T)
Base.zero(::Type{<:AbstractOperator}) = 0
Base.iszero(::AbstractOperator) = false
