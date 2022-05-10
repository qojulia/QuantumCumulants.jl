"""
    QNumber

Abstract type representing any expression involving operators.
"""
abstract type QNumber end

"""
    QSym <: QNumber

Abstract type representing fundamental operator types.
"""
abstract type QSym <: QNumber end

# Generic hash fallback for interface -- this will be slow
function Base.hash(op::T, h::UInt) where T<:QSym
    n = fieldcount(T)
    if n == 3
        # These three fields need to be defined for any QSym
        return hash(T, hash(op.hilbert, hash(op.name, hash(op.aon, h))))
    else
        # If there are more we'll need to iterate through
        h_ = copy(h)
        for k = n:-1:4
            if fieldname(typeof(op), k) !== :metadata
                h_ = hash(getfield(op, k), h_)
            end
        end
        return hash(T, hash(op.hilbert, hash(op.name, hash(op.aon, h_))))
    end
end

"""
    QTerm <: QNumber

Abstract type representing noncommutative expressions.
"""
abstract type QTerm <: QNumber end

Base.isless(a::QSym, b::QSym) = a.name < b.name

## Interface for SymbolicUtils

TermInterface.exprhead(::QNumber) = :call
TermInterface.istree(::QSym) = false
TermInterface.istree(::QTerm) = true
TermInterface.istree(::Type{T}) where {T<:QTerm} = true

# Symbolic type promotion
SymbolicUtils.promote_symtype(f, Ts::Type{<:QNumber}...) = promote_type(Ts...)
SymbolicUtils.promote_symtype(f, T::Type{<:QNumber}, Ts...) = T
SymbolicUtils.promote_symtype(f,T::Type{<:QNumber},S::Type{<:Number}) = T
SymbolicUtils.promote_symtype(f,T::Type{<:Number},S::Type{<:QNumber}) = S
SymbolicUtils.promote_symtype(f,T::Type{<:QNumber},S::Type{<:QNumber}) = promote_type(T,S)

SymbolicUtils.symtype(x::T) where T<:QNumber = T

# Standard simplify
function SymbolicUtils.simplify(x::QNumber;kwargs...)
    avg = average(x)
    avg_ = SymbolicUtils.simplify(avg;kwargs...)
    return undo_average(avg_)
end

## End of interface

## Methods
import Base: *, +, -
const SNuN = Union{<:SymbolicUtils.Symbolic{<:Number}, <:Number}

Base.:~(a::QNumber, b::QNumber) = Symbolics.Equation(a, b)

## Multiplication
"""
    QMul <: QTerm

Represent a multiplication involving [`QSym`](@ref) types.

Fields:
======

* arg_c: The commutative prefactor.
* args_nc: A vector containing all [`QSym`](@ref) types.
"""
struct QMul{M} <: QTerm
    arg_c
    args_nc::Vector{Any}
    metadata::M
    function QMul{M}(arg_c, args_nc, metadata) where {M}
        if SymbolicUtils._isone(arg_c) && length(args_nc)==1
            return args_nc[1]
        else
            return new(arg_c, args_nc, metadata)
        end
    end
end
QMul(arg_c, args_nc; metadata::M=NO_METADATA) where {M} = QMul{M}(arg_c, args_nc, metadata)
Base.hash(q::QMul, h::UInt) = hash(QMul, hash(q.arg_c, SymbolicUtils.hashvec(q.args_nc, h)))
#Base.isless(a::QMul, b::QMul) = isless(a.h, b.h)

SymbolicUtils.operation(::QMul) = (*)
SymbolicUtils.arguments(a::QMul) = vcat(a.arg_c, a.args_nc)

function SymbolicUtils.similarterm(::QMul, ::typeof(*), args; metadata=NO_METADATA, exprhead=nothing)
    args_c = filter(x->!(x isa QNumber), args)
    args_nc = filter(x->x isa QNumber, args)
    arg_c = *(args_c...)
    return QMul(arg_c, args_nc; metadata)
end

SymbolicUtils.metadata(a::QMul) = a.metadata

function Base.adjoint(q::QMul)
    args_nc = map(adjoint, q.args_nc)
    reverse!(args_nc)
    sort!(args_nc, by=acts_on)
    return QMul(conj(q.arg_c), args_nc; q.metadata)
end


function Base.isequal(a::QMul, b::QMul)
    isequal(a.arg_c, b.arg_c) || return false
    length(a.args_nc)==length(b.args_nc) || return false
    for (arg_a, arg_b) ∈ zip(a.args_nc, b.args_nc)
        isequal(arg_a,arg_b) || return false
    end
    return true
end

function *(a::QSym,b::QSym)
    check_hilbert(a, b)
    args = [a,b]
    sort!(args, by=acts_on)
    QMul(1,args)
end

function *(a::QSym, b::SNuN)
    SymbolicUtils._iszero(b) && return b
    SymbolicUtils._isone(b) && return a
    return QMul(b,[a])
end
*(b::SNuN, a::QNumber) = a*b

function *(a::QMul, b::SNuN)
    SymbolicUtils._iszero(b) && return b
    SymbolicUtils._isone(b) && return a
    arg_c = a.arg_c * b
    return QMul(arg_c, a.args_nc)
end

function *(a::QSym, b::QMul)
    check_hilbert(a, b)
    args_nc = vcat(a,b.args_nc)
    sort!(args_nc, by=acts_on)
    return merge_commutators(b.arg_c,args_nc)
end
function *(a::QMul, b::QSym)
    check_hilbert(a, b)
    args_nc = vcat(a.args_nc, b)
    sort!(args_nc, by=acts_on)
    return merge_commutators(a.arg_c,args_nc)
end

function *(a::QMul, b::QMul)
    check_hilbert(a, b)
    args_nc = vcat(a.args_nc, b.args_nc)
    sort!(args_nc, by=acts_on)
    arg_c = a.arg_c*b.arg_c
    return merge_commutators(arg_c,args_nc)
end

Base.:/(a::QNumber, b::SNuN) = (1/b) * a

function merge_commutators(arg_c,args_nc)
    i = 1
    was_merged = false
    while i<length(args_nc)
        if _ismergeable(args_nc[i], args_nc[i+1])
            args_nc[i] = *(args_nc[i], args_nc[i+1])
            iszero(args_nc[i]) && return 0
            deleteat!(args_nc, i+1)
            was_merged = true
        end
        i += 1
    end
    if was_merged
        return *(arg_c, args_nc...)
    else
        return QMul(arg_c, args_nc)
    end
end

_ismergeable(a,b) = isequal(acts_on(a),acts_on(b)) && ismergeable(a,b) && isequal(hilbert(a),hilbert(b))
ismergeable(a,b) = false

## Powers
function Base.:^(a::QNumber, n::Integer)
    iszero(n) && return 1
    isone(n) && return a
    return *((a for i=1:n)...)
end

## Addition
"""
    QAdd <: QTerm

Represent an addition involving [`QNumber`](@ref) and other types.
"""
struct QAdd <: QTerm
    arguments::Vector{Any}
end

Base.hash(q::T, h::UInt) where T<:QAdd = hash(T, SymbolicUtils.hashvec(q.arguments, h))
function Base.isequal(a::QAdd,b::QAdd)
    length(a.arguments)==length(b.arguments) || return false
    for (arg_a,arg_b) ∈ zip(a.arguments, b.arguments)
        isequal(arg_a, arg_b) || return false
    end
    return true
end

SymbolicUtils.operation(::QAdd) = (+)
SymbolicUtils.arguments(a::QAdd) = a.arguments
SymbolicUtils.similarterm(::QAdd, ::typeof(+), args; metadata=NO_METADATA, exprhead=nothing) = QAdd(args; metadata)

SymbolicUtils.metadata(q::QAdd) = q.metadata

Base.adjoint(q::QAdd) = QAdd(map(adjoint, q.arguments))

-(a::QNumber) = -1*a
-(a,b::QNumber) = a + (-b)
-(a::QNumber,b) = a + (-b)
-(a::QNumber,b::QNumber) = a + (-b)

+(a::QNumber, b::SNuN) = QAdd([a,b])
+(a::SNuN,b::QNumber) = +(b,a)
function +(a::QAdd,b::SNuN)
    SymbolicUtils._iszero(b) && return a
    args = vcat(a.arguments, b)
    return QAdd(args)
end

function +(a::QNumber, b::QNumber)
    check_hilbert(a, b)
    args = [a,b]
    return QAdd(args)
end

function +(a::QAdd,b::QNumber)
    check_hilbert(a, b)
    args = vcat(a.arguments, b)
    return QAdd(args)
end
function +(b::QNumber,a::QAdd)
    check_hilbert(a, b)
    args = vcat(a.arguments, b)
    return QAdd(args)
end
function +(a::QAdd,b::QAdd)
    check_hilbert(a, b)
    args = vcat(a.arguments, b.arguments)
    return QAdd(args)
end

function *(a::QAdd, b)
    check_hilbert(a, b)
    args = Any[a_ * b for a_ ∈ a.arguments]
    flatten_adds!(args)
    q = QAdd(args)
    return q
end
function *(a::QNumber, b::QAdd)
    check_hilbert(a, b)
    args = Any[a * b_ for b_ ∈ b.arguments]
    flatten_adds!(args)
    q = QAdd(args)
    return q
end

function *(a::QAdd, b::QAdd)
    check_hilbert(a, b)
    args = []
    for a_ ∈ a.arguments, b_ ∈ b.arguments
        push!(args, a_ * b_)
    end
    flatten_adds!(args)
    q = QAdd(args)
    return q
end

function flatten_adds!(args)
    i = 1
    while i <= length(args)
        if args[i] isa QAdd
            append!(args,args[i].arguments)
            deleteat!(args, i)
        end
        i += 1
    end
    return args
end



## Hilbert space checks
check_hilbert(a::QNumber,b::QNumber) = (hilbert(a) == hilbert(b)) || error("Incompatible Hilbert spaces $(a.hilbert) and $(b.hilbert)!")
check_hilbert(x,y) = nothing

hilbert(a::QSym) = a.hilbert
hilbert(a::QMul) = hilbert(a.args_nc[1])
function hilbert(a::QAdd)
    idx = findfirst(x->x isa QNumber, a.arguments)
    hilbert(a.arguments[idx])
end

const AonType = Union{<:Int,<:ClusterAon}
"""
    acts_on(op)

Shows on which Hilbert space `op` acts. For [`QSym`](@ref) types, this
returns an Integer, whereas for a `Term` it returns a `Vector{Int}`
whose entries specify all subspaces on which the expression acts.
"""
acts_on(op::QSym) = op.aon
function acts_on(q::QMul)
    aon = AonType[]
    for arg ∈ q.args_nc
        aon_ = acts_on(arg)
        aon_ ∈ aon || push!(aon, aon_)
    end
    return aon
end
function acts_on(q::QAdd)
    aon = AonType[]
    for arg ∈ q.arguments
        append!(aon, acts_on(arg))
    end
    unique!(aon)
    sort!(aon)
    return aon
end
acts_on(x) = Int[]

Base.one(::T) where T<:QNumber = one(T)
Base.one(::Type{<:QNumber}) = 1
Base.isone(::QNumber) = false
Base.zero(::T) where T<:QNumber = zero(T)
Base.zero(::Type{<:QNumber}) = 0
Base.iszero(::QNumber) = false

"""
    @qnumbers

Convenience macro for the construction of operators.

Examples
========
```
julia> h = FockSpace(:fock)
ℋ(fock)

julia> @qnumbers a::Destroy(h)
(a,)

julia> h = FockSpace(:one) ⊗ FockSpace(:two)
ℋ(one) ⊗ ℋ(two)

julia> @qnumbers b::Destroy(h,2)
(b,)
```
"""
macro qnumbers(qs...)
    ex = Expr(:block)
    qnames = []
    for q in qs
        @assert q isa Expr && q.head==:(::)
        q_ = q.args[1]
        @assert q_ isa Symbol
        push!(qnames, q_)
        f = q.args[2]
        @assert f isa Expr && f.head==:call
        op = _make_operator(q_, f.args...)
        ex_ = Expr(:(=), esc(q_), op)
        push!(ex.args, ex_)
    end
    push!(ex.args, Expr(:tuple, map(esc, qnames)...))
    return ex
end

function _make_operator(name, T, h, args...)
    name_ = Expr(:quote, name)
    d = source_metadata(:qnumbers, name)
    return Expr(:call, T, esc(h), name_, args..., Expr(:kw, :metadata, Expr(:quote, d)))
end
