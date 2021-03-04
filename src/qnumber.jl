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

Base.isequal(a::T,b::T) where T<:QSym = isequal(a.hilbert, b.hilbert) && isequal(a.name, b.name) && isequal(a.aon, b.aon)
Base.isless(a::QSym,b::QSym) = a.name < b.name

"""
    QTerm <: QNumber

Symbolic expression tree consisting of [`QNumber`](@ref) and `Number`
arguments.
"""
struct QTerm{F,ARGS} <: QNumber
    f::F
    arguments::ARGS
end
function Base.isequal(t1::QTerm,t2::QTerm)
    t1.f===t2.f || return false
    length(t1.arguments)==length(t2.arguments) || return false
    for (a,b) ∈ zip(t1.arguments, t2.arguments)
        isequal(a,b) || return false
    end
    return true
end
Base.hash(t::QTerm, h::UInt) = hash(t.arguments, hash(t.f, h))


### Interface for SymbolicUtils

# Symbolic type promotion
SymbolicUtils.promote_symtype(f, Ts::Type{<:QNumber}...) = promote_type(QNumber,Ts...)
SymbolicUtils.promote_symtype(f, T::Type{<:QNumber}, Ts...) = promote_type(QNumber,T)
SymbolicUtils.promote_symtype(f,T::Type{<:QNumber},S::Type{<:Number}) = QNumber
SymbolicUtils.promote_symtype(f,T::Type{<:Number},S::Type{<:QNumber}) = QNumber
SymbolicUtils.promote_symtype(f,T::Type{<:QNumber},S::Type{<:QNumber}) = QNumber

# Type promotion for average(::QNumber)::Average
SymbolicUtils.promote_symtype(f::SymbolicUtils.Sym{SymbolicUtils.FnType{Tuple{<:QNumber},T}},
                                ::Type{<:QNumber}) where T = T


SymbolicUtils.symtype(x::T) where T<:QNumber = T
SymbolicUtils.to_symbolic(x::QNumber) = x

SymbolicUtils.istree(::QSym) = false
SymbolicUtils.istree(::QTerm) = true
SymbolicUtils.arguments(t::QTerm) = t.arguments
for f in [*,+,-,/,^]
    @eval SymbolicUtils.operation(t::QTerm{<:$(typeof(f))}) = $f
end
SymbolicUtils.similarterm(t::QTerm{F}, f::F, args::A) where {F,A} = QTerm{F,A}(t.f, args)

### End of interface


### Methods

for f = [:+,:-,:*]
    @eval Base.$f(a::QNumber,b::QNumber) = (check_hilbert(a,b); QTerm($f, [a,b]))
    @eval Base.$f(a::QNumber,b::Number) = QTerm($f, [a,b])
    @eval Base.$f(a::Number,b::QNumber) = QTerm($f, [a,b])
    @eval Base.$f(a::QNumber,b::SymbolicUtils.Symbolic{<:Number}) = QTerm($f, [a,b])
    @eval Base.$f(a::SymbolicUtils.Symbolic{<:Number},b::QNumber) = QTerm($f, [a,b])
end
Base.:^(a::QNumber,b::Integer) = QTerm(^, [a,b])
Base.:/(a::QNumber,b::Number) = QTerm(/, [a,b])
Base.:/(a::QNumber,b::SymbolicUtils.Symbolic{<:Number}) = QTerm(/, [a,b])

# Variadic methods
Base.:-(x::QNumber) = -1*x
for f in [:+,:*]
    @eval Base.$f(x::QNumber) = x
    @eval Base.$f(x::QNumber, w::QNumber...) = (check_hilbert(x,w...); QTerm($f, [x;w...]))
    @eval Base.$f(x, y::QNumber, w...) = (check_hilbert(x,y,w...); QTerm($f, [x;y;w...]))
    @eval Base.$f(x::QNumber, y::QNumber, w...) = (check_hilbert(x,y,w...); QTerm($f, [x;y;w...]))
end

Base.adjoint(t::QTerm) = QTerm(t.f, adjoint.(t.arguments))
function Base.adjoint(t::QTerm{<:typeof(*)})
    args = reverse(adjoint.(t.arguments))
    is_c = iscommutative.(args)
    args_c = args[is_c]
    args_nc = sort(args[.!is_c], lt=lt_aon)
    return QTerm(t.f, [args_c;args_nc])
end

# Hilbert space checks
check_hilbert(a::QSym,b::QSym) = (a.hilbert == b.hilbert) || error("Incompatible Hilbert spaces $(a.hilbert) and $(b.hilbert)!")
function check_hilbert(a::QTerm,b::QSym)
    a_ = findfirst(x->isa(x,QNumber), a.arguments)
    return check_hilbert(a_,b)
end
function check_hilbert(a::QSym,b::QTerm)
    b_ = findfirst(x->isa(x,QNumber), b.arguments)
    return check_hilbert(a,b_)
end
function check_hilbert(a::QTerm,b::QTerm)
    a_ = findfirst(x->isa(x,QNumber), a.arguments)
    b_ = findfirst(x->isa(x,QNumber), b.arguments)
    return check_hilbert(a_,b_)
end
function check_hilbert(args...)
    for i=1:length(args)-1
        check_hilbert(args[i], args[i+1])
    end
end
check_hilbert(x,y) = true


"""
    acts_on(op::QNumber)

Shows on which Hilbert space `op` acts. For [`QSym`](@ref) types, this
returns an Integer, whereas for a [`QTerm`](@ref) it returns a `Vector{Int}`
whose entries specify all subspaces on which the expression acts.
"""
acts_on(op::QSym) = op.aon # TODO make Int[]
function acts_on(t::QTerm)
    ops = filter(SymbolicUtils.sym_isa(QNumber), t.arguments)
    aon = Int[]
    for op in ops
        append!(aon, acts_on(op))
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

function Base.copy(op::T) where T<:QSym
    fields = [getfield(op, n) for n in fieldnames(T)]
    return T(fields...)
end
function Base.copy(t::QTerm)
    return QTerm(t.f, copy.(t.arguments))
end

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
    return Expr(:call, T, esc(h), name_, args...)
end
