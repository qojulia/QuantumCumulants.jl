struct PositionSpace{S} <: HilbertSpace
    name::S
end
Base.:(==)(h1::T,h2::T) where T<:PositionSpace = (h1.name==h2.name)

abstract type AbstractPosition <: BasicOperator end
isabstractposition(s::SymbolicUtils.Sym{T}) where T<:AbstractPosition = true

abstract type AbstractMomentum <: BasicOperator end
isabstractmomentum(s::SymbolicUtils.Sym{T}) where T<:AbstractMomentum = true

struct Position{H<:HilbertSpace,S,A} <: AbstractPosition
    hilbert::H
    name::S
    aon::A
    function Position{H,S,A}(hilbert::H,name::S,aon::A) where {H,S,A}
        @assert has_hilbert(PositionSpace,hilbert,aon)
        op = new(hilbert,name,aon)
        if !haskey(OPERATORS_TO_SYMS, op)
            sym = SymbolicUtils.Sym{Position}(gensym(:Position))
            OPERATORS_TO_SYMS[op] = sym
            SYMS_TO_OPERATORS[sym] = op
        end
        return op
    end
end
isposition(s::SymbolicUtils.Sym{T}) where T<:Position = true

struct SymmetrizedPosition{H<:HilbertSpace,S,A} <: AbstractPosition
    hilbert::H
    name::S
    aon::A
    function SymmetrizedPosition{H,S,A}(hilbert::H,name::S,aon::A) where {H,S,A}
        @assert has_hilbert(PositionSpace,hilbert,aon)
        op = new(hilbert,name,aon)
        if !haskey(OPERATORS_TO_SYMS, op)
            sym = SymbolicUtils.Sym{SymmetrizedPosition}(gensym(:SymPosition))
            OPERATORS_TO_SYMS[op] = sym
            SYMS_TO_OPERATORS[sym] = op
        end
        return op
    end
end
issymmetrizedposition(s::SymbolicUtils.Sym{T}) where T<:SymmetrizedPosition = true
symmetrize(x::Position) = SymmetrizedPosition(x.hilbert,x.name,x.aon)
symmetrize(x::SymmetrizedPosition) = x
unsymmetrize(x::SymmetrizedPosition) = Position(x.hilbert,x.name,x.aon)

struct Momentum{H<:HilbertSpace,S,A} <: AbstractMomentum
    hilbert::H
    name::S
    aon::A
    function Momentum{H,S,A}(hilbert::H,name::S,aon::A) where {H,S,A}
        @assert has_hilbert(PositionSpace,hilbert,aon)
        op = new(hilbert,name,aon)
        if !haskey(OPERATORS_TO_SYMS, op)
            sym = SymbolicUtils.Sym{Momentum}(gensym(:Momentum))
            OPERATORS_TO_SYMS[op] = sym
            SYMS_TO_OPERATORS[sym] = op
        end
        return op
    end
end
ismomentum(s::SymbolicUtils.Sym{T}) where T<:Momentum = true

struct SymmetrizedMomentum{H<:HilbertSpace,S,A} <: AbstractMomentum
    hilbert::H
    name::S
    aon::A
    function SymmetrizedMomentum{H,S,A}(hilbert::H,name::S,aon::A) where {H,S,A}
        @assert has_hilbert(PositionSpace,hilbert,aon)
        op = new(hilbert,name,aon)
        if !haskey(OPERATORS_TO_SYMS, op)
            sym = SymbolicUtils.Sym{SymmetrizedMomentum}(gensym(:SymMomentum))
            OPERATORS_TO_SYMS[op] = sym
            SYMS_TO_OPERATORS[sym] = op
        end
        return op
    end
end
issymmetrizedmomentum(s::SymbolicUtils.Sym{T}) where T<:SymmetrizedMomentum = true
symmetrize(p::Momentum) = SymmetrizedMomentum(p.hilbert,p.name,p.aon)
symmetrize(p::SymmetrizedMomentum) = p
unsymmetrize(x::SymmetrizedMomentum) = Momentum(x.hilbert,x.name,x.aon)

for f in [:symmetrize,:unsymmetrize]
    @eval $(f)(s::SymbolicUtils.Symbolic) = _to_symbolic($(f)(_to_qumulants(s)))
    @eval $(f)(x) = x
    @eval function $(f)(x::OperatorTerm)
        args_ = [$(f)(arg) for arg in x.arguments]
        return x.f(args_...)
    end
end

for f in [:Position,:Momentum,:SymmetrizedPosition,:SymmetrizedMomentum]
    @eval $(f)(hilbert::H,name::S,aon::A) where {H,S,A} = $(f){H,S,A}(hilbert,name,aon)
    @eval $(f)(hilbert::PositionSpace,name) = $(f)(hilbert,name,1)
    @eval function $(f)(hilbert::ProductSpace,name)
        i = findall(x->isa(x,PositionSpace),hilbert.spaces)
        if length(i)==1
            return $(f)(hilbert,name,i[1])
        else
            isempty(i) && error("Can only create $($(f)) on PositionSpace! Not included in $(hilbert)")
            length(i)>1 && error("More than one PositionSpace in $(hilbert)! Specify on which Hilbert space $($(f)) should be created with $($(f))(hilbert,name,i)!")
        end
    end
    @eval function embed(h::ProductSpace,op::T,aon::Int) where T<:($(f))
        check_hilbert(h.spaces[aon],op.hilbert)
        op_ = $(f)(h,op.name,aon)
        return op_
    end
    @eval function Base.hash(op::T, h::UInt) where T<:($(f))
        hash(op.hilbert, hash(op.name, hash(op.aon, hash($(f), h))))
    end
    @eval Base.adjoint(op::T) where T<:($(f)) = op
end

# Commutation relation in simplification
function commute_symmetric(x::SymbolicUtils.Symbolic{<:AbstractPosition},p::SymbolicUtils.Symbolic{<:AbstractMomentum})
    r = symmetrize(x)
    q = symmetrize(p)
    return 0.5*r*q + 0.5*q*r + 0.5im
end
function commute_symmetric(p::SymbolicUtils.Symbolic{<:AbstractMomentum},x::SymbolicUtils.Symbolic{<:AbstractPosition})
    r = symmetrize(x)
    q = symmetrize(p)
    return 0.5*r*q + 0.5*q*r - 0.5im
end

function commute_nonsymmetric(x::SymbolicUtils.Symbolic{<:AbstractPosition},p::SymbolicUtils.Symbolic{<:AbstractMomentum})
    r = symmetrize(x)
    q = symmetrize(p)
    return q*r + im
end
function commute_nonsymmetric(p::SymbolicUtils.Symbolic{<:AbstractMomentum},x::SymbolicUtils.Symbolic{<:AbstractPosition})
    r = symmetrize(x)
    q = symmetrize(p)
    return r*q - im
end


# commute_xp(f,fargs) = (x,p) -> commute_xp(f, fargs, x, p)
# commute_xp(::Nothing, fargs, x, p) = im + p*x
# function commute_xp(f::typeof(cos), fargs, x, p)
#     c = *(fargs...)
#     return c*sin(c*x) + p*cos(c*x)
# end
# function commute_xp(f::typeof(sin), fargs, x, p)
#     c = *(fargs...)
#     return -1*c*cos(c*x) + p*sin(c*x)
# end
# commute_xx(f,fargs) = (x,p) -> commute_xx(f, fargs, x, p)
# for F in trig
#     @eval function commute_xx(f::typeof($(F)), fargs, x, p)
#         c = *(fargs...)
#         return p*f(c*x)
#     end
# end
