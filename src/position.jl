struct PositionSpace{S} <: HilbertSpace
    name::S
end
Base.:(==)(h1::T,h2::T) where T<:PositionSpace = (h1.name==h2.name)

struct Position{H<:HilbertSpace,S,A} <: BasicOperator
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

struct Momentum{H<:HilbertSpace,S,A} <: BasicOperator
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

for f in [:Position,:Momentum]
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
commute_xp(f,fargs) = (x,p) -> commute_xp(f, fargs, x, p)
commute_xp(::Nothing, fargs, x, p) = im + p*x
function commute_xp(f::typeof(cos), fargs, x, p)
    c = *(fargs...)
    return c*sin(c*x) + p*cos(c*x)
end
function commute_xp(f::typeof(sin), fargs, x, p)
    c = *(fargs...)
    return -1*c*cos(c*x) + p*sin(c*x)
end
commute_xx(f,fargs) = (x,p) -> commute_xx(f, fargs, x, p)
for F in trig
    @eval function commute_xx(f::typeof($(F)), fargs, x, p)
        c = *(fargs...)
        return p*f(c*x)
    end
end
