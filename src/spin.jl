"""
    SpinSpace <: HilbertSpace

[`HilbertSpace`](@ref) defining a space for spin operators.
"""
struct SpinSpace{S,N} <: HilbertSpace
    name::S
    spin::N
end
Base.:(==)(h1::T,h2::T) where T<:SpinSpace = (h1.name==h2.name && h1.spin==h2.spin)
Base.hash(h::SpinSpace, i::UInt) = hash(h.name, hash(h.spin, i))
isspin(N,h::SpinSpace,args...) = (N==h.spin)
isspin(N,h::ProductSpace,aon) = isspin(N,h.spaces[aon])

"""
    SigmaX <: BasicOperator

Operator on a [`SpinSpace`](@ref) representing the spin x component.
"""
struct SigmaX{H<:HilbertSpace,S,A} <: BasicOperator
    hilbert::H
    name::S
    aon::A
    function SigmaX{H,S,A}(hilbert::H,name::S,aon::A) where {H,S,A}
        @assert has_hilbert(SpinSpace,hilbert,aon)
        op = new(hilbert,name,aon)
        if !haskey(OPERATORS_TO_SYMS, op)
            sym = SymbolicUtils.Sym{SigmaX}(gensym(:SigmaX))
            OPERATORS_TO_SYMS[op] = sym
            SYMS_TO_OPERATORS[sym] = op
        end
        return op
    end
end
issigmax(a::SymbolicUtils.Sym{T}) where {T<:SigmaX} = true

"""
    SigmaY <: BasicOperator

Operator on a [`SpinSpace`](@ref) representing the spin y component.
"""
struct SigmaY{H<:HilbertSpace,S,A} <: BasicOperator
    hilbert::H
    name::S
    aon::A
    function SigmaY{H,S,A}(hilbert::H,name::S,aon::A) where {H,S,A}
        @assert has_hilbert(SpinSpace,hilbert,aon)
        op = new(hilbert,name,aon)
        if !haskey(OPERATORS_TO_SYMS, op)
            sym = SymbolicUtils.Sym{SigmaY}(gensym(:SigmaY))
            OPERATORS_TO_SYMS[op] = sym
            SYMS_TO_OPERATORS[sym] = op
        end
        return op
    end
end
issigmay(a::SymbolicUtils.Sym{T}) where {T<:SigmaY} = true

"""
    SigmaZ <: BasicOperator

Operator on a [`SpinSpace`](@ref) representing the spin y component.
"""
struct SigmaZ{H<:HilbertSpace,S,A} <: BasicOperator
    hilbert::H
    name::S
    aon::A
    function SigmaZ{H,S,A}(hilbert::H,name::S,aon::A) where {H,S,A}
        @assert has_hilbert(SpinSpace,hilbert,aon)
        op = new(hilbert,name,aon)
        if !haskey(OPERATORS_TO_SYMS, op)
            sym = SymbolicUtils.Sym{SigmaZ}(gensym(:SigmaZ))
            OPERATORS_TO_SYMS[op] = sym
            SYMS_TO_OPERATORS[sym] = op
        end
        return op
    end
end
issigmaz(a::SymbolicUtils.Sym{T}) where {T<:SigmaZ} = true

for f in [:SigmaX,:SigmaY,:SigmaZ]
    @eval $(f)(hilbert::H,name::S,aon::A) where {H,S,A} = $(f){H,S,A}(hilbert,name,aon)
    @eval $(f)(hilbert::SpinSpace,name) = $(f)(hilbert,name,1)
    @eval function $(f)(hilbert::ProductSpace,name)
        i = findall(x->isa(x,SpinSpace),hilbert.spaces)
        if length(i)==1
            return $(f)(hilbert,name,i[1])
        else
            isempty(i) && error("Can only create $($(f)) on SpinSpace! Not included in $(hilbert)")
            length(i)>1 && error("More than one SpinSpace in $(hilbert)! Specify on which Hilbert space $($(f)) should be created with $($(f))(hilbert,name,i)!")
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
    @eval isspin(N,op::T) where T<:($(f)) = isspin(N,op.hilbert,op.aon)
end

# Commutation relation in simplification
function commute_spin(a::SymbolicUtils.Sym{S1},b::SymbolicUtils.Sym{S2}) where {S1<:SigmaY,S2<:SigmaX}
    op = _to_qumulants(a)
    h = hilbert(op)
    aon = acts_on(op)
    z = _to_symbolic(SigmaZ(h,op.name,aon))
    return -im*z + b*a
end
function commute_spin(a::SymbolicUtils.Sym{S1},b::SymbolicUtils.Sym{S2}) where {S1<:SigmaZ,S2<:SigmaY}
    op = _to_qumulants(a)
    h = hilbert(op)
    aon = acts_on(op)
    x = _to_symbolic(SigmaX(h,op.name,aon))
    return -im*x + b*a
end
function commute_spin(a::SymbolicUtils.Sym{S1},b::SymbolicUtils.Sym{S2}) where {S1<:SigmaZ,S2<:SigmaX}
    op = _to_qumulants(a)
    h = hilbert(op)
    aon = acts_on(op)
    y = _to_symbolic(SigmaY(h,op.name,aon))
    return im*y + b*a
end
