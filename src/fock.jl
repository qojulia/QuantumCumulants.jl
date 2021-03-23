"""
    FockSpace <: HilbertSpace

[`HilbertSpace`](@ref) defining a Fock space for bosonic operators.
See also: [`Destroy`](@ref), [`Create`](@ref)
"""
struct FockSpace{S} <: HilbertSpace
    name::S
end
Base.:(==)(h1::T,h2::T) where T<:FockSpace = (h1.name==h2.name)

"""
    Destroy <: QSym

Bosonic operator on a [`FockSpace`](@ref) representing the quantum harmonic
oscillator annihilation operator.
"""
struct Destroy{H<:HilbertSpace,S,A} <: QSym
    hilbert::H
    name::S
    aon::A
    hash::UInt
    function Destroy{H,S,A}(hilbert::H,name::S,aon::A) where {H,S,A}
        @assert has_hilbert(FockSpace,hilbert,aon)
        h = hash(Destroy, hash(hilbert, hash(name, hash(aon, zero(UInt)))))
        new(hilbert,name,aon,h)
    end
end

"""
    Create <: QSym

Bosonic operator on a [`FockSpace`](@ref) representing the quantum harmonic
oscillator creation operator.
"""
struct Create{H<:HilbertSpace,S,A} <: QSym
    hilbert::H
    name::S
    aon::A
    function Create{H,S,A}(hilbert::H,name::S,aon::A) where {H,S,A}
        @assert has_hilbert(FockSpace,hilbert,aon)
        new(hilbert,name,aon)
    end
end

for f in [:Destroy,:Create]
    @eval $(f)(hilbert::H,name::S,aon::A) where {H,S,A} = $(f){H,S,A}(hilbert,name,aon)
    @eval $(f)(hilbert::FockSpace,name) = $(f)(hilbert,name,1)
    @eval function $(f)(hilbert::ProductSpace,name)
        i = findall(x->isa(x,FockSpace),hilbert.spaces)
        if length(i)==1
            return $(f)(hilbert,name,i[1])
        else
            isempty(i) && error("Can only create $($(f)) on FockSpace! Not included in $(hilbert)")
            length(i)>1 && error("More than one FockSpace in $(hilbert)! Specify on which Hilbert space $($(f)) should be created with $($(f))(hilbert,name,i)!")
        end
    end
    @eval function embed(h::ProductSpace,op::T,aon::Int) where T<:($(f))
        check_hilbert(h.spaces[aon],op.hilbert)
        op_ = $(f)(h,op.name,aon)
        return op_
    end
    @eval function Base.hash(op::T, h::UInt) where T<:($(f))
        hash(T, hash(op.hilbert, hash(op.name, hash(op.aon, h))))
    end
end

Base.adjoint(op::Destroy) = Create(op.hilbert,op.name,acts_on(op))
Base.adjoint(op::Create) = Destroy(op.hilbert,op.name,acts_on(op))

# Commutation relation in simplification
function *(a::Destroy,b::Create)
    check_hilbert(a,b)
    aon_a = acts_on(a)
    aon_b = acts_on(b)
    if aon_a == aon_b
        return b*a + 1
    elseif aon_a < aon_b
        return QMul(1, [a,b])
    else
        return QMul(1, [b,a])
    end
end
ismergeable(::Destroy,::Create) = true
