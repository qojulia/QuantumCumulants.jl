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
    Destroy <: BasicOperator

Bosonic operator on a [`FockSpace`](@ref) representing the quantum harmonic
oscillator annihilation operator.
"""
struct Destroy{H<:HilbertSpace,S,A} <: BasicOperator
    hilbert::H
    name::S
    aon::A
    function Destroy{H,S,A}(hilbert::H,name::S,aon::A) where {H,S,A}
        @assert has_hilbert(FockSpace,hilbert,aon)
        op = new(hilbert,name,aon)
        if !haskey(OPERATORS_TO_SYMS, op)
            sym = SymbolicUtils.Sym{Destroy}(gensym(:Destroy))
            OPERATORS_TO_SYMS[op] = sym
            SYMS_TO_OPERATORS[sym] = op
        end
        return op
    end
end
isdestroy(a::SymbolicUtils.Sym{T}) where {T<:Destroy} = true

"""
    Create <: BasicOperator

Bosonic operator on a [`FockSpace`](@ref) representing the quantum harmonic
oscillator creation operator.
"""
struct Create{H<:HilbertSpace,S,A} <: BasicOperator
    hilbert::H
    name::S
    aon::A
    function Create{H,S,A}(hilbert::H,name::S,aon::A) where {H,S,A}
        @assert has_hilbert(FockSpace,hilbert,aon)
        op = new(hilbert,name,aon)
        if !haskey(OPERATORS_TO_SYMS, op)
            sym = SymbolicUtils.Sym{Create}(gensym(:Create))
            OPERATORS_TO_SYMS[op] = sym
            SYMS_TO_OPERATORS[sym] = op
        end
        return op
    end
end
iscreate(a::SymbolicUtils.Sym{T}) where {T<:Create} = true

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
        hash(op.hilbert, hash(op.name, hash(op.aon, hash($(f), h))))
    end
end

Base.adjoint(op::Destroy) = Create(op.hilbert,op.name,acts_on(op))
Base.adjoint(op::Create) = Destroy(op.hilbert,op.name,acts_on(op))

# Commutation relation in simplification
function commute_bosonic(a,b)
    if acts_on(a)==acts_on(b)
        return b*a + one(a)
    else
        return nothing
    end
end
