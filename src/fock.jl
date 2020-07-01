"""
    FockSpace <: HilbertSpace

[`HilbertSpace`](@ref) defining a Fock space for bosonic operators.
See also: [`Destroy`](@ref), [`Create`](@ref)
"""
struct FockSpace{S} <: HilbertSpace
    name::S
end
Base.isequal(h1::T,h2::T) where T<:FockSpace = isequal(h1.name, h2.name)

"""
    Destroy <: BasicOperator

Bosonic operator on a [`FockSpace`](@ref) representing the quantum harmonic
oscillator annihilation operator.
"""
struct Destroy{H<:HilbertSpace,S,A,IND} <: BasicOperator
    hilbert::H
    name::S
    aon::A
    index::IND
    function Destroy{H,S,A,IND}(hilbert::H,name::S,aon::A,index::IND) where {H,S,A,IND}
        @assert has_hilbert(FockSpace,hilbert,aon)
        op = new(hilbert,name,aon,index)
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
struct Create{H<:HilbertSpace,S,A,IND} <: BasicOperator
    hilbert::H
    name::S
    aon::A
    index::IND
    function Create{H,S,A,IND}(hilbert::H,name::S,aon::A,index::IND) where {H,S,A,IND}
        @assert has_hilbert(FockSpace,hilbert,aon)
        op = new(hilbert,name,aon,index)
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
    @eval $(f)(hilbert::H,name::S,aon::A;index::IND=default_index()) where {H,S,A,IND} = $(f){H,S,A,IND}(hilbert,name,aon,index)
    @eval $(f)(hilbert::FockSpace,name;kwargs...) = $(f)(hilbert,name,1;kwargs...)
    @eval function $(f)(hilbert::ProductSpace,name;kwargs...)
        i = findall(x->isa(x,FockSpace),hilbert.spaces)
        if length(i)==1
            return $(f)(hilbert,name,i[1];kwargs...)
        else
            isempty(i) && error("Can only create $($(f)) on FockSpace! Not included in $(hilbert)")
            length(i)>1 && error("More than one FockSpace in $(hilbert)! Specify on which Hilbert space $($(f)) should be created with $($(f))(hilbert,name,i)!")
        end
    end
    @eval function embed(h::ProductSpace,op::T,aon::Int;kwargs...) where T<:($(f))
        check_hilbert(h.spaces[aon],op.hilbert)
        op_ = $(f)(h,op.name,aon;kwargs...)
        return op_
    end
    @eval function Base.hash(op::T, h::UInt) where T<:($(f))
        hash(op.hilbert, hash(op.name, hash(op.aon, hash(op.index, hash($(f), h)))))
    end
end

Base.adjoint(op::Destroy) = Create(op.hilbert,op.name,acts_on(op);index=get_index(op))
Base.adjoint(op::Create) = Destroy(op.hilbert,op.name,acts_on(op);index=get_index(op))

# Commutation relation in simplification
function commute_bosonic(a,b)
    if acts_on(a)==acts_on(b)
        idx1 = _to_symbolic(get_index(a))
        idx2 = _to_symbolic(get_index(b))
        δ = idx1==idx2
        return b*a + δ
    else
        return nothing
    end
end
