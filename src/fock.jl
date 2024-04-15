"""
    FockSpace <: HilbertSpace

[`HilbertSpace`](@ref) defining a Fock space for bosonic operators.
See also: [`Destroy`](@ref), [`Create`](@ref)
"""
struct FockSpace{S} <: ConcreteHilbertSpace
    name::S
end
Base.:(==)(h1::T,h2::T) where T<:FockSpace = (h1.name==h2.name)

"""
    Destroy <: QSym

Bosonic operator on a [`FockSpace`](@ref) representing the quantum harmonic
oscillator annihilation operator.
"""
struct Destroy{H<:HilbertSpace,S,A,M} <: QSym
    hilbert::H
    name::S
    aon::A
    metadata::M
    function Destroy{H,S,A,M}(hilbert::H, name::S, aon::A, metadata::M) where {H,S,A,M}
        @assert has_hilbert(FockSpace,hilbert,aon)
        new(hilbert,name,aon,metadata)
    end
end

"""
    Create <: QSym

Bosonic operator on a [`FockSpace`](@ref) representing the quantum harmonic
oscillator creation operator.
"""
struct Create{H<:HilbertSpace,S,A,M} <: QSym
    hilbert::H
    name::S
    aon::A
    metadata::M
    function Create{H,S,A,M}(hilbert::H, name::S, aon::A, metadata::M) where {H,S,A,M}
        @assert has_hilbert(FockSpace,hilbert,aon)
        new(hilbert,name,aon,metadata)
    end
end

for T in (:Create, :Destroy)
    @eval Base.isequal(a::$T, b::$T) = isequal(a.hilbert, b.hilbert) && isequal(a.name, b.name) && isequal(a.aon, b.aon)
end

for f in [:Destroy,:Create]
    @eval $(f)(hilbert::H, name::S, aon::A; metadata::M=NO_METADATA) where {H,S,A,M} = $(f){H,S,A,M}(hilbert,name,aon,metadata)
    @eval $(f)(hilbert::FockSpace, name; metadata=NO_METADATA) = $(f)(hilbert,name,1; metadata)
    @eval function $(f)(hilbert::H, name::S, aon::A; metadata::M=NO_METADATA) where {H<:ProductSpace,S,A<:Int,M}
        if hilbert.spaces[aon] isa ClusterSpace
            hilbert.spaces[get_i(aon)].op_name[] = name
            op = $(f)(hilbert.spaces[aon].original_space, name; metadata)
            return _cluster(hilbert, op, aon) #todo OK?
        else
            return $(f){H,S,A,M}(hilbert,name,aon,metadata)
        end
    end
    @eval function $(f)(hilbert::ProductSpace, name; metadata=NO_METADATA)
        i = findall(x->isa(x,FockSpace) || isa(x,ClusterSpace{<:FockSpace}),hilbert.spaces)
        if length(i)==1
            return $(f)(hilbert, name, i[1]; metadata)
        else
            isempty(i) && error("Can only create $($(f)) on FockSpace! Not included in $(hilbert)")
            length(i)>1 && error("More than one FockSpace in $(hilbert)! Specify on which Hilbert space $($(f)) should be created with $($(f))(hilbert,name,i)!")
        end
    end
    @eval function Base.hash(op::T, h::UInt) where T<:($(f))
        hash($(f), hash(op.hilbert, hash(op.name, hash(op.aon, h))))
    end
end

Base.adjoint(op::Destroy) = Create(op.hilbert,op.name,acts_on(op); op.metadata)
Base.adjoint(op::Create) = Destroy(op.hilbert,op.name,acts_on(op); op.metadata)

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
# ismergeable(::Destroy,::Create) = true

# TODO: test if faster; delete if and elseif in *-function above?
function ismergeable(a::Destroy,b::Create)
    aon_a = acts_on(a)
    aon_b = acts_on(b)
    return aon_a == aon_b
end