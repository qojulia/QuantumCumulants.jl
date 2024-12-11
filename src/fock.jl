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

"""
    GroundStateProjection <: QSym

Bosonic operator on a [`FockSpace`](@ref) representing the quantum harmonic
oscillator ground state projection.

Warning: The ground state projection is not meant to be used in the cumulant
 expansion.
"""
struct GroundStateProjection{H<:HilbertSpace,S,A,M} <: QSym
    hilbert::H
    name::S
    aon::A
    metadata::M
    function GroundStateProjection{H,S,A,M}(hilbert::H, name::S, aon::A, metadata::M) where {H,S,A,M}
        @assert has_hilbert(FockSpace,hilbert,aon)
        new(hilbert,name,aon,metadata)
    end
end

for T in (:Create, :Destroy, :GroundStateProjection)
    @eval Base.isequal(a::$T, b::$T) = isequal(a.hilbert, b.hilbert) && isequal(a.name, b.name) && isequal(a.aon, b.aon)
end

for f in [:Destroy,:Create, :GroundStateProjection]
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
Base.adjoint(op::GroundStateProjection) = op

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

function *(a::GroundStateProjection,b::GroundStateProjection)
    check_hilbert(a,b)
    aon_a = acts_on(a)
    aon_b = acts_on(b)
    if aon_a == aon_b
        return a
    elseif aon_a < aon_b
        return QMul(1, [a,b])
    else
        return QMul(1, [b,a])
    end
end

for (T1,T2) in ((:Destroy, :GroundStateProjection), (:GroundStateProjection, :Create))
    @eval function *(a::$T1,b::$T2)
        check_hilbert(a,b)
        aon_a = acts_on(a)
        aon_b = acts_on(b)
        if aon_a == aon_b
            return 0
        elseif aon_a < aon_b
            return QMul(1, [a,b])
        else
            return QMul(1, [b,a])
        end
    end
end

for T1 in (:Destroy, :GroundStateProjection)
    for T2 in (:Create, :GroundStateProjection)
        # ismergeable(::$T1,::$T2) = true
        # TODO: test if faster; delete if and elseif in *-function above?
        @eval function ismergeable(a::$T1,b::$T2)
            aon_a = acts_on(a)
            aon_b = acts_on(b)
            return aon_a == aon_b
        end
    end
end
