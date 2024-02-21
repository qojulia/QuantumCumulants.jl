"""
    SpinSpace <: HilbertSpace
[`HilbertSpace`](@ref) defining a Spin space for two-level atom operators.
See also: [`Sigma`](@ref), [`Create`](@ref)
"""
struct SpinSpace{S} <: ConcreteHilbertSpace
    name::S
end
Base.:(==)(h1::T,h2::T) where T<:SpinSpace = (h1.name==h2.name)

"""
    Spin <: QSym
Spin operator on a [`SpinSpace`](@ref) representing the sigma-operators σx, σy and σz for two-level spin systems. 
The field axis represents x, y and z as 1, 2 and 3, repectively.
"""
struct Sigma{H<:HilbertSpace,S,AX<:Int,A,M} <: QSym
    hilbert::H
    name::S
    axis::AX
    aon::A
    metadata::M
    function Sigma{H,S,AX,A,M}(hilbert::H, name::S, axis::AX, aon::A, metadata::M) where {H,S,AX,A,M}
        @assert has_hilbert(SpinSpace,hilbert,aon)
        new(hilbert,name,axis,aon,metadata)
    end
end
Sigma(hilbert::H, name::S, axis::AX, aon::A; metadata::M=NO_METADATA) where {H,S,AX,A,M} = Sigma{H,S,AX,A,M}(hilbert,name,axis,aon,metadata)
Sigma(hilbert::SpinSpace, name, axis; metadata=NO_METADATA) = Sigma(hilbert, name, axis, 1; metadata)
Sigma(hilbert::SpinSpace, name, axis::Symbol; metadata=NO_METADATA) = Sigma(hilbert, name, axis, 1; metadata)
function Sigma(hilbert::ProductSpace,name,axis; metadata=NO_METADATA)
    inds = findall(x->isa(x,SpinSpace) || isa(x,ClusterSpace{<:SpinSpace}),hilbert.spaces)
    if length(inds)==1
        return Sigma(hilbert,name,axis,inds[1]; metadata)
    else
        isempty(inds) && error("Can only create Sigma on SpinSpace! Not included in $(hilbert)")
        length(inds)>1 && error("More than one SpinSpace in $(hilbert)! Specify on which Hilbert space Sigma should be created with Sigma(hilbert,name,axis,acts_on)!")
    end
end
function Sigma(hilbert::HilbertSpace,name,axis::Symbol,aon; kwargs...)
    if axis in [:x,:X]
        Sigma(hilbert,name,1,aon; kwargs...)
    elseif axis in [:y,:Y]
        Sigma(hilbert,name,2,aon; kwargs...)
    elseif axis in [:z,:Z]
        Sigma(hilbert,name,3,aon; kwargs...)
    end
end
function Sigma(hilbert::HilbertSpace,name,axis::Symbol; kwargs...)
    if axis in [:x,:X]
        Sigma(hilbert,name,1; kwargs...)
    elseif axis in [:y,:Y]
        Sigma(hilbert,name,2; kwargs...)
    elseif axis in [:z,:Z]
        Sigma(hilbert,name,3; kwargs...)
    end
end
# CallableSigma not possible, due to acts_on convenience (same number of arguments)

Base.hash(s::Sigma, h::UInt) = hash(s.hilbert, hash(s.name, hash(s.axis, hash(s.aon, h))))
Base.adjoint(s::Sigma) = s
Base.isequal(s1::Sigma,s2::Sigma) = isequal(s1.hilbert, s2.hilbert) && isequal(s1.name,s2.name) && isequal(s1.axis,s2.axis) && isequal(s1.aon,s2.aon)
Base.:(==)(s1::Sigma,s2::Sigma) = isequal(s1,s2)

ismergeable(::Sigma,::Sigma) = true
function Base.:*(si::Sigma,sj::Sigma)
    check_hilbert(si,sj)
    aon_si = acts_on(si)
    aon_sj = acts_on(sj)
    if aon_si == aon_sj
        i=si.axis
        j=sj.axis
        k = filter(x->x ∉ [i,j],[1,2,3])[1]
        sk = Sigma(si.hilbert, si.name, k, aon_si)
        return (i==j) + 1im*levicivita([i,j,k])*sk    
    elseif aon_si < aon_sj
        return QMul(1, [si,sj])
    else
        return QMul(1, [sj,si])
    end
end


### Collective Spin ###

# Note: I think that the initial state defines the size of the spin. That's somewhat cool!

"""
    CollectiveSpinSpace <: HilbertSpace
[`HilbertSpace`](@ref) defining a CollectiveSpin space for N identical two-level atom operators.
See also: [`CollectiveSigma`](@ref), [`Create`](@ref)
"""
struct CollectiveSpinSpace{S} <: ConcreteHilbertSpace
    name::S
end
Base.:(==)(h1::T,h2::T) where T<:CollectiveSpinSpace = (h1.name==h2.name)

"""
    CollectiveSpin <: QSym
CollectiveSpin operator on a [`CollectiveSpinSpace`](@ref) representing the CollectiveSigma-operators Sx, Sy and Sz for collective spin systems. 
The field axis represents x, y and z as 1, 2 and 3, repectively.
"""
struct CollectiveSigma{H<:HilbertSpace,S,AX<:Int,A,M} <: QSym
    hilbert::H
    name::S
    axis::AX
    aon::A
    metadata::M
    function CollectiveSigma{H,S,AX,A,M}(hilbert::H, name::S, axis::AX, aon::A, metadata::M) where {H,S,AX,A,M}
        @assert has_hilbert(CollectiveSpinSpace,hilbert,aon)
        new(hilbert,name,axis,aon,metadata)
    end
end
CollectiveSigma(hilbert::H, name::S, axis::AX, aon::A; metadata::M=NO_METADATA) where {H,S,AX,A,M} = CollectiveSigma{H,S,AX,A,M}(hilbert,name,axis,aon,metadata)
CollectiveSigma(hilbert::CollectiveSpinSpace, name, axis; metadata=NO_METADATA) = CollectiveSigma(hilbert, name, axis, 1; metadata)
CollectiveSigma(hilbert::CollectiveSpinSpace, name, axis::Symbol; metadata=NO_METADATA) = CollectiveSigma(hilbert, name, axis, 1; metadata)
function CollectiveSigma(hilbert::ProductSpace,name,axis; metadata=NO_METADATA)
    inds = findall(x->isa(x,CollectiveSpinSpace) || isa(x,ClusterSpace{<:CollectiveSpinSpace}),hilbert.spaces)
    if length(inds)==1
        return CollectiveSigma(hilbert,name,axis,inds[1]; metadata)
    else
        isempty(inds) && error("Can only create CollectiveSigma on CollectiveSpinSpace! Not included in $(hilbert)")
        length(inds)>1 && error("More than one CollectiveSpinSpace in $(hilbert)! Specify on which Hilbert space CollectiveSigma should be created with CollectiveSigma(hilbert,name,axis,acts_on)!")
    end
end
function CollectiveSigma(hilbert::HilbertSpace,name,axis::Symbol,aon; kwargs...)
    if axis in [:x,:X]
        CollectiveSigma(hilbert,name,1,aon; kwargs...)
    elseif axis in [:y,:Y]
        CollectiveSigma(hilbert,name,2,aon; kwargs...)
    elseif axis in [:z,:Z]
        CollectiveSigma(hilbert,name,3,aon; kwargs...)
    end
end
function CollectiveSigma(hilbert::HilbertSpace,name,axis::Symbol; kwargs...)
    if axis in [:x,:X]
        CollectiveSigma(hilbert,name,1; kwargs...)
    elseif axis in [:y,:Y]
        CollectiveSigma(hilbert,name,2; kwargs...)
    elseif axis in [:z,:Z]
        CollectiveSigma(hilbert,name,3; kwargs...)
    end
end
# CallableSigma not possible, due to acts_on convenience (same number of arguments)

Base.hash(s::CollectiveSigma, h::UInt) = hash(s.hilbert, hash(s.name, hash(s.axis, hash(s.aon, h))))
Base.adjoint(s::CollectiveSigma) = s
Base.isequal(s1::CollectiveSigma,s2::CollectiveSigma) = isequal(s1.hilbert, s2.hilbert) && isequal(s1.name,s2.name) && isequal(s1.axis,s2.axis) && isequal(s1.aon,s2.aon)
Base.:(==)(s1::CollectiveSigma,s2::CollectiveSigma) = isequal(s1,s2)

# ismergeable(::CollectiveSigma,::CollectiveSigma) = true
ismergeable(si::CollectiveSigma,sj::CollectiveSigma) = (acts_on(si) == acts_on(sj) && si.axis > sj.axis)

function Base.:*(si::CollectiveSigma,sj::CollectiveSigma)
    check_hilbert(si,sj)
    aon_si = acts_on(si)
    aon_sj = acts_on(sj)
    if aon_si == aon_sj
        i=si.axis
        j=sj.axis
        
        if i==2 && j==1 #SySx
            SiSj = sj*si - 1im*CollectiveSigma(si.hilbert, si.name, 3, aon_si)
        elseif i==3 && j==2 #SzSy
            SiSj = sj*si - 1im*CollectiveSigma(si.hilbert, si.name, 1, aon_si)
        elseif i==3 && j==1 #SySx
            SiSj = sj*si + 1im*CollectiveSigma(si.hilbert, si.name, 2, aon_si)
        else # SxSx, SySy, SzSz, SxSy, Sx,Sz, SySz
            SiSj = QMul(1, [si,sj])
        end
        return SiSj
    elseif aon_si < aon_sj
        return QMul(1, [si,sj])
    else
        return QMul(1, [sj,si])
    end
end
