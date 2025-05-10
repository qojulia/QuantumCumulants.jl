"""
    PauliSpace <: HilbertSpace
[`HilbertSpace`](@ref) defining a Spin space for two-level atom Pauli operators.
See also: [`Pauli`](@ref), [`Create`](@ref)
"""
struct PauliSpace{S} <: ConcreteHilbertSpace
    name::S
end
Base.:(==)(h1::T, h2::T) where {T<:PauliSpace} = (h1.name==h2.name)

"""
    Pauli <: QSym
Pauli operator on a [`PauliSpace`](@ref) representing the Pauli operators σx, σy and σz for two-level spin systems.
The field axis represents x, y and z as 1, 2 and 3, repectively. The used rewriting rule is σj⋅σk → δjk + i⋅ϵjkl⋅σl.

Examples
=======
```
julia> h = PauliSpace("Spin-1/2")
ℋ(Spin-1/2)

julia> σx = Pauli(h,:σ,1)
σx
```
"""
struct Pauli{H<:HilbertSpace,S,AX<:Int,A,M} <: QSym
    hilbert::H
    name::S
    axis::AX
    aon::A
    metadata::M
    function Pauli{H,S,AX,A,M}(
        hilbert::H, name::S, axis::AX, aon::A, metadata::M
    ) where {H,S,AX,A,M}
        @assert has_hilbert(PauliSpace, hilbert, aon)
        new(hilbert, name, axis, aon, metadata)
    end
end
function Pauli(
    hilbert::H, name::S, axis::AX, aon::A; metadata::M=NO_METADATA
) where {H,S,AX,A,M}
    Pauli{H,S,AX,A,M}(hilbert, name, axis, aon, metadata)
end
function Pauli(hilbert::PauliSpace, name, axis; metadata=NO_METADATA)
    Pauli(hilbert, name, axis, 1; metadata)
end
function Pauli(hilbert::PauliSpace, name, axis::Symbol; metadata=NO_METADATA)
    Pauli(hilbert, name, axis, 1; metadata)
end
function Pauli(hilbert::ProductSpace, name, axis; metadata=NO_METADATA)
    inds = findall(
        x->isa(x, PauliSpace) || isa(x, ClusterSpace{<:PauliSpace}), hilbert.spaces
    )
    if length(inds)==1
        return Pauli(hilbert, name, axis, inds[1]; metadata)
    else
        isempty(inds) &&
            error("Can only create Pauli on PauliSpace! Not included in $(hilbert)")
        length(inds)>1 && error(
            "More than one PauliSpace in $(hilbert)! Specify on which Hilbert space Pauli should be created with Pauli(hilbert,name,axis,acts_on)!",
        )
    end
end
function Pauli(hilbert::HilbertSpace, name, axis::Symbol, aon; kwargs...)
    if axis in [:x, :X]
        Pauli(hilbert, name, 1, aon; kwargs...)
    elseif axis in [:y, :Y]
        Pauli(hilbert, name, 2, aon; kwargs...)
    elseif axis in [:z, :Z]
        Pauli(hilbert, name, 3, aon; kwargs...)
    end
end
function Pauli(hilbert::HilbertSpace, name, axis::Symbol; kwargs...)
    if axis in [:x, :X]
        Pauli(hilbert, name, 1; kwargs...)
    elseif axis in [:y, :Y]
        Pauli(hilbert, name, 2; kwargs...)
    elseif axis in [:z, :Z]
        Pauli(hilbert, name, 3; kwargs...)
    end
end
# CallablePauli not possible, due to acts_on convenience (same number of arguments)

Base.hash(s::Pauli, h::UInt) = hash(s.hilbert, hash(s.name, hash(s.axis, hash(s.aon, h))))
Base.adjoint(s::Pauli) = s
function Base.isequal(s1::Pauli, s2::Pauli)
    isequal(s1.hilbert, s2.hilbert) &&
        isequal(s1.name, s2.name) &&
        isequal(s1.axis, s2.axis) &&
        isequal(s1.aon, s2.aon)
end
Base.:(==)(s1::Pauli, s2::Pauli) = isequal(s1, s2)

ismergeable(::Pauli, ::Pauli) = true
function Base.:*(si::Pauli, sj::Pauli)
    check_hilbert(si, sj)
    aon_si = acts_on(si)
    aon_sj = acts_on(sj)
    if aon_si == aon_sj
        i=si.axis
        j=sj.axis
        k = filter(x->x ∉ [i, j], [1, 2, 3])[1]
        sk = Pauli(si.hilbert, si.name, k, aon_si)
        return (i==j) + 1im*Combinatorics.levicivita([i, j, k])*sk
    elseif aon_si < aon_sj
        return QMul(1, [si, sj])
    else
        return QMul(1, [sj, si])
    end
end

### Spin-N/2 ###
"""
    SpinSpace <: HilbertSpace
[`HilbertSpace`](@ref) defining a Spin space for N > 1 identical two-level atom operators.
See also: [`Spin`](@ref), [`Create`](@ref)
"""
struct SpinSpace{S} <: ConcreteHilbertSpace
    name::S
end
Base.:(==)(h1::T, h2::T) where {T<:SpinSpace} = (h1.name==h2.name)

"""
    Spin <: QSym
Spin operator on a [`SpinSpace`](@ref) representing the Spin-operators Sx, Sy and Sz for collective spin systems.
The field axis represents x, y and z as 1, 2 and 3, repectively. The operators follow the rules for angular momentum operators:
[Sj,Sk] = i⋅∑ϵjkl⋅Sl

Examples
=======
```
julia> h = SpinSpace("Spin-N/2")
ℋ(Spin-N/2)

julia> Sx = Spin(h,:S,1)
Sx
```
"""
struct Spin{H<:HilbertSpace,S,AX<:Int,A,M} <: QSym
    hilbert::H
    name::S
    axis::AX
    aon::A
    metadata::M
    function Spin{H,S,AX,A,M}(
        hilbert::H, name::S, axis::AX, aon::A, metadata::M
    ) where {H,S,AX,A,M}
        @assert has_hilbert(SpinSpace, hilbert, aon)
        new(hilbert, name, axis, aon, metadata)
    end
end
function Spin(
    hilbert::H, name::S, axis::AX, aon::A; metadata::M=NO_METADATA
) where {H,S,AX,A,M}
    Spin{H,S,AX,A,M}(hilbert, name, axis, aon, metadata)
end
function Spin(hilbert::SpinSpace, name, axis; metadata=NO_METADATA)
    Spin(hilbert, name, axis, 1; metadata)
end
function Spin(hilbert::SpinSpace, name, axis::Symbol; metadata=NO_METADATA)
    Spin(hilbert, name, axis, 1; metadata)
end
function Spin(hilbert::ProductSpace, name, axis; metadata=NO_METADATA)
    inds = findall(
        x->isa(x, SpinSpace) || isa(x, ClusterSpace{<:SpinSpace}), hilbert.spaces
    )
    if length(inds)==1
        return Spin(hilbert, name, axis, inds[1]; metadata)
    else
        isempty(inds) &&
            error("Can only create Spin on SpinSpace! Not included in $(hilbert)")
        length(inds)>1 && error(
            "More than one SpinSpace in $(hilbert)! Specify on which Hilbert space Spin should be created with Spin(hilbert,name,axis,acts_on)!",
        )
    end
end
function Spin(hilbert::HilbertSpace, name, axis::Symbol, aon; kwargs...)
    if axis in [:x, :X]
        Spin(hilbert, name, 1, aon; kwargs...)
    elseif axis in [:y, :Y]
        Spin(hilbert, name, 2, aon; kwargs...)
    elseif axis in [:z, :Z]
        Spin(hilbert, name, 3, aon; kwargs...)
    end
end
function Spin(hilbert::HilbertSpace, name, axis::Symbol; kwargs...)
    if axis in [:x, :X]
        Spin(hilbert, name, 1; kwargs...)
    elseif axis in [:y, :Y]
        Spin(hilbert, name, 2; kwargs...)
    elseif axis in [:z, :Z]
        Spin(hilbert, name, 3; kwargs...)
    end
end

Base.hash(s::Spin, h::UInt) = hash(s.hilbert, hash(s.name, hash(s.axis, hash(s.aon, h))))
Base.adjoint(s::Spin) = s
function Base.isequal(s1::Spin, s2::Spin)
    isequal(s1.hilbert, s2.hilbert) &&
        isequal(s1.name, s2.name) &&
        isequal(s1.axis, s2.axis) &&
        isequal(s1.aon, s2.aon)
end
Base.:(==)(s1::Spin, s2::Spin) = isequal(s1, s2)

# ismergeable(::Spin,::Spin) = true
ismergeable(si::Spin, sj::Spin) = (acts_on(si) == acts_on(sj) && si.axis > sj.axis)

function Base.:*(si::Spin, sj::Spin)
    check_hilbert(si, sj)
    aon_si = acts_on(si)
    aon_sj = acts_on(sj)
    if aon_si == aon_sj
        i=si.axis
        j=sj.axis

        if i==2 && j==1 #SySx
            SiSj = sj*si - 1im*Spin(si.hilbert, si.name, 3, aon_si)
        elseif i==3 && j==2 #SzSy
            SiSj = sj*si - 1im*Spin(si.hilbert, si.name, 1, aon_si)
        elseif i==3 && j==1 #SySx
            SiSj = sj*si + 1im*Spin(si.hilbert, si.name, 2, aon_si)
        else # SxSx, SySy, SzSz, SxSy, Sx,Sz, SySz
            SiSj = QMul(1, [si, sj])
        end
        return SiSj
    elseif aon_si < aon_sj
        return QMul(1, [si, sj])
    else
        return QMul(1, [sj, si])
    end
end
