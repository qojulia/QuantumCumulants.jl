"""
    HilbertSpace

Abstract type for representing Hilbert spaces.
"""
abstract type HilbertSpace end
Base.:(==)(h1::HilbertSpace,h2::HilbertSpace) = false
Base.hash(h::T, i::UInt) where T<:HilbertSpace = hash(T, hash(h.name, i))

abstract type ConcreteHilbertSpace <: HilbertSpace end

"""
    ProductSpace <: HilbertSpace

Stores a composite [`HilbertSpace`](@ref) consisting of multiple subspaces.
Generally created by computing the tensor product [`⊗`](@ref) of subspaces.
"""
struct ProductSpace{S} <: HilbertSpace
    spaces::S
    hash::UInt
    function ProductSpace(spaces::S) where S
        h = hash(ProductSpace, hash(spaces, zero(UInt)))
        new{S}(spaces, h)
    end
end
Base.:(==)(h1::T,h2::T) where T<:ProductSpace = isequal(h1.hash, h2.hash)
function Base.hash(p::ProductSpace, h::UInt)
    iszero(h) && return p.hash
    return hash(hash(p, zero(UInt)), h)
end

"""
    ⊗(spaces::HilbertSpace...)

Create a [`ProductSpace`](@ref) consisting of multiple subspaces.
Unicode `\\otimes<tab>` alias of [`tensor`](@ref)

Examples
=======
```
julia> hf = FockSpace(:f)
ℋ(f)

julia> ha = NLevelSpace(:a,2)
ℋ(a)

julia> h = hf⊗ha
ℋ(f) ⊗ ℋ(a)
```
"""
⊗(a::HilbertSpace,b::HilbertSpace) = ProductSpace([a,b])
⊗(a::HilbertSpace,b::ProductSpace) = ProductSpace([a;b.spaces])
⊗(a::ProductSpace,b::HilbertSpace) = ProductSpace([a.spaces;b])
⊗(a::ProductSpace,b::ProductSpace) = ProductSpace([a.spaces;b.spaces])
⊗(a::HilbertSpace, b::HilbertSpace, c::HilbertSpace...) = ⊗(a⊗b,c...)
⊗(a::HilbertSpace) = a

"""
    tensor(spaces::HilbertSpace...)

Create a [`ProductSpace`](@ref) consisting of multiple subspaces.
See also [`⊗`](@ref).
"""
tensor(args...) = ⊗(args...)

"""
    ClusterSpace <: HilbertSpace
    ClusterSpace(original_space,N,order)

A Hilbert space representing `N` identical copies of another Hilbert space, with
correlations up to a specified `order`.
"""
struct ClusterSpace{H<:ConcreteHilbertSpace,NType,M<:Integer} <: HilbertSpace
    original_space::H
    N::NType
    order::M
end
Base.:(==)(h1::T,h2::T) where T<:ClusterSpace = (h1.original_space==h2.original_space && isequal(h1.N,h2.N) && h1.order==h2.order)
Base.hash(c::ClusterSpace, h::UInt) = hash(c.original_space, hash(c.N, hash(c.order, h)))

"""
    ClusterAon(i,j)

When an operator acts on the Hilbert space `i` which is a [`ClusterSpace`](@ref),
the index `j` denotes which copy of the Hilbert space the operator acts on.
"""
struct ClusterAon{T<:Integer}
    i::T
    j::T
end

Base.isless(h1::HilbertSpace,h2::HilbertSpace) = isless(h1.name,h2.name)
Base.isless(h1::ProductSpace,h2::ProductSpace) = isless(h1.spaces,h2.spaces)

function Base.copy(h::T) where T<:HilbertSpace
    fields = [getfield(h, n) for n in fieldnames(T)]
    return T(fields...)
end

has_hilbert(::Type{T},::T,args...) where T<:HilbertSpace = true
has_hilbert(T::Type{<:HilbertSpace},h::ProductSpace,aon) = has_hilbert(T,h.spaces[aon])
has_hilbert(::Type{<:HilbertSpace},::HilbertSpace,args...) = false
has_hilbert(::Type{T},::ClusterSpace{<:T},args...) where T<:HilbertSpace = true
