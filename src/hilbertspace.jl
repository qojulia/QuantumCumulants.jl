abstract type HilbertSpace end
Base.isequal(h1::HilbertSpace,h2::HilbertSpace) = false
Base.hash(h::HilbertSpace, i::UInt) = hash(h.name, i)

"""
    ProductSpace <: HilbertSpace

Stores a composite [`HilbertSpace`](@ref) consisting of multiple subspaces.
Generally created by computing the tensor product [`⊗`](@ref) of subspaces.
"""
struct ProductSpace{S} <: HilbertSpace
    spaces::S
end
Base.isequal(h1::T,h2::T) where T<:ProductSpace = isequal(h1.spaces, h2.spaces)
Base.hash(p::ProductSpace, h::UInt) = hash(p.spaces, h)

"""
    ⊗(spaces::HilbertSpace...)

Create a [`ProductSpace`](@ref) consisting of multiple subspaces.

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

Base.isless(h1::HilbertSpace,h2::HilbertSpace) = isless(h1.name,h2.name)
Base.isless(h1::ProductSpace,h2::ProductSpace) = isless(h1.spaces,h2.spaces)

function Base.copy(h::T) where T<:HilbertSpace
    fields = [getfield(h, n) for n in fieldnames(T)]
    return T(fields...)
end

has_hilbert(::Type{T},::T,args...) where T<:HilbertSpace = true
has_hilbert(T::Type{<:HilbertSpace},h::ProductSpace,aon) = has_hilbert(T,h.spaces[aon])
has_hilbert(::Type{<:HilbertSpace},::HilbertSpace,args...) = false
