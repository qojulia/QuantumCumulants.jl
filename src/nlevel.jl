"""
    NLevelSpace <: HilbertSpace
    NLevelSpace(name::Symbol,levels,GS=1)

Define a [`HilbertSpace`](@ref) for an object consisting of `N` discrete energy
levels. The given `levels` must be an integer specifying the number of levels,
or an iterable collection of levels. The argument `GS` specifies which state
should be treated as ground state and is rewritten using population conservation
during simplification.
See also: [`Transition`](@ref)

Examples:
========
```
julia> ha = NLevelSpace(:a,3)
ℋ(a)

julia> ha = NLevelSpace(:a,(:g,:e))
ℋ(a)
```
"""
struct NLevelSpace{S,L,G} <: HilbertSpace
    name::S
    levels::L
    GS::G
end
NLevelSpace(name,N::Int,GS) = NLevelSpace(name,1:N,GS)
NLevelSpace(name,N::Int) = NLevelSpace(name,1:N,1)
NLevelSpace(name,levels) = NLevelSpace(name,levels,levels[1])
Base.:(==)(h1::T,h2::T) where T<:NLevelSpace = (h1.name==h2.name && h1.levels==h2.levels && h1.GS==h2.GS)
Base.hash(n::NLevelSpace, h::UInt) = hash(n.name, hash(n.levels, hash(n.GS, h)))

levels(h::NLevelSpace) = h.levels
levels(h::NLevelSpace,aon) = levels(h)
levels(h::ProductSpace,aon) = levels(h.spaces[aon])
ground_state(h::NLevelSpace) = h.GS
ground_state(h::NLevelSpace,aon) = h.GS
ground_state(h::ProductSpace,aon) = ground_state(h.spaces[aon])

"""
    Transition <: QSym
    Transition(h::NLevelSpace,name::Symbol,i,j)

Fundamental operator defining a transition from level `j` to level `i` on a
[`NLevelSpace`](@ref). The notation corresponds to Dirac notation, i.e. the
above is equivalent to `|i⟩⟨j|`.

Examples
=======
```
julia> ha = NLevelSpace(:a,(:g,:e))
ℋ(a)

julia> σ = Transition(ha,:σ,:g,:e)
σge
```
"""
struct Transition{H,S,I,A} <: QSym
    hilbert::H
    name::S
    i::I
    j::I
    aon::A
    function Transition{H,S,I,A}(hilbert::H,name::S,i::I,j::I,aon::A) where {H,S,I,A}
        @assert has_hilbert(NLevelSpace,hilbert,aon)
        @assert i∈levels(hilbert,aon) && j∈levels(hilbert,aon)
        new(hilbert,name,i,j,aon)
    end
end
Transition(hilbert::H,name::S,i::I,j::I,aon::A) where {H,S,I,A} = Transition{H,S,I,A}(hilbert,name,i,j,aon)
Transition(hilbert::NLevelSpace,name,i,j) = Transition(hilbert,name,i,j,1)
function Transition(hilbert::ProductSpace,name,i,j)
    inds = findall(x->isa(x,NLevelSpace),hilbert.spaces)
    if length(inds)==1
        return Transition(hilbert,name,i,j,inds[1])
    else
        isempty(inds) && error("Can only create Transition on NLevelSpace! Not included in $(hilbert)")
        length(inds)>1 && error("More than one NLevelSpace in $(hilbert)! Specify on which Hilbert space Transition should be created with Transition(hilbert,name,i,j,acts_on)!")
    end
end

function Base.isless(a::Transition, b::Transition)
    if a.name == b.name
        if a.i == b.i
            return a.j < b.j
        else
            return a.i < b.i
        end
    else
        return a.name < b.name
    end
end

function embed(h::ProductSpace,op::T,aon::Int) where T<:Transition
    check_hilbert(h.spaces[aon],op.hilbert)
    op_ = Transition(h,op.name,op.i,op.j,aon)
    return op_
end
levels(t::Transition,args...) = levels(t.hilbert,args...)
ground_state(t::Transition,args...) = ground_state(t.hilbert,args...)

Base.adjoint(t::Transition) = Transition(t.hilbert,t.name,t.j,t.i,acts_on(t))
Base.isequal(t1::Transition,t2::Transition) = isequal(t1.hilbert, t2.hilbert) && isequal(t1.name,t2.name) && isequal(t1.i,t2.i) && isequal(t1.j,t2.j) && isequal(t1.aon,t2.aon)
Base.hash(t::Transition, h::UInt) = hash(t.hilbert, hash(t.name, hash(t.i, hash(t.j, hash(t.aon, h)))))

"""
    CallableTransition

A [`Transition`](@ref) where no levels have been specified. This type is callable
to allow for easy construction of concrete [`Transition`](@ref) instances.

Examples
========
```
julia> h = NLevelSpace(:atom, (:g,:e))
ℋ(atom)

julia> σ = Transition(h,:σ)
σ

julia> σ(:g,:e)
σge
```
"""
struct CallableTransition{H,S,A}
    hilbert::H
    name::S
    aon::A
end

acts_on(c::CallableTransition) = c.aon

function (c::CallableTransition)(i, j)
    Transition(c.hilbert, c.name, i, j, c.aon)
end

Transition(hilbert, name, aon) = CallableTransition(hilbert, name, aon)
Transition(hilbert::NLevelSpace, name) = CallableTransition(hilbert, name, 1)
function Transition(hilbert::ProductSpace,name)
    inds = findall(x->isa(x,NLevelSpace),hilbert.spaces)
    if length(inds)==1
        return CallableTransition(hilbert,name,inds[1])
    else
        isempty(inds) && error("Can only create Transition on NLevelSpace! Not included in $(hilbert)")
        length(inds)>1 && error("More than one NLevelSpace in $(hilbert)! Specify on which Hilbert space Transition should be created with Transition(hilbert,name,i,j,acts_on)!")
    end
end

# Simplification
function merge_transitions(σ1::Transition, σ2::Transition)
    if σ1.j == σ2.i
        return Transition(σ1.hilbert,σ1.name,σ1.i,σ2.j,σ1.aon)
    else
        return 0
    end
end
function rewrite_gs(σ::Transition)
    h = σ.hilbert
    aon = acts_on(σ)
    gs = ground_state(h,aon)
    i,j = σ.i, σ.j
    if i==j==gs
        args = Any[1]
        for k in levels(h,aon)
            if k != i
                t_ = QTerm(*, [-1, Transition(h, σ.name, k, k, aon)])
                push!(args, t_)
            end
        end
        return QTerm(+, args)
    else
        return nothing
    end
end
