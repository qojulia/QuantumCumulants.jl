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
function Transition(hilbert::H,name::S,i::I,j::I,aon::A) where {H,S,I,A}
    gs = ground_state(hilbert, aon)
    if isequal(i,j) && isequal(i,gs)
        args = Any[1]
        for k∈levels(hilbert, aon)
            isequal(k,gs) && continue
            push!(args, QMul(-1, [Transition{H,S,I,A}(hilbert,name,k,k,aon)]))
        end
        return QAdd(args)
    else
        return Transition{H,S,I,A}(hilbert,name,i,j,aon)
    end
end
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
function *(a::Transition,b::Transition)
    check_hilbert(a, b)
    aon_a = acts_on(a)
    aon_b = acts_on(b)
    if aon_a == aon_b
        if isequal(a.j, b.i)
            return Transition(a.hilbert, a.name, a.i, b.j, a.aon)
        else
            return 0
        end
    elseif aon_a < aon_b
        return QMul(1, [a,b])
    else
        return QMul(1, [b,a])
    end
end
ismergeable(::Transition,::Transition) = true
