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
struct NLevelSpace{S,L,G,N} <: HilbertSpace
    name::S
    levels::L
    GS::G
    n::N
    function NLevelSpace{S,L,G,N}(name::S,levels::L,GS::G,n::N) where {S,L,G,N}
        r1 = SymbolicUtils.@rule(*(~~x::has_consecutive(istransition)) => merge_transitions(*, ~~x))
        r2 = SymbolicUtils.@rule(~x::istransition => rewrite_gs(~x))
        (r1 ∈ COMMUTATOR_RULES.rules) || push!(COMMUTATOR_RULES.rules, r1)
        (r2 ∈ COMMUTATOR_RULES.rules) || push!(COMMUTATOR_RULES.rules, r2)
        new(name,levels,GS,n)
    end
end
NLevelSpace(name::S,levels::L,GS::G,n::N) where {S,L,G,N} = NLevelSpace{S,L,G,N}(name,levels,GS,n)
NLevelSpace(name,N::Int,GS;n=1) = NLevelSpace(name,1:N,GS,n)
NLevelSpace(name,N::Int;n=1) = NLevelSpace(name,1:N,1,n)
NLevelSpace(name,levels;n=1) = NLevelSpace(name,levels,levels[1],n)
Base.:(==)(h1::T,h2::T) where T<:NLevelSpace = (h1.name==h2.name && h1.levels==h2.levels && h1.GS==h2.GS && h1.n==h2.n)

levels(h::NLevelSpace) = h.levels
levels(h::NLevelSpace,aon) = levels(h)
levels(h::ProductSpace,aon) = levels(h.spaces[aon])
ground_state(h::NLevelSpace) = h.GS
ground_state(h::NLevelSpace,aon) = h.GS
ground_state(h::ProductSpace,aon) = ground_state(h.spaces[aon])

"""
    Transition <: BasicOperator
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
struct Transition{H,S,I,A,IND} <: BasicOperator
    hilbert::H
    name::S
    i::I
    j::I
    aon::A
    index::IND
    function Transition{H,S,I,A,IND}(hilbert::H,name::S,i::I,j::I,aon::A,index::IND) where {H,S,I,A,IND}
        @assert i∈levels(hilbert,aon) && j∈levels(hilbert,aon)
        @assert 0 < index <= length(hilbert,aon)
        return new(hilbert,name,i,j,aon,index)
    end
end
Transition(hilbert::H,name::S,i::I,j::I,aon::A,index::IND=1) where {H,S,I,A,IND} = Transition{H,S,I,A,IND}(hilbert,name,i,j,aon,index)
Transition(hilbert::NLevelSpace,name,i,j) = Transition(hilbert,name,i,j,1,1)
function Transition(hilbert::ProductSpace,name,i,j;index=1)
    inds = findall(x->isa(x,NLevelSpace),hilbert.spaces)
    if length(inds)==1
        return Transition(hilbert,name,i,j,inds[1],index)
    else
        isempty(inds) && error("Can only create Transition on NLevelSpace! Not included in $(hilbert)")
        length(inds)>1 && error("More than one NLevelSpace in $(hilbert)! Specify on which Hilbert space Transition should be created with Transition(hilbert,name,i,j,acts_on)!")
    end
end

function embed(h::ProductSpace,op::T,aon::Int) where T<:Transition
    check_hilbert(h.spaces[aon],op.hilbert)
    op_ = Transition(h,op.name,op.i,op.j,aon)
    return op_
end
levels(t::Transition,args...) = levels(t.hilbert,args...)
ground_state(t::Transition,args...) = ground_state(t.hilbert,args...)

function _to_symbolic(op::T) where T<:Transition
    aon = acts_on(op) - 1 + get_index(op)
    sym = SymbolicUtils.term(Transition, op.hilbert, op.name, op.i, op.j, aon, op.index; type=Transition)
    return sym
end
function _to_qumulants(t::SymbolicUtils.Term{T}) where T<:Transition
    aon = acts_on(t) + 1 - get_index(t)
    return t.f(t.arguments...)
end
function acts_on(t::SymbolicUtils.Term{T}) where T<:Transition
    return t.arguments[5]
end

Base.adjoint(t::Transition) = Transition(t.hilbert,t.name,t.j,t.i,acts_on(t),get_index(t))
Base.:(==)(t1::Transition,t2::Transition) = (t1.hilbert==t2.hilbert && t1.name==t2.name && t1.i==t2.i && t1.j==t2.j && t1.aon==t2.aon && t1.index==t2.index)

# Simplification
istransition(x) = false
istransition(x::Union{T,SymbolicUtils.Term{T}}) where T<:Transition = true

function merge_transitions(f::Function, args)
    merged = Any[]
    i = 1
    while i <= length(args)
       if istransition(args[i])&&(i<length(args))&&istransition(args[i+1])&&(acts_on(args[i])==acts_on(args[i+1]))
           push!(merged, merge_transitions(args[i],args[i+1]))
           i += 2
       else
           push!(merged, args[i])
           i += 1
       end
   end
   return f(merged...)
end
function merge_transitions(σ1::SymbolicUtils.Term{<:Transition},σ2::SymbolicUtils.Term{<:Transition})
    i1,j1 = σ1.arguments[3:4]
    i2,j2 = σ2.arguments[3:4]
    if j1==i2
        return SymbolicUtils.term(σ1.f, σ1.arguments[1], σ1.arguments[2], i1, j2, acts_on(σ1), get_index(σ1); type=Transition)
    else
        return 0
    end
end
function rewrite_gs(t::SymbolicUtils.Term{<:Transition})
    h = t.arguments[1]
    aon = acts_on(t)
    index = get_index(t)
    aon_h = aon - index + 1
    gs = ground_state(h,aon_h)
    i,j = t.arguments[3:4]
    if i==j==gs
        args = Any[1]
        for k in levels(h,aon_h)
            (k==i) || push!(args, -1*SymbolicUtils.term(t.f, h, t.arguments[2], k, k, aon, index; type=Transition))
        end
        return +(args...)
    else
        return t
    end
end
