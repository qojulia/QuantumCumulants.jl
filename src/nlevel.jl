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
Base.isequal(h1::T,h2::T) where T<:NLevelSpace = isequal(h1.name, h2.name) && isequal(h1.levels, h2.levels) && isequal(h1.GS, h2.GS)
Base.hash(n::NLevelSpace, h::UInt) = hash(n.name, hash(n.levels, hash(n.GS, h)))

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
        @assert has_hilbert(NLevelSpace,hilbert,aon)
        @assert i∈levels(hilbert,aon) && j∈levels(hilbert,aon)
        op = new(hilbert,name,i,j,aon,index)
        if !haskey(OPERATORS_TO_SYMS, op)
            sym = SymbolicUtils.Sym{Transition}(gensym(:Transition))
            OPERATORS_TO_SYMS[op] = sym
            SYMS_TO_OPERATORS[sym] = op
        end
        return op
    end
end
Transition(hilbert::H,name::S,i::I,j::I,aon::A,index::IND=default_index()) where {H,S,I,A,IND<:Index} = Transition{H,S,I,A,IND}(hilbert,name,i,j,aon,index)
Transition(hilbert::NLevelSpace,name,i,j;index=default_index()) = Transition(hilbert,name,i,j,1,index)
function Transition(hilbert::ProductSpace,name,i,j;index=default_index())
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

Base.adjoint(t::Transition) = Transition(t.hilbert,t.name,t.j,t.i,acts_on(t),get_index(t))
Base.hash(t::Transition, h::UInt) = hash(t.hilbert, hash(t.name, hash(t.i, hash(t.j, hash(t.aon, hash(t.index, h))))))

# Simplification
istransition(x::Union{T,SymbolicUtils.Sym{T}}) where T<:Transition = true

function merge_transitions(f::Function, args)
    merged = Any[]
    i = 1
    while i <= length(args)
       if istransition(args[i])&&(i<length(args))&&istransition(args[i+1])&&(acts_on(args[i])==acts_on(args[i+1]))
           idx1 = get_index(args[i])
           idx2 = get_index(args[i+1])
           if isequal(idx1, idx2)
               push!(merged, merge_transitions(args[i],args[i+1]))
           else
               δ = _to_symbolic(idx1==idx2)
               ex1 = δ*merge_transitions(args[i],args[i+1])
               ex2 = SymbolicUtils.term(neq_inds_prod, [args[i], args[i+1]], [idx1!=idx2]; type=AbstractOperator)
               push!(merged, ex1+ex2)
           end

           i += 2
       else
           push!(merged, args[i])
           i += 1
       end
   end
   return f(merged...)
end
function merge_transitions(σ1::SymbolicUtils.Sym{<:Transition},σ2::SymbolicUtils.Sym{<:Transition})
    merge_transitions(_to_qumulants(σ1), _to_qumulants(σ2))
end
function merge_transitions(σ1::Transition, σ2::Transition)
    if acts_on(σ1)==acts_on(σ2) && isequal(σ1.index,σ2.index)
        i1,j1 = σ1.i, σ1.j
        i2,j2 = σ2.i, σ2.j
        if j1==i2
            op = Transition(σ1.hilbert,σ1.name,i1,j2,σ1.aon,σ1.index)
            return _to_symbolic(op)
        else
            return 0
        end
    end
    return nothing
end
rewrite_gs(t::SymbolicUtils.Sym{<:Transition}) = rewrite_gs(_to_qumulants(t))
function rewrite_gs(σ::Transition)
    h = σ.hilbert
    aon = acts_on(σ)
    gs = ground_state(h,aon)
    idx = get_index(σ)
    i,j = σ.i, σ.j
    if i==j==gs
        args = Any[1]
        for k in levels(h,aon)
            (k==i) || push!(args, -1*Transition(h, σ.name, k, k, aon, idx))
        end
        out = +(args...)
        return _to_symbolic(out)
    else
        return nothing
    end
end

function merge_transition_neq_prod(isthis, isthat, args)
    merged = Any[]
    i = 1
    while i <= length(args)
       if isthis(args[i])&&(i<length(args))&&isthat(args[i+1])&&(acts_on(args[i])==acts_on(args[i+1]))
           push!(merged, merge_transition_neq_prod(args[i], args[i+1]))
           i += 2
       else
           push!(merged, args[i])
           i += 1
       end
   end
   return *(merged...)
end
