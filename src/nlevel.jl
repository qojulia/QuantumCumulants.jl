struct NLevelSpace{S,L,G} <: HilbertSpace
    name::S
    levels::L
    GS::G
    function NLevelSpace{S,L,G}(name::S,levels::L,GS::G) where {S,L,G}
        r1 = SymbolicUtils.@rule(*(~~x::has_transitions) => merge_transitions(*, ~~x))
        r2 = SymbolicUtils.@rule(~x::istransition => rewrite_gs(~x))
        (r1 ∈ COMMUTATOR_RULES.rules) || push!(COMMUTATOR_RULES.rules, r1)
        (r2 ∈ COMMUTATOR_RULES.rules) || push!(COMMUTATOR_RULES.rules, r2)
        new(name,levels,GS)
    end
end
NLevelSpace(name::S,levels::L,GS::G) where {S,L,G} = NLevelSpace{S,L,G}(name,levels,GS)
NLevelSpace(name,N::Int,GS) = NLevelSpace(name,1:N,GS)
NLevelSpace(name,N::Int) = NLevelSpace(name,1:N,1)
NLevelSpace(name,levels) = NLevelSpace(name,levels,levels[1])
Base.:(==)(h1::T,h2::T) where T<:NLevelSpace = (h1.name==h2.name && h1.levels==h2.levels && h1.GS==h2.GS)

levels(h::NLevelSpace) = h.levels
levels(h::NLevelSpace,aon) = levels(h)
levels(h::ProductSpace,aon) = levels(h.spaces[aon])
ground_state(h::NLevelSpace) = h.GS
ground_state(h::NLevelSpace,aon) = h.GS
ground_state(h::ProductSpace,aon) = ground_state(h.spaces[aon])

struct Transition{H,S,I,A} <: BasicOperator
    hilbert::H
    name::S
    i::I
    j::I
    aon::A
    function Transition{H,S,I,A}(hilbert::H,name::S,i::I,j::I,aon::A) where {H,S,I,A}
        @assert i∈levels(hilbert,aon) && j∈levels(hilbert,aon)
        op = new(hilbert,name,i,j,aon)
        if !haskey(OPERATORS_TO_SYMS,op)
            sym = generate_symbolic(op)
            OPERATORS_TO_SYMS[op] = sym
            SYMS_TO_OPERATORS[sym] = op
        end
        return op
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

function embed(h::ProductSpace,op::T,aon::Int) where T<:Transition
    check_hilbert(h.spaces[aon],op.hilbert)
    op_ = Transition(h,op.name,op.i,op.j,aon)
    return op_
end
levels(t::Transition,args...) = levels(t.hilbert,args...)
ground_state(t::Transition,args...) = ground_state(t.hilbert,args...)

Base.adjoint(t::Transition) = Transition(t.hilbert,t.name,t.j,t.i,acts_on(t))
Base.:(==)(t1::Transition,t2::Transition) = (t1.hilbert==t2.hilbert && t1.name==t2.name && t1.i==t2.i && t1.j==t2.j)

generate_symbolic(op::Transition) = SymbolicUtils.Sym{SymbolicUtils.FnType{Tuple{Int},Transition}}(gensym(:Transition))

istransition(x) = false
istransition(x::Union{T,SymbolicUtils.Term{T}}) where T<:Transition = true

function has_transitions(args)
    for i=1:length(args)-1
        if istransition(args[i])&&istransition(args[i+1])&&(acts_on(args[i])==acts_on(args[i+1]))
            return true
        end
    end
    return false
end
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
function merge_transitions(σ1::Transition,σ2::Transition)
    check_hilbert(σ1,σ2)
    if σ1.j==σ2.i
        return Transition(σ1.hilbert,σ1.name,σ1.i,σ2.j,acts_on(σ1))
    else
        return zero(σ1)
    end
end
function merge_transitions(σ1::SymbolicUtils.Term{<:Transition},σ2::SymbolicUtils.Term{<:Transition})
    op1 = _to_qumulants(σ1)
    op2 = _to_qumulants(σ2)
    return _to_symbolic(merge_transitions(op1,op2))
end
function rewrite_gs(x)
    op = _to_qumulants(x)
    return _to_symbolic(rewrite_gs(op))
end
function rewrite_gs(op::Transition)
    if op.i==op.j==ground_state(op.hilbert,acts_on(op))
        args = Any[one(op)]
        for i in levels(op.hilbert,acts_on(op))
            (i==op.i) || push!(args, -1*Transition(op.hilbert,op.name,i,i,acts_on(op)))
        end
        return +(args...)
    else
        return op
    end
end
