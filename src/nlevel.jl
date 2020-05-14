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
NLevelSpace(name,levels) = NLevelSpace(name,levels,levels[1])
Base.:(==)(h1::T,h2::T) where T<:NLevelSpace = (h1.name==h2.name && h1.levels==h2.levels && h1.GS==h2.GS)

struct Transition{H<:NLevelSpace,S,I,J} <: BasicOperator
    hilbert::H
    name::S
    i::I
    j::J
    function Transition{H,S,I,J}(hilbert::H,name::S,i::I,j::J) where {H,S,I,J}
        @assert i∈hilbert.levels && j∈hilbert.levels
        op = new(hilbert,name,i,j)
        if !haskey(OPERATORS_TO_SYMS,op)
            sym = SymbolicUtils.Sym{Transition}(gensym(:Transition))
            OPERATORS_TO_SYMS[op] = sym
            SYMS_TO_OPERATORS[sym] = op
        end
        return op
    end
end
Transition(hilbert::H,name::S,i::I,j::J) where {H,S,I,J} = Transition{H,S,I,J}(hilbert,name,i,j)
Base.adjoint(t::Transition) = Transition(t.hilbert,t.name,t.j,t.i)

istransition(x) = false
istransition(x::Union{T,SymbolicUtils.Sym{T}}) where T<:Transition = true

function has_transitions(args)
    for i=1:length(args)-1
        if istransition(args[i])&&istransition(args[i+1])
            return true
        end
    end
    return false
end
function merge_transitions(f, args)
    merged = Any[]

    i = 1
    last_merge = false
    while i<length(args)
       if istransition(args[i])&&istransition(args[i+1])
           push!(merged, merge_transitions(args[i],args[i+1]))
           i += 2
           last_merge = (i==length(args))
       else
           push!(merged, args[i])
           i += 1
       end
   end
   last_merge && push!(merged, args[i])
   return SymbolicUtils.Term(f,merged)
end
function merge_transitions(σ1::Transition,σ2::Transition)
    check_hilbert(σ1,σ2)
    if σ1.j==σ2.i
        return Transition(σ1.hilbert,σ1.name,σ1.i,σ2.j)
    else
        return zero(σ1)
    end
end
function merge_transitions(σ1::SymbolicUtils.Sym{<:Transition},σ2::SymbolicUtils.Sym{<:Transition})
    op1 = _to_qumulants(σ1)
    op2 = _to_qumulants(σ2)
    return _to_symbolic(merge_transitions(op1,op2))
end
function rewrite_gs(x)
    op = _to_qumulants(x)
    if op.i==op.j==op.hilbert.GS
        op_ = one(op)
        for i in op.hilbert.levels
            if i!=op.hilbert.GS
                op_ += (-1*Transition(op.hilbert,op.name,i,i))
            end
        end
        return _to_symbolic(op_)
    else
        return x
    end
end
