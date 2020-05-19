struct FockSpace{S} <: HilbertSpace
    name::S
    function FockSpace(name::S) where S
        r = SymbolicUtils.@rule(*(~~x::has_destroy_create) => commute_bosonic(*, ~~x))
        (r ∈ COMMUTATOR_RULES.rules) || push!(COMMUTATOR_RULES.rules, r)
        new{S}(name)
    end
end
Base.:(==)(h1::T,h2::T) where T<:FockSpace = (h1.name==h2.name)


struct Destroy{H<:FockSpace,S} <: BasicOperator
    hilbert::H
    name::S
    function Destroy{H,S}(hilbert::H,name::S) where {H,S}
        op = new(hilbert,name)
        if !haskey(OPERATORS_TO_SYMS,op)
            sym = generate_symbolic(op)
            OPERATORS_TO_SYMS[op] = sym
            SYMS_TO_OPERATORS[sym] = op
        end
        return op
    end
end
Destroy(hilbert::H,name::S) where {H,S} = Destroy{H,S}(hilbert,name)
isdestroy(a) = false
isdestroy(a::SymbolicUtils.Term{T}) where {T<:Destroy} = true

struct Create{H<:FockSpace,S} <: BasicOperator
    hilbert::H
    name::S
    function Create{H,S}(hilbert::H,name::S) where {H,S}
        op = new(hilbert,name)
        if !haskey(OPERATORS_TO_SYMS,op)
            sym = SymbolicUtils.Sym{SymbolicUtils.FnType{Tuple{Int},Create}}(gensym(:Create))
            OPERATORS_TO_SYMS[op] = sym
            SYMS_TO_OPERATORS[sym] = op
        end
        return op
    end
end
Create(hilbert::H,name::S) where {H,S} = Create{H,S}(hilbert,name)
iscreate(a) = false
iscreate(a::SymbolicUtils.Term{T}) where {T<:Create} = true

Base.adjoint(op::Destroy) = Create(op.hilbert,op.name)
Base.adjoint(op::Create) = Destroy(op.hilbert,op.name)

generate_symbolic(::Destroy) =
    SymbolicUtils.Sym{SymbolicUtils.FnType{Tuple{Int},Destroy}}(gensym(:Destroy))
generate_symbolic(::Create) =
    SymbolicUtils.Sym{SymbolicUtils.FnType{Tuple{Int},Create}}(gensym(:Create))

# Commutation relation in simplification
function has_destroy_create(args)
    length(args) <= 1 && return false
    for i=1:length(args)-1
        if isdestroy(args[i])&&iscreate(args[i+1])&&(acts_on(args[i])==acts_on(args[i+1]))
            return true
        end
    end
    return false
end
function commute_bosonic(f,args)
    commuted_args = []
    i = 1
    while i <= length(args)
        if isdestroy(args[i]) && i<length(args) && iscreate(args[i+1]) && (acts_on(args[i])==acts_on(args[i+1]))
            push!(commuted_args, args[i+1]*args[i] + 1)
            i += 2
        else
            push!(commuted_args, args[i])
            i += 1
        end
    end
    return f(commuted_args...)
end
