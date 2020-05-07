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
            sym = SymbolicUtils.Sym{Destroy}(gensym(:Destroy))
            OPERATORS_TO_SYMS[op] = sym
            SYMS_TO_OPERATORS[sym] = op
        end
        return op
    end
end
Destroy(hilbert::H,name::S) where {H,S} = Destroy{H,S}(hilbert,name)
isdestroy(a) = false
isdestroy(a::Union{SymbolicUtils.Sym{T},T}) where T<:Destroy = true


struct Create{H<:FockSpace,S} <: BasicOperator
    hilbert::H
    name::S
    function Create{H,S}(hilbert::H,name::S) where {H,S}
        op = new(hilbert,name)
        if !haskey(OPERATORS_TO_SYMS,op)
            sym = SymbolicUtils.Sym{Create}(gensym(:Create))
            OPERATORS_TO_SYMS[op] = sym
            SYMS_TO_OPERATORS[sym] = op
        end
        return op
    end
end
Create(hilbert::H,name::S) where {H,S} = Create{H,S}(hilbert,name)
iscreate(a) = false
iscreate(a::Union{SymbolicUtils.Sym{T},T}) where T<:Create = true

Base.adjoint(op::Destroy) = Create(op.hilbert,op.name)
Base.adjoint(op::Create) = Destroy(op.hilbert,op.name)
Base.isone(::BasicOperator) = false
Base.iszero(::BasicOperator) = false

# Commutation relation in simplification
function has_destroy_create(args)
    length(args) <= 1 && return false
    for i=1:length(args)-1
        if isdestroy(args[i])&&iscreate(args[i+1])
            return true
        end
    end
    return false
end
function commute_bosonic(f,args)
    commuted_args = []
    i = 1
    last_commute = false
    while i < length(args)
        if isdestroy(args[i]) && iscreate(args[i+1])
            push!(commuted_args, args[i+1]*args[i] + one(args[i]))
            i += 2
            last_commute = (i==length(args))
        else
            push!(commuted_args, args[i])
            i += 1
        end
    end
    last_commute && push!(commuted_args, args[end])
    return SymbolicUtils.Term(f, commuted_args)
end
