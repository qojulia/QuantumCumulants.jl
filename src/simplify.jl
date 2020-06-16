### Conversion to SymbolicUtils

_to_symbolic(t::T) where T<:OperatorTerm = SymbolicUtils.Term{AbstractOperator}(t.f, _to_symbolic.(t.arguments))
_to_symbolic(x::Number) = x
_to_symbolic(x::SymbolicUtils.Symbolic) = x
_to_symbolic(op::BasicOperator) = OPERATORS_TO_SYMS[op]

_to_qumulants(t::SymbolicUtils.Sym{T}) where T<:BasicOperator = SYMS_TO_OPERATORS[t]
function _to_qumulants(t::SymbolicUtils.Term{T}) where T<:AbstractOperator
    return OperatorTerm(t.f, _to_qumulants.(t.arguments))
end
_to_qumulants(x::Number) = x

for f in [:acts_on, :hilbert, :levels, :get_index]
    @eval $f(s::SymbolicUtils.Sym{<:BasicOperator}, args...) = $f(_to_qumulants(s), args...)
end

# Symbolic type promotion
SymbolicUtils.promote_symtype(f, Ts::Type{<:AbstractOperator}...) = promote_type(AbstractOperator,Ts...)
SymbolicUtils.promote_symtype(f, T::Type{<:AbstractOperator}, Ts...) = promote_type(AbstractOperator,T)
for f in [+,-,*,/,^]
    @eval SymbolicUtils.promote_symtype(::$(typeof(f)),
                   T::Type{<:AbstractOperator},
                   S::Type{<:Number}) = AbstractOperator
    @eval SymbolicUtils.promote_symtype(::$(typeof(f)),
                   T::Type{<:Number},
                   S::Type{<:AbstractOperator}) = AbstractOperator
    @eval SymbolicUtils.promote_symtype(::$(typeof(f)),
                   T::Type{<:AbstractOperator},
                   S::Type{<:AbstractOperator}) = AbstractOperator#promote_type(T,S)
end

Base.one(x::SymbolicUtils.Symbolic{T}) where T<:AbstractOperator = 1
Base.zero(x::SymbolicUtils.Symbolic{T}) where T<:AbstractOperator = 0
Base.one(x::SymbolicUtils.Sym{SymbolicUtils.FnType{A,T}}) where {A,T<:AbstractOperator} = 1
Base.zero(x::SymbolicUtils.Sym{SymbolicUtils.FnType{A,T}}) where {A,T<:AbstractOperator} = 0

### End of interface

"""
    simplify_operators(op::AbstractOperator; rules=default_rules(), kwargs...)

Simplify an operator through standard algebraic rewriting, as well as using
fundamental commutation relations.

# Arguments
===========
* op: The operator expression to be simplified
* rules: The rules used. Defaults to SIMPLIFY_OPERATOR_RULES
* kwargs: Further arguments passed to `SymbolicUtils.simplify`.
"""
function simplify_operators(op::AbstractOperator; rules=default_rules(),
                commutator_rules=default_commutator_rules(),
                kwargs...)
    s = _to_symbolic(op)
    s_ = _simplify_operators(s, rules, commutator_rules; kwargs...)
    (SymbolicUtils.symtype(s_) == Any) && @warn "SymbolicUtils.simplify returned symtype Any; recursion failed!"
    return _to_qumulants(s_)
end
function _simplify_operators(s, rules, commutator_rules; kwargs...)
    s_ = SymbolicUtils.simplify(s; rules=commutator_rules)
    s_ = SymbolicUtils.simplify(s_; rules=rules, kwargs...)
    return s_
end
simplify_operators(x::Number, args...; kwargs...) = x

function expand(ex; rules=default_expand_rules(), fixpoint=true, applyall=true, recurse=true, kwargs...)
    s = _to_symbolic(ex)
    s_ = SymbolicUtils.simplify(s; rules=rules, fixpoint=fixpoint, applyall=applyall, recurse=recurse, kwargs...)
    return _to_qumulants(s_)
end

default_rules() = SIMPLIFY_OPERATOR_RULES
default_commutator_rules() = SIMPLIFY_COMMUTATOR_RULES
default_expand_rules() = EXPAND_RULES

"""
    substitute(arg, subs; simplify=true)

Substitute the symbolic argument, i.e. any subtype to [`AbstractOperator`](@ref)
or [`SymbolicNumber`](@ref) according to the substitutions stored in a `Dict`.
Also works on [`DifferentialEquation`](@ref). If `simplify=true`, the output
is simplified.

Examples
=======
```
julia> @parameters p
(p,)

julia> substitute(p, Dict(p=>2))
2
```
"""
function substitute(op::BasicOperator, dict; kwargs...)
    if haskey(dict, op)
        op_ = dict[op]
        check_hilbert(op_,op)
        return op_
    elseif haskey(dict, op')
        op_ = dict[op']
        check_hilbert(op_, op)
        return op_'
    else
        return op
    end
end
function substitute(t::OperatorTerm, dict; simplify=true, kwargs...)
    if haskey(dict, t)
        return dict[t]
    elseif haskey(dict, t')
        return dict[t']'
    else
        if simplify
            return simplify_operators(t.f([substitute(arg, dict; simplify=simplify) for arg in t.arguments]...), kwargs...)
        else
            return t.f([substitute(arg, dict; simplify=simplify) for arg in t.arguments]...)
        end
    end
end
substitute(x::Number, dict; kwargs...) = x

### Functions needed for simplification

# Handle noncommutative multiplication
iscommutative(::typeof(*), x::Union{SymbolicUtils.Symbolic{SymbolicUtils.FnType{A,T}},T}) where {A,T<:AbstractOperator} = false
iscommutative(::typeof(*), x::SymbolicUtils.Symbolic{T}) where {T<:AbstractOperator} = false
iscommutative(::typeof(*), x::Union{SymbolicUtils.Symbolic{T},T}) where {T<:Number} = true
iscommutative(f) = x -> iscommutative(f, x)

issorted_nc(f) = x -> issorted_nc(f, x)
function issorted_nc(f::typeof(*), args)
    is_c = iscommutative.(f, args)
    args_c = args[is_c]
    args_nc = args[.!is_c]
    return issorted(is_c, lt=(>)) && SymbolicUtils.issortedₑ(args_c) && issorted(args_nc, lt=lt_aon)
end

# Comparison for sorting according to Hilbert spaces
function lt_aon(t1::SymbolicUtils.Symbolic{<:AbstractOperator},t2::SymbolicUtils.Symbolic{<:AbstractOperator})
    aon1 = acts_on(t1)
    aon2 = acts_on(t2)
    if any(a1 ∈ aon2 for a1 in aon1)
        return false
    elseif any(a2 ∈ aon1 for a2 in aon2)
        return false
    elseif isempty(aon1)
        return isempty(aon2)
    elseif isempty(aon2)
        return isempty(aon1)
    else
        return maximum(aon1)<maximum(aon2)
    end
end

function acts_on(t::SymbolicUtils.Term{T}) where T<:AbstractOperator
    ops = filter(isoperator, t.arguments)
    aon = Int[]
    for op in ops
        append!(aon, acts_on(op))
    end
    unique!(aon)
    sort!(aon)
    return aon
end

function sort_args_nc(f::typeof(*), args)
    is_c = iscommutative.(f, args)
    args_c = SymbolicUtils.sort_args(f, args[is_c]).arguments
    args_nc = sort(args[.!is_c], lt=lt_aon)
    return f(args_c..., args_nc...)
end

# Expand products involving sums into sums of products, e.g. x*(y+z) => x*y + x*z; needed to apply commutator rules
has_inner(f) = x -> has_inner(f,x)
function has_inner(f,args)
    for t in args
        if t isa SymbolicUtils.Term && SymbolicUtils.operation(t) === (f)
            return true
        end
    end
    return false
end

function expand_term(f_outer, f_inner, args)
    expanded_args = []
    for i=1:length(args)
        if args[i] isa SymbolicUtils.Term && SymbolicUtils.operation(args[i]) === f_inner
            inner_args = SymbolicUtils.arguments(args[i])

            # Arguments left and right of the nested f_inner
            args_l = args[1:i-1]
            args_r = args[i+1:end]
            for inner_arg in inner_args
                push!(expanded_args, f_outer(args_l..., inner_arg, args_r...))
            end
            return f_inner(expanded_args...)
        end
    end
end

# Check if specific consecutive arguments occur
has_consecutive(isthis) = has_consecutive(isthis,isthis)
has_consecutive(isthis,isthat) = x -> has_consecutive(isthis,isthat,x)
function has_consecutive(isthis,isthat,args)
    length(args) <= 1 && return false
    for i=1:length(args)-1
        if isthis(args[i])&&isthat(args[i+1])&&(acts_on(args[i])==acts_on(args[i+1]))
            return true
        end
    end
    return false
end
