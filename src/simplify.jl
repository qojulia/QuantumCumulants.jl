### Interface for SymbolicUtils

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
                   S::Type{<:AbstractOperator}) = AbstractOperator
end

SymbolicUtils.islike(::AbstractOperator, ::Type{<:Number}) = true
SymbolicUtils.symtype(x::T) where T<:AbstractOperator = T
SymbolicUtils.to_symbolic(x::AbstractOperator) = x

SymbolicUtils.istree(::BasicOperator) = false
SymbolicUtils.istree(::OperatorTerm) = true
SymbolicUtils.arguments(t::OperatorTerm) = t.arguments
for f in [*,+,-,/,^]
    @eval SymbolicUtils.operation(t::OperatorTerm{<:$(typeof(f))}) = $f
end
SymbolicUtils.similarterm(t::OperatorTerm{F}, f::F, args::A) where {F,A} = OperatorTerm{F,A}(t.f, args)

SymbolicUtils.:<ₑ(a::AbstractOperator, b::Number) = false
SymbolicUtils.:<ₑ(a::Number,   b::AbstractOperator) = true
SymbolicUtils.:<ₑ(a::AbstractOperator, b::SymbolicUtils.Symbolic{<:Number}) = false
SymbolicUtils.:<ₑ(a::SymbolicUtils.Symbolic{<:Number}, b::AbstractOperator) = true


### End of interface

"""
    simplify_operators(op::AbstractOperator; rewriter=default_operator_simplifier(), kwargs...)

Simplify an operator through standard algebraic rewriting, as well as using
fundamental commutation relations.

# Arguments
===========
* op: The operator expression to be simplified.
* rewriter: The rewriter used.
* kwargs: Further arguments passed to `SymbolicUtils.simplify`.
"""
function simplify_operators(op; rewriter=default_operator_simplifier(),
                kwargs...)
    return SymbolicUtils.simplify(op; rewriter=rewriter, kwargs...)
end

"""
    expand(ex; rewriter=defualt_expand_simplifier(), kwargs...)

Simple wrapper around `SymbolicUtils.simplify` that uses a rewriter such
that expressions are expanded.

# Arguments
===========
* ex: The expression to be expanded.
* rewriter: The used rewriter.
* kwargs: Further arguments passed to `SymbolicUtils.simplify`.

Examples
========
```
julia> @parameters p q r
(p, q, r)

julia> ex = p*(q+r) + (q+p)*(r+q)
((p*(q+r))+((q+p)*(r+q)))

julia> expand(ex)
((p*q)+(p*r)+(q*r)+(p*r)+(q*q)+(p*q))
```
"""
function expand(ex; rewriter=default_expand_simplifier(), kwargs...)
    return SymbolicUtils.simplify(ex; rewriter=rewriter, kwargs...)
end

### Functions needed for simplification

# Handle noncommutative multiplication
iscommutative(::AbstractOperator) = false
iscommutative(::Union{SymbolicUtils.Symbolic{T},T}) where {T<:Number} = true

needs_sorting_nc(x) = (x.f === (*)) && !issorted_nc(x)
needs_sorting_nc(x::SymbolicUtils.Mul{<:Number}) = SymbolicUtils.needs_sorting(*)(x)
function issorted_nc(x)
    args = SymbolicUtils.arguments(x)
    is_c = map(iscommutative, args)
    issorted(is_c, lt=(>)) || return false
    args_c = args[is_c]
    SymbolicUtils.issortedₑ(args_c) || return false
    args_nc = args[.!is_c]
    return issorted(args_nc, lt=lt_aon)
end

# Comparison for sorting according to Hilbert spaces
function lt_aon(t1,t2)
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

using SymbolicUtils: <ₑ
function sort_args_nc(x)
    args = SymbolicUtils.arguments(x)
    is_c = iscommutative.(args)
    args_c = sort(args[is_c], lt=(<ₑ))
    args_nc = sort(args[.!is_c], lt=lt_aon)
    return SymbolicUtils.similarterm(x, *, vcat(args_c, args_nc))
end

# Apply commutation relation
function apply_commutator(fcomm, args_l, args_r, a, b)
    if acts_on(a)==acts_on(b)
        global tmp = (args_l, fcomm(a,b), args_r)
        return *(args_l..., fcomm(a, b), args_r...)
    else
        return nothing
    end
end
