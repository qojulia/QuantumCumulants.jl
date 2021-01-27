### Interface for SymbolicUtils

# Symbolic type promotion
SymbolicUtils.promote_symtype(f, Ts::Type{<:QNumber}...) = promote_type(QNumber,Ts...)
SymbolicUtils.promote_symtype(f, T::Type{<:QNumber}, Ts...) = promote_type(QNumber,T)
for f in [+,-,*,/,^]
    @eval SymbolicUtils.promote_symtype(::$(typeof(f)),
                   T::Type{<:QNumber},
                   S::Type{<:Number}) = QNumber
    @eval SymbolicUtils.promote_symtype(::$(typeof(f)),
                   T::Type{<:Number},
                   S::Type{<:QNumber}) = QNumber
    @eval SymbolicUtils.promote_symtype(::$(typeof(f)),
                   T::Type{<:QNumber},
                   S::Type{<:QNumber}) = QNumber
end

SymbolicUtils.islike(::QNumber, ::Type{<:Number}) = true
SymbolicUtils.symtype(x::T) where T<:QNumber = T
SymbolicUtils.to_symbolic(x::QNumber) = x

SymbolicUtils.istree(::QSym) = false
SymbolicUtils.istree(::QTerm) = true
SymbolicUtils.arguments(t::QTerm) = t.arguments
for f in [*,+,-,/,^]
    @eval SymbolicUtils.operation(t::QTerm{<:$(typeof(f))}) = $f
end
SymbolicUtils.similarterm(t::QTerm{F}, f::F, args::A) where {F,A} = QTerm{F,A}(t.f, args)

SymbolicUtils.:<ₑ(a::QNumber, b::Number) = false
SymbolicUtils.:<ₑ(a::Number,   b::QNumber) = true
SymbolicUtils.:<ₑ(a::QNumber, b::SymbolicUtils.Symbolic{<:Number}) = false
SymbolicUtils.:<ₑ(a::SymbolicUtils.Symbolic{<:Number}, b::QNumber) = true


### End of interface

"""
    qsimplify(op::QNumber; rewriter=default_operator_simplifier(), kwargs...)

Simplify an operator through standard algebraic rewriting, as well as using
fundamental commutation relations.

# Arguments
===========
* op: The operator expression to be simplified.
* rewriter: The rewriter used.
* kwargs: Further arguments passed to `SymbolicUtils.simplify`.
"""
function qsimplify(op; rewriter=default_operator_simplifier(),
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
julia> @params p q r
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
iscommutative(::QNumber) = false
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
        return *(args_l..., fcomm(a, b), args_r...)
    else
        return nothing
    end
end
