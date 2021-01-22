### Conversion to SymbolicUtils

_to_symbolic(t::T) where T<:OperatorTerm = SymbolicUtils.Term{AbstractOperator}(t.f, _to_symbolic.(t.arguments))
_to_symbolic(x::Number) = x
_to_symbolic(x::SymbolicUtils.Symbolic) = x
_to_symbolic(op::BasicOperator) = OPERATORS_TO_SYMS[op]
_to_symbolic(x::Nothing) = x # can be return from simplification

_to_qumulants(t::SymbolicUtils.Sym{T}) where T<:BasicOperator = SYMS_TO_OPERATORS[t]
function _to_qumulants(t::SymbolicUtils.Term{T}) where T<:AbstractOperator
    return OperatorTerm(t.f, _to_qumulants.(t.arguments))
end
_to_qumulants(x::Number) = x

for f in [:acts_on, :hilbert, :levels]
    @eval $f(s::SymbolicUtils.Sym{<:BasicOperator}, args...) = $f(_to_qumulants(s), args...)
end

# Symbolic type promotion
SymbolicUtils.promote_symtype(f, Ts::Type{<:AbstractOperator}...) = promote_type(AbstractOperator,Ts...)
SymbolicUtils.promote_symtype(f, T::Type{<:AbstractOperator}, Ts...) = promote_type(AbstractOperator,T)
SymbolicUtils.promote_symtype(f, T, S, Ts::Union{Type{<:Number},Type{<:AbstractOperator}}...) = SymbolicUtils.promote_symtype(f, SymbolicUtils.promote_symtype(f, T, S), Ts...)
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

# SymbolicUtils.assert_number(::SymbolicUtils.Symbolic{<:AbstractOperator}) = true
SymbolicUtils.islike(::SymbolicUtils.Symbolic{<:AbstractOperator}, ::Type{<:Number}) = true

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
function simplify_operators(op::AbstractOperator; rewriter=default_operator_simplifier(),
                kwargs...)
    s = _to_symbolic(op)
    s_ = SymbolicUtils.simplify(s; rewriter=rewriter, kwargs...)
    (SymbolicUtils.symtype(s_) == Any) && @warn "SymbolicUtils.simplify returned symtype Any; recursion failed!"
    return _to_qumulants(s_)
end
simplify_operators(x::Number, args...; kwargs...) = x

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
    s = _to_symbolic(ex)
    s_ = SymbolicUtils.simplify(s; rewriter=rewriter, kwargs...)
    return _to_qumulants(s_)
end

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
iscommutative(::AbstractOperator) = false
iscommutative(::SymbolicUtils.Symbolic{<:AbstractOperator}) = false
iscommutative(::Union{SymbolicUtils.Symbolic{T},T}) where {T<:Number} = true

needs_sorting_nc(x) = (x.f === (*)) && !issorted_nc(x)
function issorted_nc(x)
    args = SymbolicUtils.arguments(x)
    is_c = iscommutative.(args)
    args_c = args[is_c]
    args_nc = args[.!is_c]
    return issorted(is_c, lt=(>)) && SymbolicUtils.issortedₑ(args_c) && issorted(args_nc, lt=lt_aon)
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
