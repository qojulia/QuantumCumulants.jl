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
function qsimplify(x;
                  polynorm=false,
                  threaded=false,
                  thread_subtree_cutoff=100,
                  rewriter=nothing)
    f = if rewriter === nothing
        if threaded
            if SymbolicUtils.symtype(x) <: QNumber
                threaded_q_simplifier(thread_subtree_cutoff)
            else
                threaded_c_simplifier(thread_subtree_cutoff)
            end
        elseif polynorm
            if SymbolicUtils.symtype(x) <: QNumber
                error("Cannot use `polynormalize` on `QNumber`!")
            else
                serial_c_polynorm
            end
        else
            if SymbolicUtils.symtype(x) <: QNumber
                serial_q_simplifier
            else
                serial_c_simplifier
            end
        end
    else
        SymbolicUtils.Fixpoint(rewriter)
    end

    SymbolicUtils.PassThrough(f)(x)
end

"""
    expand(ex; rewriter=default_expand_simplifier(), kwargs...)

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
function expand(x;
                polynorm=false,
                threaded=false,
                thread_subtree_cutoff=100,
                rewriter=nothing)
    f = if rewriter === nothing
        if threaded
            threaded_expand_simplifier(thread_subtree_cutoff)
        elseif polynorm
            if SymbolicUtils.symtype(x) <: QNumber
                error("Cannot use `polynormalize` on `QNumber`!")
            else
                serial_expand_polynorm
            end
        else
            serial_expand_simplifier
        end
    else
        SymbolicUtils.Fixpoint(rewriter)
    end

    SymbolicUtils.PassThrough(f)(x)
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
