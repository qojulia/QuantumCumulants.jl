function _is_avg_leaf(x::SymbolicUtils.BasicSymbolic)
    return SQA.is_average(x) &&
        SymbolicUtils.iscall(x) &&
        SymbolicUtils.operation(x) === SQA.sym_average
end
_is_avg_leaf(::Any) = false

"""
Collect the leaf averages of `x`, left to right, into a `Vector`. A lazy iterator
is unnecessary at the bounded build-time sizes seen here.
"""
eachleaf(x::Symbolics.Num) = eachleaf(SymbolicUtils.unwrap(x))
function eachleaf(x)
    out = SymbolicUtils.BasicSymbolic[]
    _eachleaf!(out, x)
    return out
end
_eachleaf!(::Any, ::Any) = nothing
function _eachleaf!(out, x::SymbolicUtils.BasicSymbolic)
    if _is_avg_leaf(x)
        push!(out, x)
        return nothing
    end
    SymbolicUtils.iscall(x) || return nothing
    for a in SymbolicUtils.arguments(x)
        _eachleaf!(out, a)
    end
    return nothing
end

"""
Reduce `f` over the leaf averages of `x`, threading `acc`, without rebuilding the
expression tree.
"""
foldleaves(f, acc, x::Symbolics.Num) = foldleaves(f, acc, SymbolicUtils.unwrap(x))
function foldleaves(f, acc, x)
    x isa SymbolicUtils.BasicSymbolic || return acc
    _is_avg_leaf(x) && return f(acc, x)
    SymbolicUtils.iscall(x) || return acc
    for a in SymbolicUtils.arguments(x)
        acc = foldleaves(f, acc, a)
    end
    return acc
end

"""
Rewrite each leaf average of `x` via `f`, rebuilding only the branches that
actually changed.
"""
mapleaves(f, x::Symbolics.Num) = Symbolics.wrap(mapleaves(f, SymbolicUtils.unwrap(x)))
function mapleaves(f, x)
    x isa SymbolicUtils.BasicSymbolic || return x
    _is_avg_leaf(x) && return f(x)
    SymbolicUtils.iscall(x) || return x
    op = SymbolicUtils.operation(x)
    args = SymbolicUtils.arguments(x)
    new_args = Any[mapleaves(f, a) for a in args]
    all(((a, b),) -> a === b, zip(args, new_args)) && return x
    return _rebuild(op, new_args, x)
end

"""
Rebuild a node from its mapped arguments. A `complex(re, im)` call is rewritten to
the `re + im*IM` form, because SymbolicUtils does not unify a `complex(0, …)`
literal with the factored symbolic form. Falls back to `maketerm` when the
operation is not directly callable on the new arguments.
"""
function _rebuild(op, args, x)
    if op === complex && length(args) == 2
        return args[1] + args[2] * Symbolics.IM
    end
    try
        return op(args...)
    catch
        return TermInterface.maketerm(typeof(x), op, args, TermInterface.metadata(x))
    end
end
