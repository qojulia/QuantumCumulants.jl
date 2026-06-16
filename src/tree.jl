function _is_avg_leaf(x::SymbolicUtils.BasicSymbolic)
    SQA.is_average(x) || return false
    (SymbolicUtils.iscall(x) && SymbolicUtils.operation(x) === SQA.sym_average) && return true
    return SymbolicUtils.hasmetadata(x, SQA.AverageOperator)
end
_is_avg_leaf(::Any) = false

"""
The shared descent condition: recurse into a node only when it is a call and not an
average leaf. The `iscall` guard keeps read traversals from recursing into numbers and
bare symbols.
"""
_descend(x) = SymbolicUtils.iscall(x) && !_is_avg_leaf(x)

"""
Reconstruct a node from `op` and (possibly rewritten) `args`. A `complex(re, im)` call is
rewritten to `re + im*IM`, because SymbolicUtils does not unify a `complex(0, …)` literal
with the factored symbolic form. Otherwise the operation is applied directly (so SQA
operator algebra is re-run), falling back to `maketerm` with the original `meta` when the
operation is not callable on the new arguments. This is the single reconstruction path for
the whole package.
"""
function _qc_maketerm(T, op, args, meta)
    if op === complex && length(args) == 2
        return args[1] + args[2] * Symbolics.IM
    end
    try
        return op(args...)
    catch err
        # Fall back only when `op` is not callable on these args; other errors are real.
        (err isa MethodError || err isa ArgumentError) || rethrow()
        return TermInterface.maketerm(T, op, args, meta)
    end
end

"""
    rewrite(rule, x; descend=_descend, post=nothing, maketerm=_qc_maketerm)

Post-order rewrite of the expression tree `x`. At each node `rule(node)` is consulted
first: a non-`nothing` return replaces the subtree and stops descent into it (this is how
leaf transforms and identity substitutions are expressed). Otherwise, when `descend(node)`
holds, the arguments are rewritten and the node is rebuilt via `maketerm` only if an
argument changed (cheap `===` identity check). An optional `post(node)` is applied to every
descended node after reconstruction, for per-node metadata edits. `maketerm` defaults to
[`_qc_maketerm`](@ref) (re-applies operator algebra, normalizes `complex`); the scope
walkers pass [`_structural_maketerm`](@ref) to rebuild without re-simplification.
"""
rewrite(rule, x::Symbolics.Num; kw...) = Symbolics.wrap(rewrite(rule, SymbolicUtils.unwrap(x); kw...))
rewrite(rule, x; descend = _descend, post = nothing, maketerm = _qc_maketerm) =
    _rewrite(rule, x, descend, post, maketerm)
function _rewrite(rule, x, descend::D, post::P, maketerm::M) where {D, P, M}
    x isa SymbolicUtils.BasicSymbolic || return x
    r = rule(x)
    r === nothing || return r
    descend(x) || return x
    op = SymbolicUtils.operation(x)
    args = SymbolicUtils.arguments(x)
    # `===` not `isequal`: `isequal` ignores metadata and would collapse distinct
    # averages that differ only in their scope metadata.
    new_args = nothing
    for i in eachindex(args)
        a = args[i]
        b = _rewrite(rule, a, descend, post, maketerm)
        if new_args !== nothing
            @inbounds new_args[i] = b
        elseif b !== a
            new_args = Vector{Any}(undef, length(args))
            @inbounds copyto!(new_args, 1, args, 1, i - 1)
            @inbounds new_args[i] = b
        end
    end
    y = new_args === nothing ? x : maketerm(typeof(x), op, new_args, TermInterface.metadata(x))
    return post === nothing ? y : post(y)
end

"""Structural reconstruction (metadata-preserving `maketerm`, no operator re-application)."""
_structural_maketerm(T, op, args, meta) = TermInterface.maketerm(T, op, args, meta)

"""
    walk(visit, x)

Read-only pre-order traversal. `visit(node)::Bool` runs on each `BasicSymbolic` node and
returns `false` to stop descending into that node. Non-symbolic nodes are skipped.
"""
walk(visit, x::Symbolics.Num) = walk(visit, SymbolicUtils.unwrap(x))
function walk(visit, x)
    x isa SymbolicUtils.BasicSymbolic || return nothing
    visit(x) || return nothing
    SymbolicUtils.iscall(x) || return nothing
    for a in SymbolicUtils.arguments(x)
        walk(visit, a)
    end
    return nothing
end

"""
Collect the leaf averages of `x`, left to right, into a `Vector` (with multiplicity). A
lazy iterator is unnecessary at the bounded build-time sizes seen here.
"""
eachleaf(x::Symbolics.Num) = eachleaf(SymbolicUtils.unwrap(x))
function eachleaf(x)
    out = SymbolicUtils.BasicSymbolic[]
    walk(x) do n
        _is_avg_leaf(n) ? (push!(out, n); false) : true
    end
    return out
end

"""
Rewrite each leaf average of `x` via `f`, rebuilding only the branches that actually
changed.
"""
mapleaves(f, x::Symbolics.Num) = Symbolics.wrap(mapleaves(f, SymbolicUtils.unwrap(x)))
mapleaves(f, x) = rewrite(y -> _is_avg_leaf(y) ? f(y) : nothing, x)

"""
Replace every subtree of `x` that is a key of `sub` (descending through the whole tree)
with its mapped value, rebuilding the enclosing calls.
"""
_subtree_substitute(x, sub) = rewrite(y -> get(sub, y, nothing), x; descend = SymbolicUtils.iscall)
