"""
    scale!(eqs::MeanFieldEquations)
    scale(eqs::MeanFieldEquations)

Permutation-symmetry collapse for indexed mean-field systems.

For every free `Index` on a Hilbert subspace where the user declared multiple
indices (typical pattern: `i, j ∈ atom space`), rename it to that subspace's
canonical-first index. Two operators that differ only by which symmetric
atom they reference therefore collapse to one. Sum-scope `.indices` metadata
on an averaged `QAdd` collapses to a multiplicative prefactor
`(range − |constraint pairs|)` only when the bound index actually appears
in the operator (a spurious annotation otherwise; see SQA `_accumulate_with_diag!`).

`IndexedVariable(:g, i)` flattens to a scalar `Num g` (uniform coupling
under the scale symmetry).

The implementation is a thin layer on SQA primitives: `change_index` does
the substitution (and routes through `_canonicalize!`), `average` /
`undo_average` round-trip across the operator-vs-symbolic boundary, and
SQA's pipeline handles same-site collapse and NE-violated-term drop.
"""
function scale!(eqs::MeanFieldEquations)
    isempty(eqs.equations) && return eqs
    canon = _build_canonical_indices(eqs)
    new_eqs, new_states, new_ops, new_op_eqs = _do_scale(eqs, canon)
    empty!(eqs.equations);          append!(eqs.equations, new_eqs)
    empty!(eqs.operator_equations); append!(eqs.operator_equations, new_op_eqs)
    empty!(eqs.states);              append!(eqs.states, new_states)
    empty!(eqs.operators);           append!(eqs.operators, new_ops)
    return eqs
end

scale(eqs::MeanFieldEquations) = scale!(_copy(eqs))

"""
    scale!(eqs::NoiseMeanFieldEquations)
    scale(eqs::NoiseMeanFieldEquations)

Same as `scale!(::MeanFieldEquations)` but also rewrites
`eqs.noise_equations` (and the operator-form mirrors) under the same
canonicalisation.
"""
function scale!(eqs::NoiseMeanFieldEquations)
    isempty(eqs.equations) && return eqs
    canon = _build_canonical_indices(eqs)
    new_eqs, new_states, new_ops, new_op_eqs = _do_scale(eqs, canon)
    # Noise channel
    new_noise_eqs = Symbolics.Equation[]
    new_op_noise_eqs = Symbolics.Equation[]
    seen = Set{Any}()
    for (k, eq) in enumerate(eqs.equations)
        new_lhs = _scale_expr(eq.lhs, canon)
        key = _scale_state_key(new_lhs)
        key in seen && continue
        push!(seen, key)
        new_nrhs = _scale_expr(eqs.noise_equations[k].rhs, canon)
        push!(new_noise_eqs, new_lhs ~ new_nrhs)
        push!(new_op_noise_eqs, eqs.operator_noise_equations[k])
    end
    empty!(eqs.equations);                 append!(eqs.equations, new_eqs)
    empty!(eqs.noise_equations);           append!(eqs.noise_equations, new_noise_eqs)
    empty!(eqs.operator_equations);        append!(eqs.operator_equations, new_op_eqs)
    empty!(eqs.operator_noise_equations);  append!(eqs.operator_noise_equations, new_op_noise_eqs)
    empty!(eqs.states);                    append!(eqs.states, new_states)
    empty!(eqs.operators);                 append!(eqs.operators, new_ops)
    return eqs
end

scale(eqs::NoiseMeanFieldEquations) = scale!(_copy(eqs))

# Shared core: rewrite each equation, deduplicate by state key.
function _do_scale(eqs, canon)
    new_eqs = Symbolics.Equation[]
    new_states = SymbolicUtils.BasicSymbolic[]
    new_ops = QAdd[]
    new_op_eqs = Symbolics.Equation[]
    seen = Set{Any}()
    for k in eachindex(eqs.equations)
        new_lhs = _scale_expr(eqs.equations[k].lhs, canon)
        # If scaling reduced the LHS to a tautological zero (e.g. an average
        # over an empty NE-constrained sum), drop the equation entirely.
        # `0 ~ rhs` is not a useful state evolution and would crash
        # `_as_average`.
        _is_lhs_zero(new_lhs) && continue
        new_rhs = _scale_expr(eqs.equations[k].rhs, canon)
        new_lhs_avg = _as_average(new_lhs)
        key = _scale_state_key(new_lhs_avg)
        key in seen && continue
        push!(seen, key)
        push!(new_eqs, new_lhs_avg ~ new_rhs)
        push!(new_states, new_lhs_avg)
        push!(new_ops, _scale_qadd(eqs.operators[k], canon))
        push!(new_op_eqs, eqs.operator_equations[k])
    end
    return new_eqs, new_states, new_ops, new_op_eqs
end

# Tree walker. At averaged leaves: rename free indices + collapse sum scope.
# At arithmetic nodes: recurse and rebuild.
function _scale_expr(x, canon)
    x isa SymbolicUtils.BasicSymbolic || return x
    if _is_leaf_average(x)
        return _scale_avg(x, canon)
    end
    if _is_indexed_var(x)
        return _flatten_indexed_var(x)
    end
    SymbolicUtils.iscall(x) || return x
    op = SymbolicUtils.operation(x)
    args = SymbolicUtils.arguments(x)
    new_args = Any[_scale_expr(a, canon) for a in args]
    all(((a, b),) -> isequal(a, b), zip(args, new_args)) && return x
    op === complex && length(new_args) == 2 &&
        return new_args[1] + new_args[2] * Symbolics.IM
    try
        return op(new_args...)
    catch err
        err isa MethodError || err isa ArgumentError || rethrow()
        return TermInterface.maketerm(typeof(x), op, new_args, TermInterface.metadata(x))
    end
end

# Average-leaf rewrite: prefactor from sum-scope (when bound index actually
# appears in op), then strip-and-rename, then re-wrap.
function _scale_avg(avg::SymbolicUtils.BasicSymbolic, canon)
    op = SQA.undo_average(avg)
    op isa QAdd || return avg
    prefactor = _sum_scope_prefactor(op)
    new_op = _scale_qadd(op, canon)
    result = prefactor === 1 ? average(new_op) : prefactor * average(new_op)
    return _flatten_indexed_vars_in_tree(result)
end

# Operator-level rewrite: strip sum scope, then for each free Index
# appearing in operators, rename via SQA.change_index to canonical-first
# for its Hilbert subspace. SQA's _canonicalize! (invoked by change_index
# when .indices is empty) handles same-name collapse, ne-violated drop,
# and projector squashing.
function _scale_qadd(op::QAdd, canon)
    stripped = QAdd(op.arguments, SQA.Index[])
    free = _free_op_indices(stripped)
    result = stripped
    for idx in free
        target = _canon_first(canon, idx.space_index)
        target === nothing && continue
        idx == target && continue
        result = SQA.change_index(result, idx, target)
    end
    return result
end
_scale_qadd(op::SQA.QSym, canon) = op
_scale_qadd(op, _) = op

_canon_first(canon, space_index) = begin
    list = get(canon, space_index, nothing)
    (list === nothing || isempty(list)) ? nothing : first(list)
end

# Sum-scope prefactor: for each bound index `b` in `op.indices` that
# *actually* appears in some op.term.ops, contribute `(b.range - count_ne)`.
# Spurious bound indices (range scope from accumulation, not used in ops)
# contribute 1 because the sum was already eagerly scalarised by SQA.
function _sum_scope_prefactor(op::QAdd)
    isempty(op.indices) && return 1
    op_indices = Set{SQA.Index}()
    for (term, _) in op.arguments, o in term.ops
        SQA.has_index(o.index) && push!(op_indices, o.index)
    end
    prefactor = 1
    used_ne = Set{Int}()
    for b in op.indices
        b in op_indices || continue
        count_b = 0
        for (term, _) in op.arguments, (k, pair) in enumerate(term.ne)
            k in used_ne && continue
            (pair[1] == b || pair[2] == b) || continue
            count_b += 1
            push!(used_ne, k)
        end
        factor = b.range - count_b
        prefactor = prefactor isa Number && prefactor == 1 ? factor : prefactor * factor
    end
    return prefactor
end

# Dedup key: the underlying operator of the scaled state. After scaling,
# free indices all live at canonical-first names, so op-equality is enough.
function _scale_state_key(avg::SymbolicUtils.BasicSymbolic)
    op = SQA.undo_average(avg)
    return op
end

# Unwrap a 1*avg product if `_scale_expr` produced one.
_is_lhs_zero(x::Number) = iszero(x)
function _is_lhs_zero(x::SymbolicUtils.BasicSymbolic)
    SymbolicUtils.isconst(x) && return iszero(x.val)
    return false
end
_is_lhs_zero(_) = false

function _as_average(x)
    x isa SymbolicUtils.BasicSymbolic ||
        error("scale: scaled LHS must be a BasicSymbolic, got $(typeof(x))")
    _is_leaf_average(x) && return x
    if SymbolicUtils.iscall(x) && SymbolicUtils.operation(x) === *
        args = SymbolicUtils.arguments(x)
        nonone = filter(a -> !(a isa Number && a == 1), args)
        length(nonone) == 1 && return _as_average(nonone[1])
    end
    error("scale: LHS did not reduce to an average: $x")
end

# IndexedVariable detection and flattening (`g(i)` -> scalar Num `g`).
# SQA materialises `IndexedVariable(:g, i)` as a Term whose `operation` is a
# Sym{SymReal} of FnType{Tuple{Int}, Real, Nothing}.
function _is_indexed_var(x)
    x isa SymbolicUtils.BasicSymbolic || return false
    SymbolicUtils.iscall(x) || return false
    op = SymbolicUtils.operation(x)
    op isa SymbolicUtils.BasicSymbolic || return false
    st = SymbolicUtils.symtype(op)
    return st <: SymbolicUtils.FnType && _fn_returns_real(st)
end

_fn_returns_real(::Type{<:SymbolicUtils.FnType{A, R}}) where {A, R} = R <: Real
_fn_returns_real(::Type) = false

function _flatten_indexed_var(x::SymbolicUtils.BasicSymbolic)
    op = SymbolicUtils.operation(x)
    name = nameof(op)
    return SymbolicUtils.Sym{SymbolicUtils.SymReal}(name; type = Real)
end

function _flatten_indexed_vars_in_tree(x)
    x isa SymbolicUtils.BasicSymbolic || return x
    _is_indexed_var(x) && return _flatten_indexed_var(x)
    SymbolicUtils.iscall(x) || return x
    op = SymbolicUtils.operation(x)
    op === SQA.sym_average && return x
    args = SymbolicUtils.arguments(x)
    new_args = Any[_flatten_indexed_vars_in_tree(a) for a in args]
    all(((a, b),) -> isequal(a, b), zip(args, new_args)) && return x
    try
        return op(new_args...)
    catch err
        err isa MethodError || err isa ArgumentError || rethrow()
        return TermInterface.maketerm(typeof(x), op, new_args, TermInterface.metadata(x))
    end
end

# Substitute a coefficient sub_dict into a BasicSymbolic without breaking on
# Σ-shaped nodes. Used by evaluate.jl; lives here for legacy reasons.
function _safe_substitute(x, sub)
    isempty(sub) && return x
    x isa SymbolicUtils.BasicSymbolic || return x
    if haskey(sub, x)
        return sub[x]
    end
    SymbolicUtils.iscall(x) || return x
    op = SymbolicUtils.operation(x)
    op === SQA.sym_average && return x
    args = SymbolicUtils.arguments(x)
    new_args = Any[_safe_substitute(a, sub) for a in args]
    all(((a, b),) -> a === b, zip(args, new_args)) && return x
    op === complex && length(new_args) == 2 &&
        return new_args[1] + new_args[2] * Symbolics.IM
    try
        return op(new_args...)
    catch err
        err isa MethodError || err isa ArgumentError || rethrow()
        return TermInterface.maketerm(typeof(x), op, new_args, TermInterface.metadata(x))
    end
end
