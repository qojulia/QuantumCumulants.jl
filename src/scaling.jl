"""
    scale!(eqs::MeanFieldEquations; h::Vector{Int} = Int[])
    scale(eqs::MeanFieldEquations; h::Vector{Int} = Int[])

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

`h` selects which Hilbert subspaces participate in scaling, identified by
`space_index` (the 1-based tensor-product position). The empty default
scales every subspace. Indices on unselected subspaces keep their symbolic
names and sum-scope, and `IndexedVariable`s whose source index sits on an
unselected subspace are not flattened. Useful for hybrid systems (e.g.
unroll a filter cavity array via `evaluate` while collapsing a
permutation-symmetric atom ensemble via `scale`).

The implementation is a thin layer on SQA primitives: `change_index` does
the substitution (and routes through `_canonicalize!`), `average` /
`undo_average` round-trip across the operator-vs-symbolic boundary, and
SQA's pipeline handles same-site collapse and NE-violated-term drop.
"""
function scale!(eqs::MeanFieldEquations; h::Vector{Int} = Int[])
    isempty(eqs.equations) && return eqs
    canon = _build_canonical_indices(eqs)
    h_set = Set{Int}(h)
    sym_to_space = _build_sym_to_space(canon)
    new_eqs, new_states, new_ops, new_op_eqs =
        _do_scale(eqs, canon, h_set, sym_to_space)
    empty!(eqs.equations);          append!(eqs.equations, new_eqs)
    empty!(eqs.operator_equations); append!(eqs.operator_equations, new_op_eqs)
    empty!(eqs.states);              append!(eqs.states, new_states)
    empty!(eqs.operators);           append!(eqs.operators, new_ops)
    return eqs
end

scale(eqs::MeanFieldEquations; h::Vector{Int} = Int[]) = scale!(_copy(eqs); h)

"""
    scale!(eqs::NoiseMeanFieldEquations)
    scale(eqs::NoiseMeanFieldEquations)

Same as `scale!(::MeanFieldEquations)` but also rewrites
`eqs.noise_equations` (and the operator-form mirrors) under the same
canonicalisation.
"""
function scale!(eqs::NoiseMeanFieldEquations; h::Vector{Int} = Int[])
    isempty(eqs.equations) && return eqs
    canon = _build_canonical_indices(eqs)
    h_set = Set{Int}(h)
    sym_to_space = _build_sym_to_space(canon)
    new_eqs, new_states, new_ops, new_op_eqs =
        _do_scale(eqs, canon, h_set, sym_to_space)
    # Noise channel
    new_noise_eqs = Symbolics.Equation[]
    new_op_noise_eqs = Symbolics.Equation[]
    seen = Set{Any}()
    for (k, eq) in enumerate(eqs.equations)
        new_lhs = _scale_expr(eq.lhs, canon, h_set, sym_to_space)
        key = _scale_state_key(new_lhs)
        key in seen && continue
        push!(seen, key)
        new_nrhs = _scale_expr(eqs.noise_equations[k].rhs, canon, h_set, sym_to_space)
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

scale(eqs::NoiseMeanFieldEquations; h::Vector{Int} = Int[]) = scale!(_copy(eqs); h)

# Predicate: should we scale Hilbert subspace `sp`?  Empty `h_set` means "all".
_in_h(h_set::Set{Int}, sp::Int) = isempty(h_set) || sp in h_set

# Build a Symbol -> space_index lookup from the canonical index registry,
# so `IndexedVariable(:g, i)` (whose Term carries only `i.sym`) can be traced
# back to its source Hilbert subspace at flatten time.
function _build_sym_to_space(canon)
    out = Dict{Symbol, Int}()
    for (sp, indices) in canon, idx in indices
        out[idx.name] = sp
    end
    return out
end

# Shared core: rewrite each equation, deduplicate by state key.
function _do_scale(eqs, canon, h_set::Set{Int}, sym_to_space::Dict{Symbol, Int})
    new_eqs = Symbolics.Equation[]
    new_states = SymbolicUtils.BasicSymbolic[]
    new_ops = QAdd[]
    new_op_eqs = Symbolics.Equation[]
    seen = Set{Any}()
    for k in eachindex(eqs.equations)
        new_lhs = _scale_expr(eqs.equations[k].lhs, canon, h_set, sym_to_space)
        # If scaling reduced the LHS to a tautological zero (e.g. an average
        # over an empty NE-constrained sum), drop the equation entirely.
        # `0 ~ rhs` is not a useful state evolution and would crash
        # `_as_average`.
        _is_lhs_zero(new_lhs) && continue
        new_rhs = _scale_expr(eqs.equations[k].rhs, canon, h_set, sym_to_space)
        new_lhs_avg = _as_average(new_lhs)
        key = _scale_state_key(new_lhs_avg)
        key in seen && continue
        push!(seen, key)
        push!(new_eqs, new_lhs_avg ~ new_rhs)
        push!(new_states, new_lhs_avg)
        push!(new_ops, _scale_qadd(eqs.operators[k], canon, h_set))
        push!(new_op_eqs, eqs.operator_equations[k])
    end
    return new_eqs, new_states, new_ops, new_op_eqs
end

# Tree walker. At averaged leaves: rename free indices + collapse sum scope.
# At arithmetic nodes: recurse and rebuild.
function _scale_expr(x, canon, h_set::Set{Int}, sym_to_space::Dict{Symbol, Int})
    x isa SymbolicUtils.BasicSymbolic || return x
    if _is_leaf_average(x)
        return _scale_avg(x, canon, h_set, sym_to_space)
    end
    if _is_indexed_var(x) && _indexed_var_in_h(x, h_set, sym_to_space)
        return _flatten_indexed_var(x)
    end
    SymbolicUtils.iscall(x) || return x
    op = SymbolicUtils.operation(x)
    args = SymbolicUtils.arguments(x)
    new_args = Any[_scale_expr(a, canon, h_set, sym_to_space) for a in args]
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
function _scale_avg(avg::SymbolicUtils.BasicSymbolic, canon, h_set, sym_to_space)
    op = SQA.undo_average(avg)
    op isa QAdd || return avg
    prefactor = _sum_scope_prefactor(op, h_set)
    new_op = _scale_qadd(op, canon, h_set)
    result = prefactor === 1 ? average(new_op) : prefactor * average(new_op)
    return _flatten_indexed_vars_in_tree(result, h_set, sym_to_space)
end

# Operator-level rewrite: strip sum scope (only for selected subspaces),
# then for each free Index appearing in operators whose Hilbert subspace is
# in `h_set`, rename via SQA.change_index to canonical-first for that
# subspace. SQA's _canonicalize! (invoked by change_index when .indices is
# empty) handles same-name collapse, ne-violated drop, and projector
# squashing. Indices on unselected subspaces keep both their names and
# their sum-scope entries in `.indices`.
function _scale_qadd(op::QAdd, canon, h_set::Set{Int})
    kept = SQA.Index[b for b in op.indices if !_in_h(h_set, b.space_index)]
    stripped = QAdd(op.arguments, kept)
    free = _free_op_indices(stripped)
    result = stripped
    for idx in free
        _in_h(h_set, idx.space_index) || continue
        target = _canon_first(canon, idx.space_index)
        target === nothing && continue
        idx == target && continue
        result = SQA.change_index(result, idx, target)
    end
    return result
end
_scale_qadd(op::SQA.QSym, canon, h_set::Set{Int}) = op
_scale_qadd(op, _, _) = op

_canon_first(canon, space_index) = begin
    list = get(canon, space_index, nothing)
    (list === nothing || isempty(list)) ? nothing : first(list)
end

# Sum-scope prefactor: for each bound index `b` in `op.indices` that
# *actually* appears in some op.term.ops and lives on a selected subspace,
# contribute `(b.range - count_ne)`. Spurious bound indices (range scope
# from accumulation, not used in ops) contribute 1 because the sum was
# already eagerly scalarised by SQA. Bound indices on unselected subspaces
# also contribute 1, because their sum stays symbolic (handled by
# `_scale_qadd`, which keeps them in `.indices`).
function _sum_scope_prefactor(op::QAdd, h_set::Set{Int})
    isempty(op.indices) && return 1
    op_indices = Set{SQA.Index}()
    for (term, _) in op.arguments, o in term.ops
        SQA.has_index(o.index) && push!(op_indices, o.index)
    end
    prefactor = 1
    used_ne = Set{Int}()
    for b in op.indices
        b in op_indices || continue
        _in_h(h_set, b.space_index) || continue
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
    # `c + avg` shape: SQA's `_canonicalize!` produces a constant offset when
    # a same-index commutator fires (e.g. `a_m * a_m' = a_m' * a_m + 1` after
    # renaming two distinct cavity indices to the canonical one). The state
    # `⟨a*a'⟩ = ⟨a'*a⟩ + 1` differs from `⟨a'*a⟩` by a constant only, so the
    # time derivative is the same; drop the additive constant here and let
    # the dedup key collapse the redundant equation against its normal-ordered
    # sibling. Numeric constant + (scaled) average leaves give the LHS;
    # ignore non-numeric coefficients (they're real RHS terms, not LHS shape).
    if SymbolicUtils.iscall(x) && SymbolicUtils.operation(x) === +
        args = SymbolicUtils.arguments(x)
        non_const = filter(a -> !_is_numeric(a), args)
        if length(non_const) == 1
            return _as_average(non_const[1])
        end
    end
    error("scale: LHS did not reduce to an average: $x")
end

_is_numeric(x::Number) = true
_is_numeric(x::SymbolicUtils.BasicSymbolic) = SymbolicUtils.isconst(x)
_is_numeric(_) = false

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

# True when the IndexedVariable's source index sits on a selected subspace.
# Unknown / unrecognised arg shapes fall back to "selected" so the default
# `h = []` flattens uniformly as before.
function _indexed_var_in_h(
        x::SymbolicUtils.BasicSymbolic, h_set::Set{Int},
        sym_to_space::Dict{Symbol, Int}
    )
    isempty(h_set) && return true
    SymbolicUtils.iscall(x) || return true
    args = SymbolicUtils.arguments(x)
    length(args) == 1 || return true
    a = args[1]
    a isa SymbolicUtils.BasicSymbolic || return true
    SymbolicUtils.iscall(a) && return true
    sp = get(sym_to_space, Base.nameof(a), 0)
    sp == 0 && return true
    return sp in h_set
end

function _flatten_indexed_vars_in_tree(
        x, h_set::Set{Int},
        sym_to_space::Dict{Symbol, Int}
    )
    x isa SymbolicUtils.BasicSymbolic || return x
    if _is_indexed_var(x) && _indexed_var_in_h(x, h_set, sym_to_space)
        return _flatten_indexed_var(x)
    end
    SymbolicUtils.iscall(x) || return x
    op = SymbolicUtils.operation(x)
    op === SQA.sym_average && return x
    args = SymbolicUtils.arguments(x)
    new_args = Any[_flatten_indexed_vars_in_tree(a, h_set, sym_to_space) for a in args]
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
