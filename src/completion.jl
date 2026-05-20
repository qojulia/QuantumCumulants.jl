"""
    find_missing(eqs::AbstractMeanFieldEquations; filter_func=nothing,
                 get_adjoints=true)

Return the unique averages appearing on RHSes of `eqs.equations` whose
identity is not already present in `eqs.states`.

An average ⟨X⟩ and its conjugate ⟨X†⟩ are treated as the same state for
the purpose of "already covered" detection — if either is a state, neither
counts as missing. When `get_adjoints=false`, only one of each
conjugate-pair is emitted into `missing_states` (the conjugate is marked
seen so it isn't returned again on a later RHS scan).
"""
function find_missing(
        eqs::AbstractMeanFieldEquations; filter_func = nothing,
        get_adjoints::Bool = true
    )
    canon = _build_canonical_indices(eqs)
    seen_keys = Set{SQA.QAdd}()
    for s in eqs.states
        push!(seen_keys, _canonical_key(s, canon))
        push!(seen_keys, _canonical_key(_avg_conj(s), canon))
    end
    missing_states = SymbolicUtils.BasicSymbolic[]
    for eq in eqs.equations
        _collect_missing!(missing_states, seen_keys, canon, eq.rhs, get_adjoints)
    end
    # For noise equations, scan the stochastic RHS too so missing states
    # introduced only by the dW/dt term get derived as well.
    if eqs isa NoiseMeanFieldEquations
        for eq in eqs.noise_equations
            _collect_missing!(missing_states, seen_keys, canon, eq.rhs, get_adjoints)
        end
    end
    if filter_func !== nothing
        filter!(filter_func, missing_states)
    end
    return missing_states
end

function _collect_missing!(
        missing_states, seen_keys, canon, x, get_adjoints::Bool = true
    )
    x isa SymbolicUtils.BasicSymbolic || return
    if _is_leaf_average(x)
        key = _canonical_key(x, canon)
        key in seen_keys && return
        conj_key = _canonical_key(_avg_conj(x), canon)
        push!(seen_keys, key)
        push!(seen_keys, conj_key)
        # Push canonical-form averages: the missing set must be invariant
        # of RHS scan order, so two passes that encounter different members
        # of a conjugate pair first still seed identical states.
        push!(missing_states, average(key))
        if get_adjoints && conj_key != key
            push!(missing_states, average(conj_key))
        end
        return
    end
    SymbolicUtils.iscall(x) || return
    for a in SymbolicUtils.arguments(x)
        _collect_missing!(missing_states, seen_keys, canon, a, get_adjoints)
    end
    return
end

# Two leaf averages are alpha-equivalent, and therefore the same physical
# state, when their operator structures coincide after renaming free or
# bound indices on each Hilbert subspace to a deterministic canonical name.
# The canonical name comes from the user's own index vocabulary: for each
# subspace, the indices the user constructed in declaration order. State
# identity never depends on a name the algebra invented.
const _CanonIndex = Dict{Int, Vector{SQA.Index}}

function _build_canonical_indices(eqs::AbstractMeanFieldEquations)
    canon = _CanonIndex()
    # Collective-decay `jumps` / `jumps_dagger` may be nested vectors; flatten
    # one level so each source is a single `QField` SQA can call `get_indices`
    # on.
    flat_jumps = _flatten_jumps(eqs.jumps)
    flat_jdag = _flatten_jumps(eqs.jumps_dagger)
    sources = Iterators.flatten(
        (
            (eqs.hamiltonian,), eqs.operators, flat_jumps, flat_jdag,
        )
    )
    for src in sources, idx in SQA.get_indices(src)
        v = get!(canon, idx.space_index, SQA.Index[])
        idx in v || push!(v, idx)
    end
    for v in values(canon)
        sort!(v, by = idx -> idx.name)
    end
    return canon
end

_flatten_jumps(js::AbstractVector{<:QField}) = js
function _flatten_jumps(js::AbstractVector)
    isempty(js) && return QField[]
    eltype(js) <: AbstractVector || return js
    out = QField[]
    for jk in js
        append!(out, jk)
    end
    return out
end

function _canonical_key(x::SymbolicUtils.BasicSymbolic, canon::_CanonIndex)
    op_raw = SQA.undo_average(x)
    op_raw isa SQA.QAdd || return op_raw
    # Strip sum-scope `.indices` for the dedup key: states are stored in the
    # per-atom template form (see `_undo_for_derivation`), so a sum-scoped RHS
    # leaf `Σ_i ⟨σ_{i,22}⟩` must dedup against the per-atom state
    # `⟨σ_{i,22}⟩`. `evaluate` later re-materialises the sum at codegen by
    # enumerating the concrete atom states.
    op = isempty(op_raw.indices) ? op_raw :
        SQA.QAdd(op_raw.arguments, SQA.Index[])
    encountered = _free_op_indices(op)
    pos_by_space = Dict{Int, Int}()
    rename = Dict{SQA.Index, SQA.Index}()
    for idx in encountered
        pos = get(pos_by_space, idx.space_index, 0) + 1
        pos_by_space[idx.space_index] = pos
        space_canon = get(canon, idx.space_index, nothing)
        space_canon === nothing && continue
        pos <= length(space_canon) || continue
        target = space_canon[pos]
        target == idx && continue
        rename[idx] = target
    end
    if isempty(rename)
        return op === op_raw ? op : SQA.QAdd(op.arguments, SQA.Index[])
    end
    # Batched `change_index` does the rename simultaneously, so a target
    # name that coincides with another encountered index (e.g. swapping
    # `k <-> j` under canon `[j, k]`) does not fuse them destructively
    # mid-rename.
    result = SQA.change_index(op, rename)
    return SQA.QAdd(result.arguments, SQA.Index[])
end

# Collect distinct free indices appearing in operators inside `op`, in
# encounter order. Used by both `_canonical_key` (alpha-equivalence) and
# `_scale_qadd` (permutation collapse).
function _free_op_indices(op::SQA.QAdd)
    out = SQA.Index[]
    for (term, _) in op.arguments, o in term.ops
        SQA.has_index(o.index) || continue
        o.index in out || push!(out, o.index)
    end
    return out
end

# Drop sum-scope `.indices` so the derived operator represents the
# index-parametrized family rather than a sum-wrapped scalar.
function _undo_for_derivation(m::SymbolicUtils.BasicSymbolic)
    op = SQA.undo_average(m)
    op isa SQA.QAdd || return op
    return SQA.QAdd(
        op.arguments, SQA.Index[]
    )
end

# Conjugate of a leaf average: ⟨op⟩ ↦ ⟨op'⟩. We build it via SQA's `average`
# constructor rather than calling `conj` so the result lands in the same
# canonical form (sym_average head, same hashing) as states pushed by
# `meanfield`.
function _avg_conj(x::SymbolicUtils.BasicSymbolic)
    _is_leaf_average(x) || return x
    op = SQA.undo_average(x)
    return average(adjoint(op))
end

# `SQA.is_average` returns true for the whole AvgSym
# expression tree (e.g. ⟨a⟩*⟨σ₂₂⟩, which is a product of averages, still has
# symtype === AvgSym). For closure detection we want only *leaf* averages,
# the ones produced by `average(op)` directly, where the head is `sym_average`.
function _is_leaf_average(x::SymbolicUtils.BasicSymbolic)
    SQA.is_average(x) || return false
    SymbolicUtils.iscall(x) || return false
    return SymbolicUtils.operation(x) === SQA.sym_average
end
_is_leaf_average(::Any) = false

"""
    complete!(eqs::AbstractMeanFieldEquations; max_iter=200, simplify=true,
              filter_func=nothing, mix_choice=maximum)

Iteratively derive equations for missing averages until the system is closed
(or `max_iter` is reached). Mutates `eqs` in place.

The order used for cumulant expansion is the one already stored in
`eqs.order` (set at `meanfield(...; order=...)` time, or via
`cumulant_expansion(eqs, order)`). The non-mutating [`complete`](@ref)
variant additionally accepts an `order` keyword.
"""
function complete!(
        eqs::MeanFieldEquations; max_iter::Int = 200,
        simplify::Bool = true, filter_func = nothing,
        mix_choice = maximum, get_adjoints::Bool = true
    )
    for _ in 1:max_iter
        missing_states = find_missing(eqs; filter_func, get_adjoints)
        if isempty(missing_states)
            filter_func !== nothing && _filter_rhs!(eqs, filter_func)
            return eqs
        end
        new_ops = QField[_undo_for_derivation(m) for m in missing_states]
        new_eqs = _derive_for(eqs, new_ops; simplify, mix_choice)
        if filter_func !== nothing
            _filter_rhs!(new_eqs, filter_func)
        end
        _append!(eqs, new_eqs)
    end
    error("complete!: did not close within $max_iter iterations")
end

function _filter_rhs!(eqs::MeanFieldEquations, filter_func)
    for (i, eq) in enumerate(eqs.equations)
        new_rhs = _filter_expr(eq.rhs, filter_func)
        eqs.equations[i] = eq.lhs ~ new_rhs
    end
    return eqs
end

function _filter_rhs!(eqs::NoiseMeanFieldEquations, filter_func)
    for (i, eq) in enumerate(eqs.equations)
        new_rhs = _filter_expr(eq.rhs, filter_func)
        eqs.equations[i] = eq.lhs ~ new_rhs
    end
    for (i, eq) in enumerate(eqs.noise_equations)
        new_rhs = _filter_expr(eq.rhs, filter_func)
        eqs.noise_equations[i] = eq.lhs ~ new_rhs
    end
    return eqs
end

function _filter_expr(x, filter_func)
    if x isa SymbolicUtils.BasicSymbolic
        if SQA.is_average(x)
            return filter_func(x) ? x : 0
        end
        if SymbolicUtils.iscall(x) && _has_average(x)
            op = SymbolicUtils.operation(x)
            args = SymbolicUtils.arguments(x)
            new_args = Any[_filter_expr(a, filter_func) for a in args]
            return op(new_args...)
        end
    end
    return x
end

"""
    complete(eqs::AbstractMeanFieldEquations; order=nothing, kw...)

Non-mutating variant. Returns a new `MeanFieldEquations`. If `order` is
given, the returned system is first cumulant-expanded to that order before
the completion loop runs.
"""
function MTK.complete(
        eqs::MeanFieldEquations; order = nothing,
        mix_choice = maximum, simplify::Bool = true, kw...
    )
    eqs_copy = order === nothing ?
        _copy(eqs) :
        cumulant_expansion(_copy(eqs), order; simplify, mix_choice)
    return complete!(eqs_copy; simplify, mix_choice, kw...)
end

function _derive_for(
        eqs::MeanFieldEquations, new_ops;
        simplify::Bool, mix_choice = maximum
    )
    return _meanfield_deterministic(
        eqs.direction, new_ops, eqs.hamiltonian, eqs.jumps,
        eqs.jumps_dagger, eqs.rates, eqs.order,
        simplify, mix_choice, eqs.iv,
    )
end

function _derive_for(
        eqs::NoiseMeanFieldEquations, new_ops;
        simplify::Bool, mix_choice = maximum
    )
    return _meanfield_noise(
        eqs.direction, new_ops, eqs.hamiltonian, eqs.jumps,
        eqs.jumps_dagger, eqs.rates, eqs.efficiencies,
        eqs.order, simplify, mix_choice, eqs.iv
    )
end

function complete!(
        eqs::NoiseMeanFieldEquations; max_iter::Int = 200,
        simplify::Bool = true, filter_func = nothing,
        mix_choice = maximum, get_adjoints::Bool = true
    )
    for _ in 1:max_iter
        missing_states = find_missing(eqs; filter_func, get_adjoints)
        if isempty(missing_states)
            filter_func !== nothing && _filter_rhs!(eqs, filter_func)
            return eqs
        end
        new_ops = QField[_undo_for_derivation(m) for m in missing_states]
        new_eqs = _derive_for(eqs, new_ops; simplify, mix_choice)
        if filter_func !== nothing
            _filter_rhs!(new_eqs, filter_func)
        end
        _append!(eqs, new_eqs)
    end
    error("complete!: did not close within $max_iter iterations")
end

function MTK.complete(
        eqs::NoiseMeanFieldEquations; order = nothing,
        mix_choice = maximum, simplify::Bool = true, kw...
    )
    order === nothing || throw(ArgumentError(
        "`complete(::NoiseMeanFieldEquations; order=...)` is not supported; pass `order` to `meanfield` instead."
    ))
    eqs_copy = _copy(eqs)
    return complete!(eqs_copy; simplify, mix_choice, kw...)
end

function _append!(a::MeanFieldEquations, b::MeanFieldEquations)
    append!(a.equations, b.equations)
    append!(a.operator_equations, b.operator_equations)
    append!(a.states, b.states)
    append!(a.operators, b.operators)
    return a
end

function _append!(a::NoiseMeanFieldEquations, b::NoiseMeanFieldEquations)
    append!(a.equations, b.equations)
    append!(a.noise_equations, b.noise_equations)
    append!(a.operator_equations, b.operator_equations)
    append!(a.operator_noise_equations, b.operator_noise_equations)
    append!(a.states, b.states)
    append!(a.operators, b.operators)
    return a
end

function _copy(eqs::MeanFieldEquations)
    return MeanFieldEquations(
        copy(eqs.equations),
        copy(eqs.operator_equations),
        copy(eqs.states),
        copy(eqs.operators),
        eqs.hamiltonian,
        copy(eqs.jumps),
        copy(eqs.jumps_dagger),
        copy(eqs.rates),
        eqs.iv, eqs.order, eqs.direction,
    )
end

function _copy(eqs::NoiseMeanFieldEquations)
    return NoiseMeanFieldEquations(
        copy(eqs.equations),
        copy(eqs.noise_equations),
        copy(eqs.operator_equations),
        copy(eqs.operator_noise_equations),
        copy(eqs.states),
        copy(eqs.operators),
        eqs.hamiltonian,
        copy(eqs.jumps),
        copy(eqs.jumps_dagger),
        copy(eqs.rates),
        copy(eqs.efficiencies),
        eqs.iv, eqs.order, eqs.direction,
    )
end
