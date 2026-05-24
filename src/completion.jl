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

function _free_indices_by_space(op::SQA.QAdd)
    out = Dict{Int, Vector{SQA.Index}}()
    for (term, _) in op.arguments, o in term.ops
        SQA.has_index(o.index) || continue
        v = get!(out, o.index.space_index, SQA.Index[])
        o.index in v || push!(v, o.index)
    end
    return out
end

# Predicate: should Hilbert subspace `sp` participate in scaling? Empty
# `h_set` means "all subspaces" (the default), matching `scale(; h=Int[])`.
_in_h(h_set::Set{Int}, sp::Int) = isempty(h_set) || sp in h_set

function _saturate_and_canonicalize(op::SQA.QAdd, h_set::Set{Int})
    out = SQA.QTermDict()
    for (term, c) in op.arguments
        ne = _maximal_pairwise_ne(term.ops, term.ne, h_set)
        SQA._canonicalize!(out, copy(term.ops), c, ne)
    end
    return SQA.QAdd(out, op.indices)
end

function _maximal_pairwise_ne(
        ops::Vector{<:SQA.QSym},
        existing::Vector{SQA.NonEqualPair}, h_set::Set{Int},
    )
    by_space = Dict{Int, Vector{SQA.Index}}()
    for o in ops
        SQA.has_index(o.index) || continue
        _in_h(h_set, o.index.space_index) || continue
        v = get!(by_space, o.index.space_index, SQA.Index[])
        o.index in v || push!(v, o.index)
    end
    out = copy(existing)
    seen = Set{Tuple{SQA.Index, SQA.Index}}()
    for (α, β) in existing
        push!(seen, (α, β))
        push!(seen, (β, α))
    end
    for (_, idxs) in by_space, a in idxs, b in idxs
        a < b || continue
        (a, b) in seen && continue
        push!(out, (a, b))
        push!(seen, (a, b))
        push!(seen, (b, a))
    end
    return out
end

# Strip NE entries among same-space indices on selected subspaces (the
# canonical slot order has absorbed them); preserve NE that touches an
# unselected subspace or crosses subspaces.
function _strip_all_ne(op::SQA.QAdd, h_set::Set{Int})
    out = SQA.QTermDict()
    for (term, c) in op.arguments
        kept = SQA.NonEqualPair[]
        for pair in term.ne
            sp = pair[1].space_index
            sp == pair[2].space_index || (push!(kept, pair); continue)
            _in_h(h_set, sp) && continue
            push!(kept, pair)
        end
        new_term = length(kept) == length(term.ne) ? term :
            SQA.QTerm(copy(term.ops), kept)
        out[new_term] = c
    end
    return SQA.QAdd(out, op.indices)
end

# Enumerate slot permutations (Cartesian product across Hilbert subspaces)
# and pick the lexicographically smallest rendered form. Each space's K
# free indices map to its K slot reps (slot 1 = user's first-declared
# index; slot k > 1 minted as `Symbol(first.name, "_", k)`). For two-atom
# correlations under permutation symmetry, this identifies
# `⟨σ_a^{12} σ_b^{21}⟩` with `⟨σ_a^{21} σ_b^{12}⟩` (same state under
# `a ↔ b` swap), matching master's numeric-label canonicalization.
function _min_slot_assignment(
        op::SQA.QAdd, free_by_space::Dict{Int, Vector{SQA.Index}},
        canon::_CanonIndex, h_set::Set{Int},
    )
    spaces = sort!(collect(keys(free_by_space)))
    slot_reps = Dict{Int, Vector{SQA.Index}}()
    for sp in spaces
        K = length(free_by_space[sp])
        canon_list = get(canon, sp, SQA.Index[])
        first_idx = isempty(canon_list) ? first(free_by_space[sp]) : first(canon_list)
        slot_reps[sp] = SQA.Index[first_idx(k) for k in 1:K]
    end
    perms_per_space = [collect(Combinatorics.permutations(free_by_space[sp])) for sp in spaces]
    best_op = nothing
    best_key = nothing
    for combo in Iterators.product(perms_per_space...)
        sub = Dict{SQA.Index, SQA.Index}()
        for (i, sp) in enumerate(spaces)
            perm = combo[i]
            for (k, idx) in enumerate(perm)
                target = slot_reps[sp][k]
                target == idx && continue
                sub[idx] = target
            end
        end
        candidate = isempty(sub) ? op : SQA.change_index(op, sub)
        candidate = _strip_all_ne(candidate, h_set)
        key = _serialize_for_compare(candidate)
        if best_key === nothing || key < best_key
            best_key = key
            best_op = candidate
        end
    end
    return best_op === nothing ? _strip_all_ne(op, h_set) : best_op
end

function _serialize_for_compare(op::SQA.QAdd)
    terms = Vector{Any}()
    for (term, c) in op.arguments
        push!(terms, (String[string(o) for o in term.ops], string(c)))
    end
    sort!(terms)
    return terms
end

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
    complete!(eqs::AbstractMeanFieldEquations; max_iter=200,
              filter_func=nothing, mix_choice=maximum)

Iteratively derive equations for missing averages until the system is closed
(or `max_iter` is reached). Mutates `eqs` in place. The newly derived RHS
expressions are left in their raw, unsimplified form; apply
`SymbolicUtils.simplify` yourself if you want a canonical representation.

The order used for cumulant expansion is the one already stored in
`eqs.order` (set at `meanfield(...; order=...)` time, or via
`cumulant_expansion(eqs, order)`). The non-mutating [`complete`](@ref)
variant additionally accepts an `order` keyword.
"""
function complete!(
        eqs::MeanFieldEquations; max_iter::Int = 200,
        filter_func = nothing,
        mix_choice = maximum, get_adjoints::Bool = true
    )
    for _ in 1:max_iter
        missing_states = find_missing(eqs; filter_func, get_adjoints)
        if isempty(missing_states)
            filter_func !== nothing && _filter_rhs!(eqs, filter_func)
            return eqs
        end
        new_ops = QField[_undo_for_derivation(m) for m in missing_states]
        new_eqs = _derive_for(eqs, new_ops; mix_choice)
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
        mix_choice = maximum, kw...
    )
    eqs_copy = order === nothing ?
        _copy(eqs) :
        cumulant_expansion(_copy(eqs), order; mix_choice)
    return complete!(eqs_copy; mix_choice, kw...)
end

function _derive_for(
        eqs::MeanFieldEquations, new_ops; mix_choice = maximum
    )
    return _meanfield_deterministic(
        eqs.direction, new_ops, eqs.hamiltonian, eqs.jumps, eqs.jumps_dagger,
        eqs.rates, eqs.order, mix_choice, eqs.iv,
    )
end

function _derive_for(
        eqs::NoiseMeanFieldEquations, new_ops; mix_choice = maximum
    )
    return _meanfield_noise(
        eqs.direction, new_ops, eqs.hamiltonian, eqs.jumps, eqs.jumps_dagger,
        eqs.rates, eqs.efficiencies, eqs.order, mix_choice, eqs.iv,
    )
end

function complete!(
        eqs::NoiseMeanFieldEquations; max_iter::Int = 200,
        filter_func = nothing,
        mix_choice = maximum, get_adjoints::Bool = true
    )
    for _ in 1:max_iter
        missing_states = find_missing(eqs; filter_func, get_adjoints)
        if isempty(missing_states)
            filter_func !== nothing && _filter_rhs!(eqs, filter_func)
            return eqs
        end
        new_ops = QField[_undo_for_derivation(m) for m in missing_states]
        new_eqs = _derive_for(eqs, new_ops; mix_choice)
        if filter_func !== nothing
            _filter_rhs!(new_eqs, filter_func)
        end
        _append!(eqs, new_eqs)
    end
    error("complete!: did not close within $max_iter iterations")
end

function MTK.complete(
        eqs::NoiseMeanFieldEquations; order = nothing,
        mix_choice = maximum, kw...
    )
    order === nothing || throw(
        ArgumentError(
            "`complete(::NoiseMeanFieldEquations; order=...)` is not supported; pass `order` to `meanfield` instead."
        )
    )
    eqs_copy = _copy(eqs)
    return complete!(eqs_copy; mix_choice, kw...)
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
