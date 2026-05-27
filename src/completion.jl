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
        get_adjoints::Bool = true,
        canon = nothing,
        bound = nothing,
    )
    canon = canon === nothing ? _build_canonical_indices(eqs) : canon
    bound = bound === nothing ? _bound_indices(eqs) : bound
    seen_keys = Set{Any}()
    for s in eqs.states
        push!(seen_keys, _dedup_key_strip_free_ne(_canonical_key(s, canon), bound))
        push!(seen_keys, _dedup_key_strip_free_ne(_canonical_key(_avg_conj(s), canon), bound))
    end
    missing_states = SymbolicUtils.BasicSymbolic[]
    for eq in eqs.equations
        _collect_missing!(
            missing_states, seen_keys, canon, bound, eq.rhs, get_adjoints,
        )
    end
    if eqs isa NoiseMeanFieldEquations
        for eq in eqs.noise_equations
            _collect_missing!(
                missing_states, seen_keys, canon, bound, eq.rhs, get_adjoints,
            )
        end
    end
    if filter_func !== nothing
        filter!(filter_func, missing_states)
    end
    return missing_states
end

function _collect_missing!(
        missing_states, seen_keys, canon, bound, x, get_adjoints::Bool = true,
    )
    x isa SymbolicUtils.BasicSymbolic || return
    if _is_leaf_average(x)
        key_full = _canonical_key(x, canon)
        dedup = _dedup_key_strip_free_ne(key_full, bound)
        dedup in seen_keys && return
        conj_full = _canonical_key(_avg_conj(x), canon)
        dedup_conj = _dedup_key_strip_free_ne(conj_full, bound)
        push!(seen_keys, dedup)
        push!(seen_keys, dedup_conj)
        push!(missing_states, average(key_full))
        if get_adjoints && dedup_conj != dedup
            push!(missing_states, average(conj_full))
        end
        return
    end
    SymbolicUtils.iscall(x) || return
    for a in SymbolicUtils.arguments(x)
        _collect_missing!(missing_states, seen_keys, canon, bound, a, get_adjoints)
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
    # Restrict canon slots to indices that are NOT bound by a sum in H or
    # claimed as a collective-jump index. Slots beyond this filtered set are
    # minted on demand inside `_canonical_key` via `Index(...)(k)`, which
    # preserves the user's naming root while keeping every canonical slot
    # disjoint from H/J bound names. If the bound names appeared as canon
    # slots, `_derive_for` would call `meanfield` on a LHS whose free index
    # re-clashes with H's bound scope, silently dropping the cross-atom
    # commutator term.
    bound = _bound_indices(eqs)
    for (space, v) in canon
        original = copy(v)
        filter!(idx -> !(idx in bound), v)
        # If every user-declared index on this space is bound (e.g. user
        # reused both `i` and `j` for two sum scopes and left no free LHS
        # name on the atom space), mint a single canonical slot as the
        # successor of the lex-first declared name. Without this fallback,
        # `complete!` would discover atom states with whatever name the
        # commutator happened to produce, mixing `i_2`, `j`, etc. across
        # iterations and breaking dedup.
        if isempty(v) && !isempty(original)
            sort!(original, by = idx -> idx.name)
            push!(v, original[1](2))
        end
        sort!(v, by = idx -> idx.name)
    end
    return canon
end

# Indices the Liouvillian treats as bound: either explicitly bound by a sum
# scope inside H, or appearing on a jump operator (collective decay sums over
# the jump's index). When `complete!` derives equations for a state whose
# free index name coincides with one of these, the inner commutator sees a
# spurious clash and drops cross-atom terms. `_derive_for` alpha-renames the
# LHS away from these names before calling `meanfield`.
function _bound_indices(eqs::AbstractMeanFieldEquations)
    out = Set{SQA.Index}()
    _collect_indices_from_qadd_bound!(out, eqs.hamiltonian)
    for j in _flatten_jumps(eqs.jumps)
        _collect_indices_from_qadd_bound!(out, j)
    end
    for j in _flatten_jumps(eqs.jumps_dagger)
        _collect_indices_from_qadd_bound!(out, j)
    end
    return out
end

# Bound from H's perspective: sum-scope `.indices` and free indices on any
# atom inside H (the Liouvillian treats J entries as collective sums over
# their carried index).
_collect_indices_from_qadd_bound!(::Set{SQA.Index}, ::Any) = nothing
function _collect_indices_from_qadd_bound!(out::Set{SQA.Index}, q::SQA.QAdd)
    for idx in q.indices
        push!(out, idx)
    end
    for (term, _) in q.arguments, o in term.ops
        SQA.has_index(o.index) && push!(out, o.index)
    end
    return nothing
end
function _collect_indices_from_qadd_bound!(out::Set{SQA.Index}, q::SQA.QSym)
    SQA.has_index(q.index) && push!(out, q.index)
    return nothing
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

function _canonical_key(
        x::SymbolicUtils.BasicSymbolic, canon::_CanonIndex,
        bound::Set{SQA.Index} = Set{SQA.Index}(),
    )
    op_raw = SQA.undo_average(x)
    op_raw isa SQA.QAdd || return op_raw
    op = _strip_irrelevant_metadata(op_raw)
    encountered = _free_op_indices(op)
    pos_by_space = Dict{Int, Int}()
    rename = Dict{SQA.Index, SQA.Index}()
    for idx in encountered
        pos = get(pos_by_space, idx.space_index, 0) + 1
        pos_by_space[idx.space_index] = pos
        space_canon = get(canon, idx.space_index, nothing)
        (space_canon === nothing || isempty(space_canon)) && continue
        # Slot 1 = first user-declared free index on this space; slot k>1 is
        # minted from the same vocabulary via `Index(...)(k)` (e.g. `j` ↦
        # `j_2`). Minting on demand keeps canon names disjoint from H/J
        # bound names regardless of how many free indices a state has.
        target = pos <= length(space_canon) ? space_canon[pos] :
            space_canon[1](pos)
        target == idx && continue
        rename[idx] = target
    end
    if isempty(rename)
        return SQA.QAdd(op.arguments, SQA.Index[])
    end
    # Batched `change_index` does the rename simultaneously, so a target
    # name that coincides with another encountered index (e.g. swapping
    # `k <-> j` under canon `[j, k]`) does not fuse them destructively
    # mid-rename.
    result = SQA.change_index(op, rename)
    return SQA.QAdd(result.arguments, SQA.Index[])
end

function _strip_irrelevant_metadata(op_raw::SQA.QAdd)
    new_args = SQA.QTermDict()
    for (term, c) in op_raw.arguments
        op_indices = Set{SQA.Index}()
        for o in term.ops
            SQA.has_index(o.index) && push!(op_indices, o.index)
        end
        kept = SQA.NonEqualPair[]
        for pair in term.ne
            (pair[1] in op_indices && pair[2] in op_indices) || continue
            push!(kept, pair)
        end
        new_term = length(kept) == length(term.ne) ? term :
            SQA.QTerm(copy(term.ops), kept)
        new_args[new_term] = c
    end
    return SQA.QAdd(new_args, SQA.Index[])
end

# Build a dedup key for `find_missing` from `_canonical_key`'s output by
# stripping NE pairs whose BOTH indices are free LHS slots (not bound by
# an H/J sum). Two distinct LHS slot indices already designate distinct
# atoms by name; the NE constraint adds no physical identity. Without
# this strip, the same physical state appears twice in `seen_keys` (once
# with NE, once without) and the closure doubles, breaking unique_squeezing
# and other indexed cross-atom systems. Filter-cavity-style NE pairs are
# preserved because at least one index is H-bound (diagonal-vs-offdiagonal
# matters).
function _dedup_key_strip_free_ne(q::SQA.QAdd, bound::Set{SQA.Index})
    out = SQA.QTermDict()
    for (term, c) in q.arguments
        isempty(term.ne) && (out[term] = c; continue)
        kept = SQA.NonEqualPair[]
        for pair in term.ne
            (pair[1] in bound || pair[2] in bound) || continue
            push!(kept, pair)
        end
        new_term = length(kept) == length(term.ne) ? term :
            SQA.QTerm(copy(term.ops), kept)
        out[new_term] = c
    end
    return SQA.QAdd(out, q.indices)
end
_dedup_key_strip_free_ne(x, _) = x

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
    bound = _bound_indices(eqs)
    fresh_ops, undo = _alpha_rename_away(new_ops, bound)
    derived = _meanfield_deterministic(
        eqs.direction, fresh_ops, eqs.hamiltonian, eqs.jumps, eqs.jumps_dagger,
        eqs.rates, eqs.order, mix_choice, eqs.iv,
    )
    return isempty(undo) ? derived : _apply_undo(derived, undo)
end

function _derive_for(
        eqs::NoiseMeanFieldEquations, new_ops; mix_choice = maximum
    )
    bound = _bound_indices(eqs)
    fresh_ops, undo = _alpha_rename_away(new_ops, bound)
    distinct = _distinct_atom_indices(fresh_ops)
    derived = _meanfield_noise(
        eqs.direction, fresh_ops, eqs.hamiltonian, eqs.jumps, eqs.jumps_dagger,
        eqs.rates, eqs.efficiencies, eqs.order, mix_choice, eqs.iv;
        distinct_indices = distinct,
    )
    return isempty(undo) ? derived : _apply_undo(derived, undo)
end

# Free LHS indices on `new_ops` that live on an N-level (atom) subspace,
# deduped. Each is a slot minted by `_canonical_key` to represent a
# physically distinct atom in a multi-atom cross moment; asserting them
# distinct enables SQA's same-site collapse via `expand_completeness`
# (`σ^{gg} = 1 - Σ σ^{kk}`). The distinctness assumption is restricted to
# Transition-carrying indices because Fock-space (filter, mode) indices
# refer to physically distinct sites by user construction and do not need
# (nor benefit from) algebraic same-site collapse; injecting NE on them
# breaks indexed-filter examples like filter-cavity_indexed.
function _distinct_atom_indices(new_ops)
    by_space = Dict{Int, Vector{SQA.Index}}()
    for op in new_ops
        _collect_atom_space_indices_by_space!(by_space, op)
    end
    # Only assert distinctness when this op-set carries 2+ indices on the
    # SAME atom space (the configuration that benefits from same-site
    # collapse). A single atom-space index per derivation (single-atom
    # moment) does not need NE injection, and triggering SQA's
    # canonicalisation on such terms changes drift in examples like
    # unique_squeezing that already produce correct results without it.
    out = SQA.Index[]
    seen = Set{SQA.Index}()
    for (_, idxs) in by_space
        length(idxs) >= 2 || continue
        for idx in idxs
            idx in seen && continue
            push!(seen, idx); push!(out, idx)
        end
    end
    return out
end

_collect_atom_space_indices_by_space!(_, ::Any) = nothing
function _collect_atom_space_indices_by_space!(by_space, op::SQA.QSym)
    op isa SQA.Transition || return nothing
    SQA.has_index(op.index) || return nothing
    v = get!(by_space, op.index.space_index, SQA.Index[])
    op.index in v || push!(v, op.index)
    return nothing
end
function _collect_atom_space_indices_by_space!(by_space, op::SQA.QAdd)
    for (term, _) in op.arguments, o in term.ops
        _collect_atom_space_indices_by_space!(by_space, o)
    end
    return nothing
end

# Pick fresh names for any LHS free index that clashes with an H/J bound name.
# The fresh names are minted from the same user-declared index via the
# `Index(...)(k)` convention so they stay in the user's vocabulary. The undo
# map flips the rename so derived equations can be reported in canonical
# (user-facing) names.
function _alpha_rename_away(ops::AbstractVector, bound::Set{SQA.Index})
    isempty(bound) && return ops, Dict{SQA.Index, SQA.Index}()
    rename = Dict{SQA.Index, SQA.Index}()
    for op in ops
        for idx in _all_indices(op)
            (idx in bound) || continue
            haskey(rename, idx) && continue
            # Mint a successor name from the same root index that is not in
            # `bound` nor already a target. We expand k=2,3,... so the chosen
            # fresh name is deterministic across runs.
            k = 2
            local target::SQA.Index
            while true
                target = idx(k)
                (target in bound) || (target in values(rename)) || break
                k += 1
            end
            rename[idx] = target
        end
    end
    isempty(rename) && return ops, Dict{SQA.Index, SQA.Index}()
    fresh_ops = [_rename_op(op, rename) for op in ops]
    undo = Dict{SQA.Index, SQA.Index}(v => k for (k, v) in rename)
    return fresh_ops, undo
end

_rename_op(op, rename) = op
_rename_op(op::SQA.QSym, rename) = SQA.change_index(op, rename)
function _rename_op(op::SQA.QAdd, rename)
    # If any sum-bound index name collides with a rename target, mint a
    # successor for the bound index first so the undo of `fresh → orig`
    # doesn't fuse a free `orig` with the sum's bound scope.
    isempty(op.indices) && return SQA.change_index(op, rename)
    bound_renames = Dict{SQA.Index, SQA.Index}()
    targets = Set{SQA.Index}(values(rename))
    for bidx in op.indices
        bidx in targets || continue
        k = 2
        local fresh::SQA.Index
        while true
            fresh = bidx(k)
            (fresh in targets) || (fresh in keys(rename)) || (fresh in op.indices) || break
            k += 1
        end
        bound_renames[bidx] = fresh
    end
    if isempty(bound_renames)
        return SQA.change_index(op, rename)
    end
    relabelled = SQA.change_index(op, bound_renames)
    return SQA.change_index(relabelled, rename)
end

_all_indices(op) = SQA.Index[]
_all_indices(op::SQA.QSym) = SQA.has_index(op.index) ? SQA.Index[op.index] : SQA.Index[]
function _all_indices(op::SQA.QAdd)
    out = SQA.Index[]
    for idx in op.indices
        idx in out || push!(out, idx)
    end
    for (term, _) in op.arguments, o in term.ops
        SQA.has_index(o.index) || continue
        o.index in out || push!(out, o.index)
    end
    return out
end

# Walk a derived MeanFieldEquations and substitute back from fresh names to
# the user-visible canonical names so eqs.states / eqs.equations / etc.
# remain in the user's vocabulary.
function _apply_undo(eqs::MeanFieldEquations, undo::Dict{SQA.Index, SQA.Index})
    new_equations = [
        _undo_in_expr(eq.lhs, undo) ~ _undo_in_expr(eq.rhs, undo) for eq in eqs.equations
    ]
    new_op_eqs = [
        _undo_in_op_eq(oe, undo) for oe in eqs.operator_equations
    ]
    new_states = [_undo_in_expr(s, undo) for s in eqs.states]
    new_operators = [_rename_op(o, undo) for o in eqs.operators]
    return MeanFieldEquations(
        new_equations, new_op_eqs, new_states, new_operators,
        eqs.hamiltonian, eqs.jumps, eqs.jumps_dagger, eqs.rates,
        eqs.iv, eqs.order, eqs.direction,
    )
end

function _apply_undo(eqs::NoiseMeanFieldEquations, undo::Dict{SQA.Index, SQA.Index})
    new_equations = [
        _undo_in_expr(eq.lhs, undo) ~ _undo_in_expr(eq.rhs, undo) for eq in eqs.equations
    ]
    new_noise_eqs = [
        _undo_in_expr(eq.lhs, undo) ~ _undo_in_expr(eq.rhs, undo) for eq in eqs.noise_equations
    ]
    new_op_eqs = [
        _undo_in_op_eq(oe, undo) for oe in eqs.operator_equations
    ]
    new_op_noise_eqs = [
        _undo_in_op_eq(oe, undo) for oe in eqs.operator_noise_equations
    ]
    new_states = [_undo_in_expr(s, undo) for s in eqs.states]
    new_operators = [_rename_op(o, undo) for o in eqs.operators]
    return NoiseMeanFieldEquations(
        new_equations, new_noise_eqs, new_op_eqs, new_op_noise_eqs,
        new_states, new_operators,
        eqs.hamiltonian, eqs.jumps, eqs.jumps_dagger, eqs.rates,
        eqs.efficiencies, eqs.iv, eqs.order, eqs.direction,
    )
end

_undo_in_op_eq(oe, undo) = _rename_op(oe.lhs, undo) ~ _rename_op(oe.rhs, undo)

function _undo_in_expr(x, undo::Dict{SQA.Index, SQA.Index})
    isempty(undo) && return x
    return _undo_walk(x, undo)
end

function _undo_walk(x, undo)
    x isa SymbolicUtils.BasicSymbolic || return x
    if _is_leaf_average(x)
        op = SQA.undo_average(x)
        new_op = _rename_op(op, undo)
        new_op === op && return x
        return average(new_op)
    end
    SymbolicUtils.iscall(x) || return x
    f = SymbolicUtils.operation(x)
    args = SymbolicUtils.arguments(x)
    new_args = Any[_undo_walk(a, undo) for a in args]
    if all(==(true), (a === b for (a, b) in zip(args, new_args)))
        return x
    end
    return f(new_args...)
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
