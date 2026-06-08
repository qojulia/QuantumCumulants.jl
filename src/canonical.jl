"""
Canonicalisation context shared by the moment-keying functions. Holds the
per-subspace vocabulary of canonical index names and the policy flags that
`canon_key`, `orbit_key` and `canonical_rep` read.
"""
struct CanonCtx
    vocab::Dict{Int, Vector{SQA.Index}}   # per-subspace canonical free-index reps, disjoint from H/J-bound names
    symmetric::Set{Int}                   # subspaces carrying an Index (permutation-symmetric families)
    selected::Set{Int}                    # subspaces in scope for scale/evaluate; empty means all
    population::Bool                      # closure policy: dephasing/population system vs concrete-site
end

"""
Per-subspace representation regime of a moment under canonicalisation:
- `Free`: a symbolic index, with no permutation or concrete quotient applied.
- `Scaled`: quotiented by the permutation symmetry ``S_n`` of the atoms.
- `Concrete`: instantiated to a fixed site `1..M`, per-site, with no alpha-renaming.
"""
@enum Coordinate Free Scaled Concrete

"""
Coordinate map placing every symmetric subspace of `ctx` in the `Free` regime.
Scalar systems (no symmetric subspace) yield an empty map, which the keying
functions read as all-`Free`.
"""
all_free_coords(ctx::CanonCtx) = Dict{Int, Coordinate}(sp => Free for sp in ctx.symmetric)

"""
Return a copy of `coords` with subspace `sp` set to coordinate `c`.
"""
function with_coord(coords::Dict{Int, Coordinate}, sp::Int, c::Coordinate)
    out = copy(coords)
    out[sp] = c
    return out
end

"""
Coordinate of subspace `sp` in `coords`, defaulting to `Free` when unlisted.
"""
coord_of(coords::Dict{Int, Coordinate}, sp::Int) = get(coords, sp, Free)

"""
Construct the canonicalisation context for a system. The canonical index
vocabulary is derived from the indices the user declared on the operators,
Hamiltonian and jumps, excluding any name bound by a Hamiltonian sum or carried by
a jump (those must stay free). When a subspace needs more distinct representatives
than the user declared, extra ones `rep(k)` are minted from its first declared name.
"""
function build_ctx(
        ops::AbstractVector, H::QField, J::AbstractVector, Jdagger::AbstractVector;
        selected::AbstractVector{Int} = Int[],
    )
    flat_j = _flatten_jumps(J)
    flat_jd = _flatten_jumps(Jdagger)
    sources = Iterators.flatten(((H,), ops, flat_j, flat_jd))

    vocab = Dict{Int, Vector{SQA.Index}}()
    symmetric = Set{Int}()
    for src in sources, idx in SQA.get_indices(src)
        push!(symmetric, idx.space_index)
        v = get!(vocab, idx.space_index, SQA.Index[])
        idx in v || push!(v, idx)
    end

    bound = Set{SQA.Index}()
    _collect_indices_from_qadd_bound!(bound, H)
    for src in Iterators.flatten((flat_j, flat_jd))
        _collect_indices_from_qadd_bound!(bound, src)
    end
    for (space, v) in vocab
        original = copy(v)
        filter!(idx -> !(idx in bound), v)
        if isempty(v) && !isempty(original)
            sort!(original, by = idx -> idx.name)
            push!(v, original[1](2))
        end
        sort!(v, by = idx -> idx.name)
    end

    population = _has_dephasing_channel(_flatten_jumps(J))
    return CanonCtx(vocab, symmetric, Set{Int}(selected), population)
end

"""
Alpha-canonical identity key of `op`, used as the node key during completion. Drops
the sum scope, fixes the order of commuting operators, alpha-renames free indices to
the context vocabulary and drops all non-equal constraints. Conjugate pairs are *not*
folded; use `canonical_rep` for that. Equivalent to `canonical_rep` with every
subspace `Free`; constants pass through unchanged.
"""
canon_key(op::QAdd, ctx::CanonCtx) = _coord_key(op, ctx, all_free_coords(ctx))
canon_key(op, ::CanonCtx) = op   # non-QAdd (constant) passthrough

"""
Conjugation-folded identity key of `op` under the per-subspace `coords`. Each
subspace is normalised according to its coordinate: `Free` alpha-renames to the
vocabulary representatives, `Scaled` additionally folds under permutation symmetry
via `symmetric_min`, and `Concrete` keeps its literal index names. The operator and
its adjoint are compared and the lexicographically smaller is returned as `key`,
with `is_conjugate` reporting whether `op` itself is the conjugate side.
"""
function canonical_rep(op::QAdd, ctx::CanonCtx; coords::Dict{Int, Coordinate} = all_free_coords(ctx))
    k = _coord_key(op, ctx, coords)
    kc = _coord_key(adjoint(op), ctx, coords)
    return _serialize(kc) < _serialize(k) ? (kc, true) : (k, false)
end
canonical_rep(op, ::CanonCtx; kw...) = (op, false)

"""
Non-conjugate-folded identity key of `op`. Brings `op` to normal form (sum scope
dropped, commuting operators ordered, NE constraints removed), then applies each
subspace's coordinate quotient from `coords`. Shared building block of `canon_key`,
`orbit_key` and `canonical_rep`.
"""
function _coord_key(op::QAdd, ctx::CanonCtx, coords::Dict{Int, Coordinate})
    scaled = Set{Int}(sp for sp in keys(coords) if coords[sp] == Scaled)
    concrete = Set{Int}(sp for sp in keys(coords) if coords[sp] == Concrete)
    # Alpha-rename everything not Concrete (Free and Scaled relabel to vocab reps).
    rename_spaces = Set{Int}()
    for sp in ctx.symmetric
        sp in concrete || push!(rename_spaces, sp)
    end

    base = _reorder_commuting(_drop_scope_ne(SQA.QAdd(op.arguments, SQA.Index[])))   # totality on all subspaces
    base = _alpha_rename_spaces(base, ctx, rename_spaces)
    key = _drop_all_ne(SQA.QAdd(base.arguments, SQA.Index[]))
    isempty(scaled) && return key
    return symmetric_min(key, ctx, scaled)
end
_coord_key(op, ::CanonCtx, ::Dict{Int, Coordinate}) = op

"""
Alpha-rename the free indices of `op` to the context vocabulary, but only on the
subspaces listed in `spaces`. Indices on other subspaces (notably `Concrete` ones)
keep their names.
"""
function _alpha_rename_spaces(op::QAdd, ctx::CanonCtx, spaces::Set{Int})
    encountered = _free_op_indices(op)
    rename = Dict{SQA.Index, SQA.Index}()
    pos_by_space = Dict{Int, Int}()
    for idx in encountered
        idx.space_index in spaces || continue
        pos = get(pos_by_space, idx.space_index, 0) + 1
        pos_by_space[idx.space_index] = pos
        reps = get(ctx.vocab, idx.space_index, SQA.Index[])
        isempty(reps) && continue
        target = pos <= length(reps) ? reps[pos] : reps[1](pos)
        target == idx && continue
        rename[idx] = target
    end
    return isempty(rename) ? op : SQA.change_index(op, rename)
end

"""
Identity key for already-materialised systems. Like `canon_key` but *without*
alpha-renaming, so the concrete per-site atoms (`σ_1`, `σ_2`, …) produced by
`evaluate` stay distinct instead of collapsing onto one vocabulary slot. Still drops
the sum scope, fixes commuting-op order and drops NE constraints; constants pass
through unchanged.
"""
function literal_key(op::QAdd)
    base = _reorder_commuting(SQA.QAdd(op.arguments, SQA.Index[]))
    return _drop_all_ne(SQA.QAdd(base.arguments, SQA.Index[]))
end
literal_key(op) = op

"""
Canonicalise the order of commuting same-subspace factors in `op` by asserting their
free indices mutually distinct, which lets SQA sort them into a fixed order. The
asserted non-equal constraints are temporary and are stripped from the final key by
`_drop_all_ne`.
"""
function _reorder_commuting(op::QAdd)
    by_space = Dict{Int, Vector{SQA.Index}}()
    for (term, _) in op.arguments, o in term.ops
        SQA.has_index(o.index) || continue
        v = get!(by_space, o.index.space_index, SQA.Index[])
        o.index in v || push!(v, o.index)
    end
    pairs = Tuple{SQA.Index, SQA.Index}[]
    for (_, idxs) in by_space, n in 1:length(idxs), m in (n + 1):length(idxs)
        push!(pairs, (idxs[n], idxs[m]))
    end
    return isempty(pairs) ? op : SQA.assume_distinct_index(op, pairs)
end

"""
Substitution mapping each free operator index of `op`, in first-encounter order per
subspace, to its canonical vocabulary representative (minting `rep(k)` beyond the
declared list).
"""
function _alpha_rename(op::QAdd, ctx::CanonCtx)
    encountered = _free_op_indices(op)
    rename = Dict{SQA.Index, SQA.Index}()
    pos_by_space = Dict{Int, Int}()
    for idx in encountered
        pos = get(pos_by_space, idx.space_index, 0) + 1
        pos_by_space[idx.space_index] = pos
        reps = get(ctx.vocab, idx.space_index, SQA.Index[])
        isempty(reps) && continue
        target = pos <= length(reps) ? reps[pos] : reps[1](pos)
        target == idx && continue
        rename[idx] = target
    end
    return rename
end

"""
Drop non-equal pairs whose index is carried by no operator in the term; these are
sum-scope constraints (e.g. `Σ_{j≠i₂}`), not part of the moment's identity. Keeping
them would let a later alpha-rename collapse the index onto its partner and
annihilate the term. Operator-internal NE pairs are preserved.
"""
function _drop_scope_ne(q::QAdd)
    out = SQA.QTermDict()
    changed = false
    for (term, c) in q.arguments
        opidx = Set{SQA.Index}()
        for o in term.ops
            SQA.has_index(o.index) && push!(opidx, o.index)
        end
        kept = SQA.NonEqualPair[p for p in term.ne if p[1] in opidx && p[2] in opidx]
        length(kept) == length(term.ne) || (changed = true)
        nt = SQA.QTerm(copy(term.ops), kept)
        out[nt] = get(out, nt, zero(c)) + c
    end
    return changed ? SQA.QAdd(out, q.indices) : q
end

"""
Strip every non-equal constraint from `q`, merging any terms that become identical
once the constraints are gone.
"""
function _drop_all_ne(q::QAdd)
    out = SQA.QTermDict()
    for (term, c) in q.arguments
        new_term = isempty(term.ne) ? term : SQA.QTerm(copy(term.ops), SQA.NonEqualPair[])
        # Accumulate: terms collapsing to the same ops must add their coefficients.
        out[new_term] = get(out, new_term, zero(c)) + c
    end
    return SQA.QAdd(out, q.indices)
end

"""
Permutation-symmetry identity key of `op`, used as the node key during `scale`. It is
`canon_key` plus a fold over the permutation symmetry of the `selected` subspaces
(the rest stay `Free`), and it folds conjugate pairs. Constants pass through unchanged.
"""
function orbit_key(op::QAdd, ctx::CanonCtx; selected = ctx.symmetric)
    coords = all_free_coords(ctx)
    for sp in selected
        coords[sp] = Scaled
    end
    return _coord_key(op, ctx, coords)
end

"""
Identity key appropriate to a system's current `coords`. Once any subspace is
`Concrete` (a prior `evaluate` materialised it), the conjugation-folded
`canonical_rep` is used so the closed set matches what `evaluate` and codegen
produce; while every subspace is still symbolic, the unfolded `_coord_key` is used.
"""
function _materialised_key(op::QAdd, ctx::CanonCtx, coords::Dict{Int, Coordinate})
    any(==(Concrete), values(coords)) && return canonical_rep(op, ctx; coords)[1]
    return _coord_key(op, ctx, coords)
end
_materialised_key(op, ::CanonCtx, ::Dict{Int, Coordinate}) = op
orbit_key(op, ::CanonCtx; kw...) = op

"""
Representative of `op` under the permutation symmetry of the `selected` subspaces.
Brute-forces every relabelling of the (few) atom slots and the commuting-op order
each induces, returning the lexicographically smallest rendered form. The number of
slots is small (2 to 4) and this runs at build time.
"""
function symmetric_min(op::QAdd, ctx::CanonCtx, selected)
    free = _free_by_space(op)
    spaces = sort!([sp for sp in keys(free) if sp in selected])
    isempty(spaces) && return op
    perms_per_space = [collect(Combinatorics.permutations(free[sp])) for sp in spaces]
    best_op = nothing
    best_key = nothing
    for combo in Iterators.product(perms_per_space...)
        sub = Dict{SQA.Index, SQA.Index}()
        for (n, sp) in enumerate(spaces)
            reps = get(ctx.vocab, sp, SQA.Index[])
            for (k, idx) in enumerate(combo[n])
                # Mint concrete slots (as `evaluate` does) so scale and evaluate agree.
                target = isempty(reps) ? idx : reps[1](k)
                target == idx && continue
                sub[idx] = target
            end
        end
        cand = isempty(sub) ? op : SQA.change_index(op, sub)
        cand = _reorder_then_drop_ne(cand, spaces)
        key = _serialize(cand)
        if best_key === nothing || key < best_key
            best_key = key
            best_op = cand
        end
    end
    return best_op === nothing ? op : best_op
end

"""
Free operator indices of `op` grouped by subspace, in first-encounter order.
"""
function _free_by_space(op::QAdd)
    out = Dict{Int, Vector{SQA.Index}}()
    for (term, _) in op.arguments, o in term.ops
        SQA.has_index(o.index) || continue
        v = get!(out, o.index.space_index, SQA.Index[])
        o.index in v || push!(v, o.index)
    end
    return out
end

"""
Canonicalise the order of commuting cross-atom factors on `spaces` (by asserting the
slots distinct so SQA sorts them), then strip the asserted non-equal constraints from
the result.
"""
function _reorder_then_drop_ne(op::QAdd, spaces)
    by_space = Dict{Int, Vector{SQA.Index}}()
    for (term, _) in op.arguments, o in term.ops
        SQA.has_index(o.index) || continue
        o.index.space_index in spaces || continue
        v = get!(by_space, o.index.space_index, SQA.Index[])
        o.index in v || push!(v, o.index)
    end
    pairs = Tuple{SQA.Index, SQA.Index}[]
    for (_, idxs) in by_space, n in 1:length(idxs), m in (n + 1):length(idxs)
        push!(pairs, (idxs[n], idxs[m]))
    end
    reordered = isempty(pairs) ? op : SQA.assume_distinct_index(op, pairs)
    return _drop_all_ne(SQA.QAdd(reordered.arguments, SQA.Index[]))
end

"""
Stable, order-independent string signature of `op`, used to compare candidate
representatives lexicographically.
"""
function _serialize(op::QAdd)
    terms = Vector{Tuple{Vector{String}, String}}()
    for (term, c) in op.arguments
        push!(terms, (String[string(o) for o in term.ops], string(c)))
    end
    sort!(terms)
    return terms
end

"""
Conjugation orbit representative for already-materialised systems, the `literal_key`
analogue of `canonical_rep`. Returns the lexicographically smaller of `literal_key(op)`
and `literal_key(op')` as `key`, with `is_conjugate` reporting whether `op` itself is
the conjugate side. Constants pass through unchanged.
"""
function concrete_rep(op::QAdd)
    k = literal_key(op)
    kc = literal_key(adjoint(op))
    return _serialize(kc) < _serialize(k) ? (kc, true) : (k, false)
end
concrete_rep(op) = (op, false)   # non-QAdd (constant) passthrough
