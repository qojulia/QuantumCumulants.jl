# Moment identity (spec Section 3.1). Two levels: `canon_key` (alpha-canonical,
# the complete/closure node key) and `orbit_key` (canon_key then symmetric_min,
# the scale node key). Identity carries NO NE (spec Section 6).

struct CanonCtx
    vocab::Dict{Int, Vector{SQA.Index}}   # per-subspace canonical free-index reps, disjoint from H/J-bound names
    symmetric::Set{Int}                   # subspaces carrying an Index (permutation-symmetric families)
    selected::Set{Int}                    # subspaces in scope for scale/evaluate; empty means all
    population::Bool                       # closure policy: dephasing/population system vs concrete-site
end

# Per-subspace coordinate state of a moment representation (spec "Coordinates").
# Free: symbolic index, no permutation/concrete quotient. Scaled: permutation
# (Sₙ) quotient. Concrete: instantiated to 1..M, per-site, no alpha-rename.
@enum Coordinate Free Scaled Concrete

# Default coordinate map: every symmetric subspace Free. The scalar systems (no
# indexed subspace) get an empty map, which `canonical_rep` reads as all-Free.
all_free_coords(ctx::CanonCtx) = Dict{Int, Coordinate}(sp => Free for sp in ctx.symmetric)

# Functional update: set subspace `sp` to coordinate `c`.
function with_coord(coords::Dict{Int, Coordinate}, sp::Int, c::Coordinate)
    out = copy(coords)
    out[sp] = c
    return out
end

# Read a subspace's coordinate, defaulting to Free for subspaces not listed.
coord_of(coords::Dict{Int, Coordinate}, sp::Int) = get(coords, sp, Free)

# Single replacement for the scattered `_build_canonical_indices` calls. Derives
# `vocab` from the user's declared indices on ops/H/J, filtered to exclude any
# index bound by an H sum or carried by a jump (those names must stay free for
# `derive` to not re-clash). Mints `rep1(k)` slots when a subspace needs more
# reps than the user declared.
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

# canon_key = drop sum scope, canonical commuting-op order, alpha-canonical
# free-index rename, drop ALL NE. The commuting-op order step (totality) makes
# the SAME physical moment built in different operator orders map to one key; it
# does NOT fold conjugate pairs (that is orbit_key's symmetric_min, which also
# minimizes over slot relabelings). See spec Section 3.1.
# canon_key = canonical_rep configured with every subspace Free.
canon_key(op::QAdd, ctx::CanonCtx) = _coord_key(op, ctx, all_free_coords(ctx))
canon_key(op, ::CanonCtx) = op   # non-QAdd (constant) passthrough

# The single source of truth (spec "The single source of truth"). Returns
# `(representative_key, is_conjugate)`. Per subspace:
#   Free     -> alpha-rename to vocab reps (as canon_key does)
#   Scaled   -> alpha-rename then permutation-fold via symmetric_min
#   Concrete -> keep concrete names (as literal_key does), no alpha-rename
# Conjugation is ALWAYS part of identity: fold {O, O†} to the min over the total
# `_serialize` order, returning the sign bit. The four legacy keys are
# configurations of this one function (see the wrappers below).
function canonical_rep(op::QAdd, ctx::CanonCtx; coords::Dict{Int, Coordinate} = all_free_coords(ctx))
    k = _coord_key(op, ctx, coords)
    kc = _coord_key(adjoint(op), ctx, coords)
    return _serialize(kc) < _serialize(k) ? (kc, true) : (k, false)
end
canonical_rep(op, ::CanonCtx; kw...) = (op, false)

# The non-conjugate-folded key for a coordinate map: normal-form, then per
# subspace apply the coordinate's quotient. Built from the existing pieces so
# the behavior matches the legacy keys where coordinates coincide.
function _coord_key(op::QAdd, ctx::CanonCtx, coords::Dict{Int, Coordinate})
    scaled = Set{Int}(sp for sp in keys(coords) if coords[sp] == Scaled)
    concrete = Set{Int}(sp for sp in keys(coords) if coords[sp] == Concrete)
    # Subspaces absent from `coords` default to Free. The set of subspaces to
    # alpha-rename is everything NOT Concrete (Free and Scaled both relabel to
    # vocab reps; Scaled then folds via symmetric_min, which re-mints slots).
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

# Like `_alpha_rename`, but only renames indices whose subspace is in `spaces`.
# Concrete subspaces are skipped so their per-site names survive (literal_key).
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

# literal_key = canon_key WITHOUT the alpha-rename. Drops sum scope, fixes
# commuting-op order, drops NE, but KEEPS concrete index names. This is the
# post-`specialize` identity: after evaluate, indices are concrete per-site
# atoms (`σ_1`, `σ_2`), which must stay distinct, so alpha-renaming them to one
# vocab slot (as `canon_key` does) would wrongly merge different atoms. See the
# CLAUDE.md note "after evaluate the dedup key is the operator itself".
function literal_key(op::QAdd)
    base = _reorder_commuting(SQA.QAdd(op.arguments, SQA.Index[]))
    return _drop_all_ne(SQA.QAdd(base.arguments, SQA.Index[]))
end
literal_key(op) = op

# Assert same-space free indices pairwise distinct so SQA canonicalizes the order
# of the commuting operators they sit on, BEFORE alpha-rename (which assigns slots
# by encounter order, so the op order must be fixed first). A multi-free-index
# moment ⟨σ_i σ_j⟩ means distinct atoms i ≠ j, so the assertion is physically
# correct; the NE it adds is dropped from the final key by `_drop_all_ne`. This is
# the totality step (fix order); it does NOT minimize over relabelings, so it does
# not fold conjugates (unlike symmetric_min).
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

# Map each free operator index, in first-encounter order per subspace, to its
# canonical vocab rep (minting `rep1(k)` past the declared list).
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

# NE pairs referencing an index carried by no operator are sum-scope constraints
# (e.g. `Σ_{j≠i₂}`: j bound, i₂ external/LHS), not part of the moment's identity.
# Drop them before alpha-rename, else the bound index renames onto the external
# partner's vocab slot and `change_index` annihilates the term (i₂≠i₂ → 0). Op-
# internal NE (genuine multi-atom ⟨σ_aσ_b⟩, a≠b) is kept; `_reorder_commuting` /
# `_drop_all_ne` handle it. This extends "scope is transparent to identity" from
# `.indices` (already stripped) to `.ne`.
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

function _drop_all_ne(q::QAdd)
    out = SQA.QTermDict()
    for (term, c) in q.arguments
        new_term = isempty(term.ne) ? term : SQA.QTerm(copy(term.ops), SQA.NonEqualPair[])
        # Accumulate (not overwrite): if two terms collapse to the same ops after
        # NE is dropped, their coefficients must add. (Node keys are single-term,
        # so this is hardening, not a currently-triggered path.)
        out[new_term] = get(out, new_term, zero(c)) + c
    end
    return SQA.QAdd(out, q.indices)
end

# orbit_key = canon_key then the permutation-symmetry quotient. The ONLY pass
# that reorders commuting cross-atom ops, so it deliberately folds conjugate
# pairs (spec Section 3.1, Correction 2).
# orbit_key = canonical_rep with the `selected` subspaces Scaled, rest Free.
function orbit_key(op::QAdd, ctx::CanonCtx; selected = ctx.symmetric)
    coords = all_free_coords(ctx)
    for sp in selected
        coords[sp] = Scaled
    end
    return _coord_key(op, ctx, coords)
end

# Node key for a graph in coordinate `coords`. Once any subspace is Concrete the
# system is materialised (a prior `evaluate` ran), and the identity is the
# conjugate-folded `canonical_rep` (the same key the codegen resolver and
# `specialize`'s `concrete_rep` dedup use), so `evaluate`-then-`scale` lands the
# same closed set as `scale`-then-`evaluate`. While every subspace is still
# symbolic (Free/Scaled) the unfolded `_coord_key` is used, preserving the
# `get_adjoints=true` convention (both conjugate partners tracked) that the
# symbolic complete/scale counts lock.
function _materialised_key(op::QAdd, ctx::CanonCtx, coords::Dict{Int, Coordinate})
    any(==(Concrete), values(coords)) && return canonical_rep(op, ctx; coords)[1]
    return _coord_key(op, ctx, coords)
end
_materialised_key(op, ::CanonCtx, ::Dict{Int, Coordinate}) = op
orbit_key(op, ::CanonCtx; kw...) = op

# Brute force over (relabelings of the K slots per selected symmetric subspace)
# times (the canonical commuting-op order each relabel induces). Returns the
# lexicographically-minimal rendered form. K is tiny (2 to 4) and build-time.
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
                # Mint concrete slots `rep(1), rep(2), …` (same as `evaluate`'s
                # `_concrete_index`), so scale's representative for slot 1 is the
                # same `i_1` that evaluate materialises, not the bare rep `i`.
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

# Free indices grouped per subspace, in first-encounter order.
function _free_by_space(op::QAdd)
    out = Dict{Int, Vector{SQA.Index}}()
    for (term, _) in op.arguments, o in term.ops
        SQA.has_index(o.index) || continue
        v = get!(out, o.index.space_index, SQA.Index[])
        o.index in v || push!(v, o.index)
    end
    return out
end

# Assert same-space slots mutually distinct so SQA._canonicalize! orders the
# commuting cross-atom factors, then drop the asserted NE off the key again.
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

# Stable, order-independent string signature for lexicographic comparison.
function _serialize(op::QAdd)
    terms = Vector{Tuple{Vector{String}, String}}()
    for (term, c) in op.arguments
        push!(terms, (String[string(o) for o in term.ops], string(c)))
    end
    sort!(terms)
    return terms
end

# Concrete-layer conjugation orbit representative: the `literal_key` analogue of
# `canonical_rep`. After `specialize` the indices are concrete (or scaled
# orbit-reps), so the inner key is `literal_key` (KEEPS names, no alpha-rename),
# not `canon_key`. Returns `(orbit_min_key, is_conjugate)`: the
# lexicographically-minimal of `{literal_key(op), literal_key(adjoint(op))}` and
# whether `op` sits on the conjugate side of that orbit.
#
# This is the single conjugation-aware identity for the concrete layer (spec
# 2026-06-04). It is adjoint-equivariant BY CONSTRUCTION (a `min` over the total
# `_serialize` order, which commutes with adjoint's order-reversal), so it
# round-trips for the mixed scaled-orbit-rep + concrete-mode states where
# `literal_key` alone does not. `specialize`'s dedup, the codegen leaf resolver,
# and `get_solution` all route conjugate folding through this, so no consumer
# assumes `literal_key` round-trips under `adjoint`.
function concrete_rep(op::QAdd)
    k = literal_key(op)
    kc = literal_key(adjoint(op))
    return _serialize(kc) < _serialize(k) ? (kc, true) : (k, false)
end
concrete_rep(op) = (op, false)   # non-QAdd (constant) passthrough
