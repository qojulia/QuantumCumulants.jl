"""
A treatment fingerprint: the per-subspace `(space_index, treatment)` pairs in sorted
order. It identifies a `treatments` map cheaply for use as a memo key. The map has one
entry per symmetric subspace (typically one to three), so the fingerprint is small to
build and hash.
"""
const TreatmentFP = Vector{Tuple{Int, SubspaceTreatment}}
const _EMPTY_FP = Tuple{Int, SubspaceTreatment}[]
treatment_fp(treatments::Dict{Int, SubspaceTreatment}) =
    isempty(treatments) ? _EMPTY_FP :
    sort!(Tuple{Int, SubspaceTreatment}[(sp, t) for (sp, t) in treatments])

"""
Per-`ctx` identity memo for the canonicalisation engine. `_treatment_key` and
`canonical_rep` are pure in `(op, ctx, treatments)`; with `ctx` fixed they are cached
under the operator and the treatment fingerprint, so a label computed during closure is
reused at `System`-build and spectrum time. The fingerprint keeps the all-`Free`,
`Scaled`, and `Concrete` configurations of one shared `ctx` in distinct slots. The memo
is monotonic and pure: it only ever grows and never needs invalidation.
"""
struct CanonCache
    key::Dict{Tuple{QAdd, TreatmentFP}, QAdd}              # _treatment_key result
    rep::Dict{Tuple{QAdd, TreatmentFP}, Tuple{QAdd, Bool}} # canonical_rep -> (rep, is_conjugate)
    CanonCache() = new(
        Dict{Tuple{QAdd, TreatmentFP}, QAdd}(),
        Dict{Tuple{QAdd, TreatmentFP}, Tuple{QAdd, Bool}}(),
    )
end

"""
Shared state for the moment-key functions (`canon_key`, `scaled_key`, `canonical_rep`,
…). A *key* assigns each operator average ⟨A⟩ a canonical label, identical for two
averages exactly when they are the same expectation value, so the equations of motion
close on a non-redundant set. Two averages can coincide through relabelling of free
atom indices, the permutation symmetry of identical atoms, or Hermitian conjugation
(⟨A†⟩ = ⟨A⟩*); the functions differ in which of these they identify. `CanonCtx` holds
the index vocabulary used for relabelling, which subspaces are symmetric/selected, and a
per-context identity memo (`cache`). Every `Int` key here and in the `treatments` maps is
an SQA `space_index`: the position of a Hilbert-space factor (subspace) in the system's
`ProductSpace`, as returned by `acts_on`.
"""
struct CanonCtx
    # canonical free-index reps per subspace (space_index), disjoint from H/J-bound names
    vocab::Dict{Int, Vector{SQA.Index}}
    # subspaces (space_index) carrying an Index: permutation-symmetric atom families
    symmetric::Set{Int}
    # subspaces (space_index) in scope for scale/evaluate; empty means all
    selected::Set{Int}
    # per-context identity memo, shared across a transform lineage; excluded from identity
    cache::CanonCache
    CanonCtx(vocab, symmetric, selected) = new(vocab, symmetric, selected, CanonCache())
end

"""
Per-subspace treatment map marking every symmetric subspace of `ctx` as `Free`. Scalar
systems (no symmetric subspace) yield an empty map, which the keying functions read
as all-`Free`.
"""
all_free_treatments(ctx::CanonCtx) = Dict{Int, SubspaceTreatment}(sp => Free for sp in ctx.symmetric)

"""
Return a copy of `treatments` with subspace `sp` set to treatment `c`.
"""
function with_treatment(treatments::Dict{Int, SubspaceTreatment}, sp::Int, c::SubspaceTreatment)
    out = copy(treatments)
    out[sp] = c
    return out
end

"""
Treatment of subspace `sp` in `treatments`, defaulting to `Free` when unlisted.
"""
treatment_of(treatments::Dict{Int, SubspaceTreatment}, sp::Int) = get(treatments, sp, Free)

"""
The per-subspace treatment map recorded on `eqs`, falling back to all-`Free` (over the
symmetric subspaces of `ctx`) when `eqs` carries none. Single source for the consumers
that re-key a stored system (`find_missing`, the MTK `System` build, `evaluate`, correlation/spectrum).
"""
_treatments(eqs::AbstractMeanfieldEquations, ctx::CanonCtx) =
    isempty(eqs.graph.treatments) ? all_free_treatments(ctx) : eqs.graph.treatments

"""
Collect the indices the Liouvillian treats as bound: sum-scope `.indices` and any
free index carried by an atom inside the source, since the Liouvillian sums
collective jumps over their carried index. `build_ctx` excludes these names from
the canonical free-index vocabulary so `derive` does not re-clash.
"""
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

"""
Distinct free operator indices of a `QAdd`, in first-encounter order. Used by the
relabelling helpers below and by `evaluate`'s index unrolling.
"""
function _free_op_indices(op::SQA.QAdd)
    out = SQA.Index[]
    for (term, _) in op.arguments, o in term.ops
        SQA.has_index(o.index) || continue
        o.index in out || push!(out, o.index)
    end
    return out
end

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

    return CanonCtx(vocab, symmetric, Set{Int}(selected))
end

"""
The stored canonicalisation context of an equation set: its graph already carries the
`ctx` built at derivation, so every consumer reads the one source instead of recomputing.
"""
build_ctx(eqs::AbstractMeanfieldEquations) = eqs.graph.ctx

"""
Canonical label of the average ⟨`op`⟩ as it enters the cumulant hierarchy: two averages
share it exactly when they are equal after relabelling free atom indices to the context
vocabulary. Hermitian conjugates are *not* identified (use `canonical_rep` for that).
This is `_treatment_key` with every subspace left as a free symbolic index.
"""
canon_key(op::QAdd, ctx::CanonCtx) = _treatment_key(op, ctx, all_free_treatments(ctx))
canon_key(op, ::CanonCtx) = op   # non-QAdd (constant) passthrough

"""
Canonical label that also accounts for Hermitian conjugation. Returns
`(key, is_conjugate)`: `key` is the smaller of the labels of ⟨`op`⟩ and ⟨`op`†⟩, and
`is_conjugate` flags whether `op` was the conjugate side, so an average and its
conjugate ⟨A†⟩ = ⟨A⟩* collapse to one dynamical variable. This is the label used when
generating the numerical code and when deduplicating the hierarchy. `treatments` chooses
the per-subspace treatment exactly as in `_treatment_key`.
"""
function canonical_rep(op::QAdd, ctx::CanonCtx; treatments::Dict{Int, SubspaceTreatment} = all_free_treatments(ctx))
    return get!(ctx.cache.rep, (op, treatment_fp(treatments))) do
        _conjugation_rep(op, o -> _treatment_key(o, ctx, treatments))
    end
end
canonical_rep(op, ::CanonCtx; kw...) = (op, false)

"""
The engine behind `canon_key`, `scaled_key` and `canonical_rep`. Puts ⟨`op`⟩ into a
comparable normal form (drops the summation scope, fixes the order of commuting
factors, removes the non-equal index constraints), then treats each atomic
subspace according to `treatments`: a `Free` or `Scaled` subspace has its indices
relabelled to the canonical reps, `Scaled` additionally reduces under the permutation
symmetry of the identical atoms (via `symmetric_min`), and `Concrete` keeps its fixed
site labels. Hermitian conjugation is not applied.
"""
function _treatment_key(op::QAdd, ctx::CanonCtx, treatments::Dict{Int, SubspaceTreatment})
    return get!(ctx.cache.key, (op, treatment_fp(treatments))) do
        scaled = Set{Int}()
        concrete = Set{Int}()
        for (sp, t) in treatments
            t == Scaled && push!(scaled, sp)
            t == Concrete && push!(concrete, sp)
        end
        # Relabel everything not Concrete (Free and Scaled relabel to vocab reps).
        rename_spaces = Set{Int}()
        for sp in ctx.symmetric
            sp in concrete || push!(rename_spaces, sp)
        end

        base = _reorder_commuting(_drop_scope_non_equal(SQA.QAdd(op.arguments, SQA.Index[])))   # totality on all subspaces
        base = _relabel_spaces(base, ctx, rename_spaces)
        key = _drop_all_non_equal(SQA.QAdd(base.arguments, SQA.Index[]))
        isempty(scaled) && return key
        return symmetric_min(key, ctx, scaled)
    end
end
_treatment_key(op, ::CanonCtx, ::Dict{Int, SubspaceTreatment}) = op

"""
The `k`-th canonical index slot of a subspace, minted from its first declared vocabulary
representative (`reps[1](k)`, per SQA's naming policy). This is the single definition of a
minted slot: `symmetric_min` (scale) and `_concrete_index` (evaluate) both mint through it
so a scaled moment and its materialised form name the same atoms, and `_relabel_spaces`
uses it for slots beyond the declared vocabulary.
"""
nth_index(reps::Vector{SQA.Index}, k::Int) = reps[1](k)

"""
Relabel the free indices of `op` to the context vocabulary, but only on the subspaces
listed in `spaces`. Indices on other subspaces (notably `Concrete` ones) keep their
names.
"""
function _relabel_spaces(op::QAdd, ctx::CanonCtx, spaces::Set{Int})
    encountered = _free_op_indices(op)
    rename = Dict{SQA.Index, SQA.Index}()
    pos_by_space = Dict{Int, Int}()
    for idx in encountered
        idx.space_index in spaces || continue
        pos = get(pos_by_space, idx.space_index, 0) + 1
        pos_by_space[idx.space_index] = pos
        reps = get(ctx.vocab, idx.space_index, SQA.Index[])
        isempty(reps) && continue
        target = pos <= length(reps) ? reps[pos] : nth_index(reps, pos)
        target == idx && continue
        rename[idx] = target
    end
    return isempty(rename) ? op : SQA.change_index(op, rename)
end

"""
Like `canon_key` but without relabelling indices, so the fixed per-site operators
(`σ_1`, `σ_2`, …) created by `evaluate` stay distinct instead of collapsing onto a
single representative atom. The label used once a system has been materialised to
explicit sites.
"""
function concrete_key(op::QAdd)
    base = _reorder_commuting(SQA.QAdd(op.arguments, SQA.Index[]))
    return _drop_all_non_equal(SQA.QAdd(base.arguments, SQA.Index[]))
end
concrete_key(op) = op

"""
Canonicalise the order of commuting same-subspace factors in `op` by asserting their
free indices mutually distinct, which lets SQA sort them into a fixed order. The
asserted non-equal constraints are temporary and are stripped from the final key by
`_drop_all_non_equal`.
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
Drop non-equal pairs whose index is carried by no operator in the term; these are
sum-scope constraints (e.g. `Σ_{j≠i₂}`), not part of the moment's identity. Keeping
them would let a later relabelling collapse the index onto its partner and
annihilate the term. Operator-internal non-equal index pairs are preserved.
"""
function _drop_scope_non_equal(q::QAdd)
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
function _drop_all_non_equal(q::QAdd)
    out = SQA.QTermDict()
    for (term, c) in q.arguments
        new_term = isempty(term.ne) ? term : SQA.QTerm(copy(term.ops), SQA.NonEqualPair[])
        # Accumulate: terms collapsing to the same ops must add their coefficients.
        out[new_term] = get(out, new_term, zero(c)) + c
    end
    return SQA.QAdd(out, q.indices)
end

"""
Canonical label used when scaling to identical-atom ensembles: `canon_key` further
reduced under the permutation symmetry of the `selected` atomic subspaces (the rest
stay free indices). Hermitian conjugates are not identified here; the scale step folds
those separately through `canonical_rep`. This is `_treatment_key` with the `selected`
subspaces `Scaled`.
"""
function scaled_key(op::QAdd, ctx::CanonCtx; selected = ctx.symmetric)
    treatments = all_free_treatments(ctx)
    for sp in selected
        treatments[sp] = Scaled
    end
    return _treatment_key(op, ctx, treatments)
end

"""
The label matching a system's current `treatments`. Once any subspace has been fixed to
explicit sites by a previous `evaluate`, the conjugation-aware `canonical_rep` is used
so the closed set agrees with the generated code; while every subspace is still a
symbolic index the plain `_treatment_key` suffices.
"""
function _materialised_key(op::QAdd, ctx::CanonCtx, treatments::Dict{Int, SubspaceTreatment})
    any(==(Concrete), values(treatments)) && return canonical_rep(op, ctx; treatments)[1]
    return _treatment_key(op, ctx, treatments)
end
_materialised_key(op, ::CanonCtx, ::Dict{Int, SubspaceTreatment}) = op
scaled_key(op, ::CanonCtx; kw...) = op

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
                target = isempty(reps) ? idx : nth_index(reps, k)
                target == idx && continue
                sub[idx] = target
            end
        end
        cand = isempty(sub) ? op : SQA.change_index(op, sub)
        cand = _reorder_then_drop_non_equal(cand, spaces)
        key = SQA.qadd_order_key(cand)
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
function _reorder_then_drop_non_equal(op::QAdd, spaces)
    by_space = Dict(k => v for (k, v) in _free_by_space(op) if k in spaces)
    pairs = Tuple{SQA.Index, SQA.Index}[]
    for (_, idxs) in by_space, n in 1:length(idxs), m in (n + 1):length(idxs)
        push!(pairs, (idxs[n], idxs[m]))
    end
    reordered = isempty(pairs) ? op : SQA.assume_distinct_index(op, pairs)
    return _drop_all_non_equal(SQA.QAdd(reordered.arguments, SQA.Index[]))
end

"""
Conjugation-aware representative shared by `canonical_rep` and `concrete_rep`. Given a
`basekey` labelling of ⟨op⟩ (per-treatment or concrete), returns `(key, is_conjugate)`
from the smaller of the labels of ⟨`op`⟩ and ⟨`op`†⟩ under SQA's `qadd_order_key`, so an
average and its conjugate ⟨A†⟩ = ⟨A⟩* collapse to one representative. The structural
order lives in SQA (`qadd_order_key`); QC only chooses which conjugation side wins.
"""
function _conjugation_rep(op::QAdd, basekey)
    k = basekey(op)
    kc = basekey(adjoint(op))
    return SQA.qadd_order_key(kc) < SQA.qadd_order_key(k) ? (kc, true) : (k, false)
end

"""
Conjugation-aware label for materialised systems, the `concrete_key` counterpart of
`canonical_rep`: returns `(key, is_conjugate)` from the smaller of the labels of
⟨`op`⟩ and ⟨`op`†⟩.
"""
concrete_rep(op::QAdd) = _conjugation_rep(op, concrete_key)
concrete_rep(op) = (op, false)

# ---- moment lookup -----------------------------------------------------------

"""
Canonical-representative lookup mapping each moment ⟨op⟩ to a payload of type `V`, keyed by
its Hermitian-conjugate representative (`canonical_rep` under `treatments`). An average, its
conjugate, and any symmetry- or index-relabelled form resolve to the same entry; the stored
`Bool` records which conjugation side the representative came from, so a query on the
opposite side can recover the conjugate payload. The single moment↔quantity matcher shared
by the state registry, `System`, `get_solution`, the spectrum kernel and the correlation
steady-state lookup, in place of separate hand-rolled `rep → (payload, side)` dictionaries.
"""
struct MomentMap{V}
    ctx::CanonCtx
    treatments::Dict{Int, SubspaceTreatment}
    by_rep::Dict{QAdd, Tuple{V, Bool}}
end

"""
Build a `MomentMap` keying `payloads[k]` by the representative of `ops[k]`. A non-`QAdd` op is
skipped together with its paired payload; the first payload to reach a representative wins,
matching the `get!` dedup that collapses conjugate/symmetry-equivalent duplicates.
"""
function MomentMap(
        ctx::CanonCtx, treatments::Dict{Int, SubspaceTreatment}, ops, payloads::AbstractVector{V},
    ) where {V}
    by_rep = Dict{QAdd, Tuple{V, Bool}}()
    for (op, p) in zip(ops, payloads)
        op isa QAdd || continue
        rep, side = canonical_rep(op, ctx; treatments)
        get!(by_rep, rep, (p, side))
    end
    return MomentMap{V}(ctx, treatments, by_rep)
end

"""
Resolve ⟨`op`⟩ to its `(payload, same_side)` entry, `same_side` being `true` when the query and
the stored representative share a conjugation side. `nothing` for a non-`QAdd` `op` or a
missing representative.
"""
function match_moment(m::MomentMap, op)
    op isa QAdd || return nothing
    rep, side = canonical_rep(op, m.ctx; treatments = m.treatments)
    haskey(m.by_rep, rep) || return nothing
    payload, rep_side = m.by_rep[rep]
    return (payload, side == rep_side)
end

"""
Resolve ⟨`op`⟩ to its symbolic payload, conjugated (`term(conj, …)`, which does not fold to
identity on `Real` symtype) when the query sits on the opposite conjugation side from the
stored representative. `nothing` when absent. The symbolic-payload form of `match_moment`.
"""
function resolve_moment_sym(m::MomentMap, op)
    r = match_moment(m, op)
    r === nothing && return nothing
    sym, same = r
    u = SymbolicUtils.unwrap(sym)
    return same ? u : SymbolicUtils.term(conj, u; type = Number)
end
