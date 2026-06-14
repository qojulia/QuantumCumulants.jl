# ---- limits parsing ----------------------------------------------------------

"""
Parse the user's `limits` (a `Pair`, a tuple of `Pair`s, or a `Dict`) into a map from
each symbolic range bound to its concrete integer size (the number of atoms/sites).
"""
function _limits_dict(limits)
    sub = Dict{Any, Int}()
    limits === nothing && return sub
    if limits isa Pair
        sub[SymbolicUtils.unwrap(first(limits))] = Int(last(limits))
    elseif limits isa Tuple || limits isa AbstractDict
        for p in limits
            sub[SymbolicUtils.unwrap(first(p))] = Int(last(p))
        end
    else
        throw(ArgumentError("`limits` must be a Pair, Tuple of Pairs, or Dict"))
    end
    return sub
end

_in_h(hset::Set{Int}, sp::Int) = isempty(hset) || sp in hset
_targeted(idx::SQA.Index, sub, hset) =
    SymbolicUtils.unwrap(idx.range) in keys(sub) && _in_h(hset, idx.space_index)

"""
Mint a concrete index from the subspace's first-declared vocabulary index, suffixed
with the position `k` (SQA naming policy: trace back to the user's vocabulary).
"""
function _concrete_index(b::SQA.Index, k::Int, ctx::CanonCtx)
    reps = get(ctx.vocab, b.space_index, SQA.Index[])
    isempty(reps) && return b(k)
    return reps[1](k)
end

"""Free indices on `op` whose range is targeted by `limits` and whose subspace is selected by `h`."""
_free_limited_indices(op::SQA.QSym, sub, hset) =
    (SQA.has_index(op.index) && _targeted(op.index, sub, hset)) ? SQA.Index[op.index] : SQA.Index[]
function _free_limited_indices(op::QAdd, sub, hset)
    out = SQA.Index[]
    bound = Set(op.indices)
    for (term, _) in op.arguments, o in term.ops
        idx = o.index
        SQA.has_index(idx) || continue
        idx in bound && continue
        _targeted(idx, sub, hset) || continue
        idx in out || push!(out, idx)
    end
    return out
end
_free_limited_indices(_, _, _) = SQA.Index[]

"""
Distinct-site cumulant semantics: two free indices on the same subspace within one
moment are distinct (i ≠ j), so only injective slot assignments are enumerated. The
diagonal is a lower-order moment owned by a separate node, reached via the drift's own
diagonal split rather than enumerated here.
"""
function _distinct_within_subspace(free::Vector{SQA.Index}, tup)
    @inbounds for p in 1:length(free), q in (p + 1):length(free)
        free[p].space_index == free[q].space_index || continue
        tup[p] == tup[q] && return false
    end
    return true
end

_used_op_indices(op::QAdd) = _free_op_indices(op)   # all op indices (incl. bound), first-encounter

# ---- materialisation ---------------------------------------------------------

_apply_free(op::QAdd, idx_sub) = isempty(idx_sub) ? op :
    foldl((q, p) -> SQA.change_index(q, p[1], p[2]), idx_sub; init = op)
_apply_free(op, _) = op

"""
Rewrite an averaged RHS for the concrete system, walking the expression tree: average
leaves go through `_materialise_leaf`, bare leaf `Sym`s through `sym_sub` (the
Symbolics-side mirror of the operator `change_index`), and a scoped `*` node through
`_materialise_scoped` so an `IndexedVariable` coefficient `g(i)` unrolls in lockstep
with its sum leaf.
"""
_materialise(x, idx_sub, sub, ctx, hset, sym_sub = _EMPTY_SYM_SUB) =
    Symbolics.Num(_materialise_walk(SymbolicUtils.unwrap(x), idx_sub, sub, ctx, hset, sym_sub))

const _EMPTY_SYM_SUB = Dict{Any, Any}()

function _materialise_walk(x, idx_sub, sub, ctx, hset, sym_sub)
    return rewrite(x; descend = SymbolicUtils.iscall) do y
        _is_avg_leaf(y) && return _materialise_leaf(y, idx_sub, sub, ctx, hset)
        # leaf Sym: mirror the operator index rename
        SymbolicUtils.iscall(y) || return get(sym_sub, y, y)
        if SymbolicUtils.operation(y) === (*) && SymbolicUtils.hasmetadata(y, SQA.SumIndices)
            bound = SymbolicUtils.getmetadata(y, SQA.SumIndices)
            bound isa Vector{SQA.Index} && !isempty(bound) &&
                return _materialise_scoped(y, bound, idx_sub, sub, ctx, hset, sym_sub)
        end
        return nothing
    end
end

"""
Unroll a scoped `*` node over the bound indices it references (those targeted by
`sub`/`h` and not already in `idx_sub`). Each enumeration extends `idx_sub` (operator
side) and `sym_sub` (coefficient side) before recursing into the scope-stripped subtree.
"""
function _materialise_scoped(x, bound::Vector{SQA.Index}, idx_sub, sub, ctx, hset, sym_sub)
    substituted = Set{SQA.Index}()
    for (from, to) in idx_sub
        push!(substituted, from); push!(substituted, to)
    end
    to_unroll = SQA.Index[
        b for b in bound
            if SymbolicUtils.unwrap(b.range) in keys(sub) &&
            _in_h(hset, b.space_index) &&
            !(b in substituted) &&
            _references_index(x, b)
    ]
    stripped = SymbolicUtils.setmetadata(x, SQA.SumIndices, SQA.Index[])
    isempty(to_unroll) && return _materialise_walk(stripped, idx_sub, sub, ctx, hset, sym_sub)
    ranges = Int[sub[SymbolicUtils.unwrap(b.range)] for b in to_unroll]
    terms = Any[]
    for tup in Iterators.product((1:r for r in ranges)...)
        local_sub = copy(idx_sub)
        local_sym = copy(sym_sub)
        for (b, v) in zip(to_unroll, tup)
            fresh = _concrete_index(b, v, ctx)
            local_sub[b] = fresh
            local_sym[SymbolicUtils.unwrap(b.sym)] = SymbolicUtils.unwrap(fresh.sym)
        end
        push!(terms, _materialise_walk(stripped, local_sub, sub, ctx, hset, local_sym))
    end
    isempty(terms) && return 0
    length(terms) == 1 && return terms[1]
    # Build the Add via maketerm: `+` enforces a numeric-symtype guard that rejects
    # the Any-typed `g(i.sym)` callable children.
    proto = findfirst(t -> t isa SymbolicUtils.BasicSymbolic, terms)
    proto === nothing && return sum(terms)
    return TermInterface.maketerm(typeof(terms[proto]), +, terms, nothing)
end

"""
True if any leaf in `x` references `b`, via the Symbolics variable `b.sym` or a QField
carrying `b`. Ignores `SumIndices` metadata, which SQA stamps uniformly across a scoped
sum's terms even where they do not depend on `b`.
"""
function _references_index(x, b::SQA.Index)
    x isa SymbolicUtils.BasicSymbolic || return false
    if _is_avg_leaf(x)
        return _qfield_references_index(undo_average(x), b)
    end
    SymbolicUtils.iscall(x) && return any(a -> _references_index(a, b), SymbolicUtils.arguments(x))
    return isequal(x, SymbolicUtils.unwrap(b.sym))
end

function _qfield_references_index(op::QAdd, b::SQA.Index)
    for (term, _) in op.arguments, o in term.ops
        SQA.has_index(o.index) && o.index === b && return true
    end
    return false
end
_qfield_references_index(op::SQA.QSym, b::SQA.Index) = SQA.has_index(op.index) && op.index === b
_qfield_references_index(_, ::SQA.Index) = false

"""
Lift `SumIndices`/`SumNonEqual` metadata from an average-leaf factor up onto its
enclosing `*` node, so a coefficient sibling (`g(i)`) unrolls together with the leaf.
"""
function _lift_scope_node(y)
    (SymbolicUtils.iscall(y) && SymbolicUtils.operation(y) === (*)) || return y
    for a in SymbolicUtils.arguments(y)
        a isa SymbolicUtils.BasicSymbolic || continue
        SymbolicUtils.hasmetadata(a, SQA.SumIndices) || continue
        y = SymbolicUtils.setmetadata(y, SQA.SumIndices, SymbolicUtils.getmetadata(a, SQA.SumIndices))
        SymbolicUtils.hasmetadata(a, SQA.SumNonEqual) &&
            (y = SymbolicUtils.setmetadata(y, SQA.SumNonEqual, SymbolicUtils.getmetadata(a, SQA.SumNonEqual)))
        break
    end
    return y
end
_lift_sum_scope(x) = rewrite(
    Returns(nothing), x; descend = SymbolicUtils.iscall,
    post = _lift_scope_node, maketerm = _structural_maketerm
)

"""
Materialise one average leaf ⟨op⟩: apply the free-index substitution, drop targeted sum
scope that is not being unrolled (it is no longer a real sum), then expand each remaining
bound index over its concrete range `1:n`. Non-`QAdd` leaves pass through.
"""
function _materialise_leaf(avg, idx_sub, sub, ctx, hset)
    op = undo_average(avg)
    op isa QAdd || return avg
    op2 = _apply_free(op, idx_sub)
    targets = Set(values(idx_sub))
    used = _used_op_indices(op2)
    bound = SQA.Index[
        b for b in op2.indices
            if _targeted(b, sub, hset) && !(b in targets) && (b in used)
    ]
    # Strip targeted scope that is NOT being unrolled (no longer a real sum).
    kept = SQA.Index[
        b for b in op2.indices
            if !(_targeted(b, sub, hset)) || (b in bound)
    ]
    length(kept) == length(op2.indices) || (op2 = QAdd(op2.arguments, kept))
    isempty(bound) && return average(op2)
    cur = op2
    for b in bound
        cur = _unroll_one(cur, b, sub[SymbolicUtils.unwrap(b.range)], ctx)
    end
    return cur
end

"""
Unroll one bound index over `1:n`: strip it from `.indices` and sum
`change_index(bare, b -> concrete_k)` over `k` (SQA handles projector squashing,
commuting-op sorting and non-equal-violated drops on the now scope-free `QAdd`).
"""
_unroll_one(x::Number, _, _, _) = x
function _unroll_one(qadd::QAdd, b::SQA.Index, n::Int, ctx)
    bare = QAdd(qadd.arguments, SQA.Index[i for i in qadd.indices if i != b])
    terms = Any[average(SQA.change_index(bare, b, _concrete_index(b, k, ctx))) for k in 1:n]
    isempty(terms) && return 0
    return reduce(+, terms)
end
function _unroll_one(x::SymbolicUtils.BasicSymbolic, b::SQA.Index, n::Int, ctx)
    return rewrite(x) do y
        _is_avg_leaf(y) || return nothing
        op = undo_average(y)
        (op isa QAdd && b in op.indices) || return y
        return _unroll_one(op, b, n, ctx)
    end
end

_is_zero_qadd(op::QAdd) = isempty(op.arguments)
_is_zero_qadd(_) = false

# ---- the pass ----------------------------------------------------------------

"""
Build a `MomentGraph` from a completed equation set using its stored drifts. Reuses the
recorded RHSs (preserving any filter applied during completion) rather than re-deriving.
"""
function _graph_from_stored(eqs::AbstractMeanfieldEquations)
    ctx = build_ctx(eqs)
    eff = eqs isa NoiseMeanfieldEquations ? eqs.efficiencies : nothing
    sys = SystemSpec(
        eqs.hamiltonian, collect(eqs.jumps), collect(eqs.jumps_dagger),
        collect(eqs.rates), eff, eqs.iv, eqs.order, maximum, eqs.direction,
    )
    treatments = _treatments(eqs, ctx)
    nodes = OrderedCollections.OrderedDict{NodeKey, NodeData}()
    for i in eachindex(eqs.operators)
        op = eqs.operators[i]
        opqa = op isa QAdd ? op : op * 1
        # Key in the STORED treatment, not bare `canon_key`: a subspace made Concrete
        # by a prior `evaluate` must not be alpha-renamed and collapsed.
        k = _materialised_key(opqa, ctx, treatments)
        haskey(nodes, k) && continue
        drift = Symbolics.Num(_lift_sum_scope(SymbolicUtils.unwrap(eqs.equations[i].rhs)))
        op_drift = opqa
        noise = eff === nothing ? nothing :
            Symbolics.Num(_lift_sum_scope(SymbolicUtils.unwrap(eqs.noise_equations[i].rhs)))
        nodes[k] = NodeData(drift, op_drift, noise, nothing, get_order(opqa), SQA.acts_on(opqa))
    end
    return MomentGraph(nodes, sys, ctx, treatments)
end

"""
Specialise a symbolic `MomentGraph` to a concrete number of atoms/sites. For each node,
enumerate its free indices over the ranges in `limits` (distinct sites, `i ≠ j`),
materialise the drift (and noise) at each assignment, and keep one node per Hermitian
conjugate pair (⟨X†⟩ = conj⟨X⟩). Every unrolled Hilbert subspace is marked `Concrete`;
`h` restricts the expansion to selected subspaces. A no-op when `limits` is empty.
"""
function specialize(g::MomentGraph, limits; h::Vector{Int} = Int[])
    sub = _limits_dict(limits)
    isempty(sub) && return g
    ctx = g.ctx
    hset = Set{Int}(h)
    nodes = OrderedCollections.OrderedDict{NodeKey, NodeData}()
    seen = Set{QAdd}()
    for (k, nd) in g.nodes
        free = _free_limited_indices(k, sub, hset)
        ranges = Int[sub[SymbolicUtils.unwrap(idx.range)] for idx in free]
        for tup in Iterators.product((1:r for r in ranges)...)
            _distinct_within_subspace(free, tup) || continue
            idx_sub = Dict{SQA.Index, SQA.Index}(
                b => _concrete_index(b, v, ctx) for (b, v) in zip(free, tup)
            )
            new_k = _apply_free(k, idx_sub)
            _is_zero_qadd(new_k) && continue
            # Deduplicate by the Hermitian conjugation (`concrete_rep`): `⟨X†⟩ = conj⟨X⟩` is an
            # exact redundancy, so keep one rep per conjugate pair (resolver recovers the partner).
            rep, _ = concrete_rep(new_k)
            rep in seen && continue
            push!(seen, rep)
            lk = concrete_key(new_k)
            # Symbolics-side mirror of the free-index substitution, so coefficient
            # references (e.g. `g(i)`) rename in lockstep with the operators.
            sym_sub = isempty(idx_sub) ? _EMPTY_SYM_SUB :
                Dict{Any, Any}(
                    SymbolicUtils.unwrap(b.sym) => SymbolicUtils.unwrap(t.sym)
                    for (b, t) in idx_sub
                )
            drift = _materialise(nd.drift, idx_sub, sub, ctx, hset, sym_sub)
            noise = nd.noise === nothing ? nothing : _materialise(nd.noise, idx_sub, sub, ctx, hset, sym_sub)
            nodes[lk] = NodeData(drift, _apply_free(nd.op_drift, idx_sub), noise, nd.op_noise, nd.order, nd.aon)
        end
    end
    # A subspace is now Concrete if any of its vocabulary indices' ranges were
    # targeted by `limits` and the subspace is selected by `h`.
    treatments = copy(g.treatments)
    for (sp, idxs) in ctx.vocab
        _in_h(hset, sp) || continue
        any(idx -> SymbolicUtils.unwrap(idx.range) in keys(sub), idxs) && (treatments[sp] = Concrete)
    end
    return MomentGraph(nodes, g.sys, ctx, treatments)
end

# ---- indexed-variable -> Symbolics-array materialisation ---------------------

_is_fntype(::Type) = false
_is_fntype(::Type{<:SymbolicUtils.FnType}) = true

"""Bucket every callable-`Sym` use `f(slot_syms...)` in `x` by function name and slot tuple."""
function _collect_callable_uses!(uses, x)
    walk(x) do n
        SymbolicUtils.iscall(n) || return true
        op = SymbolicUtils.operation(n)
        args = SymbolicUtils.arguments(n)
        if op isa SymbolicUtils.BasicSymbolic && !SymbolicUtils.iscall(op) &&
                _is_fntype(SymbolicUtils.symtype(op)) &&
                all(a -> a isa SymbolicUtils.BasicSymbolic && !SymbolicUtils.iscall(a), args)
            slots = [SQA.index_slot(a) for a in args]
            if all(!isnothing, slots)
                bucket = get!(() -> Dict{Vector{Int}, Vector{SymbolicUtils.BasicSymbolic}}(), uses, Base.nameof(op))
                push!(get!(bucket, Int[s for s in slots], SymbolicUtils.BasicSymbolic[]), n)
                return false
            end
        end
        return true
    end
    return
end

"""
Build a substitution replacing the callable `Sym`s left after unrolling (e.g.
`g(i_2_3)`, an `IndexedVariable` at a concrete atom) with `getindex` into a freshly
minted Symbolics array parameter of the right shape. MTK otherwise trips a symtype
assertion on the bare callables; `parameter_map(eqs, pairs)` later matches user values
to the array by name.
"""
function _build_callable_to_array_sub(eqs::Vector{Symbolics.Equation}, states)
    uses = Dict{Symbol, Dict{Vector{Int}, Vector{SymbolicUtils.BasicSymbolic}}}()
    for eq in eqs
        _collect_callable_uses!(uses, SymbolicUtils.unwrap(eq.lhs))
        _collect_callable_uses!(uses, SymbolicUtils.unwrap(eq.rhs))
    end
    for s in states
        _collect_callable_uses!(uses, SymbolicUtils.unwrap(s))
    end
    isempty(uses) && return Dict{Any, Any}()
    sub = Dict{Any, Any}()
    for (name, slot_map) in uses
        n_args = length(first(keys(slot_map)))
        max_per_dim = [maximum(k[d] for k in keys(slot_map)) for d in 1:n_args]
        # Mint a PROPER Symbolics array variable (symtype `Vector`/`Matrix{Real}`) so
        # MTK scalarizes and binds it; a scalar `Sym` with `getindex` cannot scalarize.
        m1 = max_per_dim[1]
        arr = if n_args == 1
            SymbolicUtils.unwrap(first(@variables $(name)[1:m1]))
        elseif n_args == 2
            m2 = max_per_dim[2]
            SymbolicUtils.unwrap(first(@variables $(name)[1:m1, 1:m2]))
        else
            throw(ArgumentError("indexed coefficient with $n_args indices is not supported"))
        end
        for (slots, terms) in slot_map
            # `type = Real` so the element stays Number-symtyped (a bare `maketerm`
            # getindex infers `Any` and breaks the `Num` wrap).
            indexed = SymbolicUtils.term(getindex, arr, slots...; type = Real)
            for t in terms
                sub[t] = indexed
            end
        end
    end
    return sub
end

"""Replace exact-match callable subtrees in `x` with their `getindex` array form."""
_apply_callable_sub(x, sub) = Symbolics.Num(_apply_callable_walk(SymbolicUtils.unwrap(x), sub))
_apply_callable_walk(x, sub) = _subtree_substitute(x, sub)

"""
Rewrite the per-site `IndexedVariable` coefficients left after unrolling (e.g. `g(i_2)`,
a coupling at one concrete atom) into `getindex` on a minted Symbolics array parameter,
in place, so MTK can scalarise and bind them. A no-op when none remain.
"""
function _arrayize_indexed_params!(eqs::AbstractMeanfieldEquations)
    arr_sub = _build_callable_to_array_sub(eqs.equations, eqs.states)
    isempty(arr_sub) && return eqs
    for k in eachindex(eqs.equations)
        eqs.equations[k] = _apply_callable_sub(eqs.equations[k].lhs, arr_sub) ~
            _apply_callable_sub(eqs.equations[k].rhs, arr_sub)
    end
    if eqs isa NoiseMeanfieldEquations
        for k in eachindex(eqs.noise_equations)
            eqs.noise_equations[k] = _apply_callable_sub(eqs.noise_equations[k].lhs, arr_sub) ~
                _apply_callable_sub(eqs.noise_equations[k].rhs, arr_sub)
        end
    end
    return eqs
end

"""
    evaluate(eqs::AbstractMeanfieldEquations; limits=nothing, h=Int[])
    evaluate(c::CorrelationFunction; limits=nothing, h=Int[])

Materialise a symbolic indexed system into a concrete, fixed-size one by inserting the
index ranges and unrolling the sums. The result is ready for
`ModelingToolkitBase.System`; its operators carry no remaining sum scope.

# Keyword arguments
* `limits=nothing`: the concrete size of each symbolic index range, given as a `Pair`
  (`N => 3`), a tuple of pairs, or a `Dict`. Required whenever a range bound is symbolic.
* `h=Int[]`: restrict unrolling to these Hilbert subspaces (by `acts_on` index). Empty
  unrolls every subspace covered by `limits`.
"""
function evaluate(eqs::AbstractMeanfieldEquations; limits = nothing, h::Vector{Int} = Int[], kwargs...)
    sub = _limits_dict(limits)
    isempty(sub) && return _copy(eqs)
    out = assemble_equations(specialize(_graph_from_stored(eqs), limits; h))
    return _arrayize_indexed_params!(out)
end
