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

_in_h(hset::Set{Int}, sp::Integer) = isempty(hset) || sp in hset
_targeted(idx::SQA.Index, sub, hset) =
    SymbolicUtils.unwrap(SQA.index_range(idx)) in keys(sub) && _in_h(hset, idx.space_index)

"""
Mint a concrete index from the subspace's first-declared vocabulary index, suffixed
with the position `k` (SQA naming policy: trace back to the user's vocabulary).
"""
function _concrete_index(b::SQA.Index, k::Int, ctx::CanonCtx)
    reps = get(ctx.vocab, b.space_index, SQA.Index[])
    isempty(reps) && return b(k)
    return nth_index(reps, k)
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
function _apply_free(op::QAdd, idx_sub::AbstractDict)
    isempty(idx_sub) && return op
    _renames_any(op, idx_sub) || return op
    return SQA.change_index(op, idx_sub)
end
_apply_free(op, _) = op

function _renames_any(op::QAdd, idx_sub::AbstractDict)
    for (term, _) in op.arguments, o in term.ops
        haskey(idx_sub, o.index) && return true
    end
    return false
end

"""
Rewrite an averaged RHS for the concrete system, walking the expression tree: plain average
leaves go through `_materialise_leaf` (a free-index rename), and any node carrying an
indexed-sum over a targeted index through `_materialise_scoped`, which unrolls that index
across the whole node, the sum body's operators (`change_index`) and any sibling coefficient
sharing the index (`sym_sub`) in lockstep. Bare leaf `Sym`s fall back to `sym_sub`, the
Symbolics-side mirror of the operator rename.
"""
_materialise(x, idx_sub, sub, ctx, hset, sym_sub = _EMPTY_SYM_SUB) =
    Symbolics.Num(_materialise_walk(SymbolicUtils.unwrap(x), idx_sub, sub, ctx, hset, sym_sub))

const _EMPTY_SYM_SUB = Dict{Any, Any}()

function _materialise_walk(x, idx_sub, sub, ctx, hset, sym_sub)
    return rewrite(x; descend = SymbolicUtils.iscall) do y
        _is_avg_leaf(y) && return _materialise_leaf(y, idx_sub, sub, ctx, hset)
        bound, ne = _scoped_unroll(y, idx_sub, sub, hset)
        isempty(bound) ||
            return _materialise_scoped(y, bound, ne, idx_sub, sub, ctx, hset, sym_sub)
        # leaf Sym: mirror the operator index rename for free-index coefficients
        SymbolicUtils.iscall(y) || return get(sym_sub, y, y)
        return nothing
    end
end

"""Materialise a plain average leaf ⟨op⟩ by renaming its free indices. Non-`QAdd` leaves pass through."""
function _materialise_leaf(avg, idx_sub, sub, ctx, hset)
    op = undo_average(avg)
    op isa QAdd || return avg
    return average(_apply_free(op, idx_sub))
end

# The indexed-sum factors `y` carries directly: `y` itself when it is a sum, else its
# `*`-factors that are sums.
function _sum_factors(y::SymbolicUtils.BasicSymbolic)
    SQA.is_indexed_sum(y) && return SymbolicUtils.BasicSymbolic[y]
    (SymbolicUtils.iscall(y) && SymbolicUtils.operation(y) === (*)) || return SymbolicUtils.BasicSymbolic[]
    return SymbolicUtils.BasicSymbolic[
        a for a in SymbolicUtils.arguments(y) if a isa SymbolicUtils.BasicSymbolic && SQA.is_indexed_sum(a)
    ]
end
_sum_factors(_) = SymbolicUtils.BasicSymbolic[]

"""
Targeted summation indices `y` carries (the scope of its indexed-sum factors, restricted to
ranges named by `limits`, subspaces selected by `h`, and indices not already pinned by
`idx_sub`), with the non-equal constraints those sums carry.
"""
function _scoped_unroll(y, idx_sub, sub, hset)
    sums = _sum_factors(y)
    isempty(sums) && return (SQA.Index[], Tuple{SQA.Index, SQA.Index}[])
    pinned = Set{SQA.Index}()
    for (f, t) in idx_sub
        push!(pinned, f); push!(pinned, t)
    end
    bound = SQA.Index[]
    ne = Tuple{SQA.Index, SQA.Index}[]
    for s in sums
        for b in SQA.get_sum_indices(s)
            (_targeted(b, sub, hset) && !(b in pinned) && !(b in bound)) && push!(bound, b)
        end
        for p in SQA.get_sum_non_equal(s)
            p in ne || push!(ne, p)
        end
    end
    return (bound, ne)
end

"""Replace each indexed-sum factor of `y` with its summand body, exposing the bound indices inline."""
function _strip_sums(y)
    SQA.is_indexed_sum(y) && return SymbolicUtils.arguments(y)[1]
    args = SymbolicUtils.arguments(y)
    new_args = Any[
        (a isa SymbolicUtils.BasicSymbolic && SQA.is_indexed_sum(a)) ? SymbolicUtils.arguments(a)[1] : a
            for a in args
    ]
    return TermInterface.maketerm(typeof(y), SymbolicUtils.operation(y), new_args, TermInterface.metadata(y))
end

"""
Unroll the indexed-sum scope `bound` over the concrete ranges in `sub`, summing the
scope-stripped node at each assignment. Each enumeration pins the bound indices on both the
operator side (`idx_sub`, applied by `change_index` in `_materialise_leaf`) and the
coefficient side (`sym_sub`), so an index-dependent coefficient `Ω(i, k)` unrolls in lockstep
with its sum. Assignments violating a non-equal constraint are skipped.
"""
function _materialise_scoped(y, bound, ne, idx_sub, sub, ctx, hset, sym_sub)
    stripped = _strip_sums(y)
    ranges = Int[sub[SymbolicUtils.unwrap(SQA.index_range(b))] for b in bound]
    terms = Any[]
    for tup in Iterators.product((1:r for r in ranges)...)
        _assignment_ok(bound, tup, ne, idx_sub) || continue
        local_sub = copy(idx_sub)
        local_sym = copy(sym_sub)
        for (b, v) in zip(bound, tup)
            fresh = _concrete_index(b, v, ctx)
            local_sub[b] = fresh
            local_sym[SymbolicUtils.unwrap(SQA.index_sym(b))] = SymbolicUtils.unwrap(SQA.index_sym(fresh))
        end
        push!(terms, _materialise_walk(stripped, local_sub, sub, ctx, hset, local_sym))
    end
    isempty(terms) && return 0
    length(terms) == 1 && return terms[1]
    # Build the Add via maketerm: `+` enforces a numeric-symtype guard that rejects the
    # Any-typed `Ω(i, k)` callable children.
    proto = findfirst(t -> t isa SymbolicUtils.BasicSymbolic, terms)
    proto === nothing && return sum(terms)
    return TermInterface.maketerm(typeof(terms[proto]), +, terms, nothing)
end

"""
Whether a concrete `bound`-index assignment respects every non-equal constraint in `ne`,
comparing concrete positions (bound indices from `tup`, already-pinned ones from `idx_sub`).
"""
function _assignment_ok(bound, tup, ne, idx_sub)
    isempty(ne) && return true
    pos = Dict{SQA.Index, Int}()
    for (b, v) in zip(bound, tup)
        pos[b] = v
    end
    for (f, t) in idx_sub
        s = SQA.index_slot(SymbolicUtils.unwrap(SQA.index_sym(t)))
        s === nothing || (pos[f] = s)
    end
    for (a, b) in ne
        (haskey(pos, a) && haskey(pos, b) && pos[a] == pos[b]) && return false
    end
    return true
end

_is_zero_qadd(op::QAdd) = isempty(op.arguments)
_is_zero_qadd(_) = false

# ---- the pass ----------------------------------------------------------------

"""
Check that every range bound named in `limits` is a declared index range of `ctx`. A key
matching none means the user targeted a range the system does not have (a typo, or the
wrong symbol), in which case the unrolling would silently leave the system unchanged.
Lists the system's known symbolic ranges in the error.
"""
function _check_limit_keys(sub, ctx::CanonCtx)
    known = Set{Any}()
    for (_, idxs) in ctx.vocab, idx in idxs
        r = SymbolicUtils.unwrap(SQA.index_range(idx))
        r isa Number || push!(known, r)
    end
    unknown = Any[k for k in keys(sub) if !(k in known)]
    isempty(unknown) && return
    known_str = isempty(known) ? "none (the system has no symbolic index ranges)" :
        join(sort!(string.(collect(known))), ", ")
    throw(
        ArgumentError(
            "`limits` targets unknown index range(s) $(join(string.(unknown), ", ")); the \
        system's symbolic index ranges are: $known_str. Check that the range bound passed \
        to `Index(...)` matches the symbol used in `limits`.",
        )
    )
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
    _check_limit_keys(sub, g.ctx)
    ctx = g.ctx
    hset = Set{Int}(h)
    nodes = OrderedCollections.OrderedDict{NodeKey, NodeData}()
    seen = Set{QAdd}()
    for (k, nd) in g.nodes
        free = _free_limited_indices(k, sub, hset)
        ranges = Int[sub[SymbolicUtils.unwrap(SQA.index_range(idx))] for idx in free]
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
                    SymbolicUtils.unwrap(SQA.index_sym(b)) => SymbolicUtils.unwrap(SQA.index_sym(t))
                    for (b, t) in idx_sub
                )
            drift = _materialise(nd.drift, idx_sub, sub, ctx, hset, sym_sub)
            noise = nd.noise === nothing ? nothing :
                _materialise(nd.noise, idx_sub, sub, ctx, hset, sym_sub)
            nodes[lk] = NodeData(drift, _apply_free(nd.op_drift, idx_sub), noise, nd.op_noise, nd.order, nd.aon)
        end
    end
    # A subspace is now Concrete if any of its vocabulary indices' ranges were
    # targeted by `limits` and the subspace is selected by `h`.
    treatments = copy(g.treatments)
    for (sp, idxs) in ctx.vocab
        _in_h(hset, sp) || continue
        any(idx -> SymbolicUtils.unwrap(SQA.index_range(idx)) in keys(sub), idxs) && (treatments[sp] = Concrete)
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
            throw(
                ArgumentError(
                    "indexed coefficient `$name` has $n_args indices; `evaluate` only \
                materialises coefficients with 1 or 2 indices into Symbolics array \
                parameters. Reformulate the coefficient to carry at most two indices.",
                )
            )
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
Rewrite the per-site `IndexedVariable` coefficients left after unrolling into `getindex` on
a minted Symbolics array parameter, as a graph drift-rewrite. A no-op when none remain.
"""
function arrayize_graph(g::MomentGraph)
    ks = collect(keys(g.nodes))
    eqs = Symbolics.Equation[average(k) ~ g.nodes[k].drift for k in ks]
    states = SymbolicUtils.BasicSymbolic[average(k) for k in ks]
    arr_sub = _build_callable_to_array_sub(eqs, states)
    isempty(arr_sub) && return g
    return map_drifts(g, (_, d) -> _apply_callable_sub(d, arr_sub))
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
    g = arrayize_graph(specialize(eqs.graph, limits; h))
    return assemble_equations(g)
end
