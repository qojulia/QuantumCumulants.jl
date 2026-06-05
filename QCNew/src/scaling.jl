# Scaling (Layer 5). `scale` is the permutation-symmetry quotient: re-key every
# node on `orbit_key` (canon_key + symmetric_min over the selected symmetric
# subspaces), dedup, and rewrite each drift/noise leaf to its orbit representative
# with a sum over a symmetric orbit collapsed to `(range - |NE|)` times the rep.
# The operator reduction is `orbit_key` (which already relabels slots, reorders
# commuting ops, and folds conjugate pairs); the only scale-specific arithmetic is
# the sum-scope prefactor.

# Resolve the subspaces that participate in the quotient: the symmetric (indexed)
# subspaces, optionally restricted by the user's `h`. Empty `h` means all of them.
_scale_selected(ctx::CanonCtx, h::Vector{Int}) =
    isempty(h) ? ctx.symmetric : intersect(ctx.symmetric, Set(h))

# The graph pass. Returns a new (quotiented) MomentGraph.
function quotient(g::MomentGraph; h::Vector{Int} = Int[])
    ctx = g.ctx
    selected = _scale_selected(ctx, h)
    sym_to_space = _build_sym_to_space(ctx)
    # Result coordinate: the graph's coordinate with the selected subspaces now
    # Scaled. Key nodes in THIS coordinate (not bare `orbit_key`, which forces
    # every non-selected subspace Free and would alpha-rename an already-Concrete
    # subspace, collapsing e.g. filter modes unrolled by a prior `evaluate`). For an
    # all-Free input graph `_coord_key(·, merged)` equals `orbit_key(·; selected)`,
    # so the common scale-first path is unchanged.
    coords = copy(g.coords)
    for sp in selected
        coords[sp] = Scaled
    end
    nodes = OrderedCollections.OrderedDict{NodeKey, NodeData}()
    for (k, nd) in g.nodes
        ok = _materialised_key(k, ctx, coords)
        haskey(nodes, ok) && continue   # symmetric image / folded conjugate already kept
        drift = _scale_expr(nd.drift, ctx, selected, sym_to_space)
        noise = nd.noise === nothing ? nothing : _scale_expr(nd.noise, ctx, selected, sym_to_space)
        nodes[ok] = NodeData(drift, nd.op_drift, noise, nd.op_noise, nd.order, nd.aon)
    end
    return MomentGraph(nodes, g.sys, ctx, coords)
end

# Symbol -> space_index map from the canonical vocabulary, so an
# `IndexedVariable(:g, i)` coefficient (which carries only `i.sym`) can be
# traced to its source subspace when deciding whether to flatten it.
function _build_sym_to_space(ctx::CanonCtx)
    out = Dict{Symbol, Int}()
    for (sp, indices) in ctx.vocab, idx in indices
        out[idx.name] = sp
    end
    return out
end

# Rewrite an averaged RHS for scale: each leaf ⟨X⟩ becomes
# `prefactor * ⟨orbit_rep(X)⟩` (mapleaves over the average leaves), then flatten
# every `IndexedVariable(:g, i)` coefficient to its scalar `g` (a separate pass
# over the non-average part of the tree, since under the scale symmetry all atoms
# share the same coupling). Wrap with the `Symbolics.Num` constructor (as `derive`
# does), not `Symbolics.wrap`: a `prefactor * ⟨X⟩` product carries the average
# symtype (`Number`, not `Real`), which `wrap` would leave unwrapped and
# `NodeData`'s field convert would reject.
function _scale_expr(x, ctx, selected, sym_to_space)
    # `_graph_from_stored` lifted each leaf's `SumIndices` onto its enclosing `*`
    # node (without removing it from the leaf), so the scope lives in BOTH places.
    # `_scale_leaf` manages the scope on the leaf (collapsing selected-subspace
    # sums to a prefactor, preserving non-selected sums), so strip the duplicate
    # `*`-node scope first to avoid a doubled `Σ`.
    stripped = _strip_mul_sum_scope(SymbolicUtils.unwrap(x))
    reduced = mapleaves(l -> _scale_leaf(l, ctx, selected), stripped)
    return Symbolics.Num(_flatten_indexed_vars_in_tree(reduced, selected, sym_to_space))
end

# Remove `SumIndices`/`SumNonEqual` metadata from non-leaf `*`/`+` nodes (the
# copy that `_lift_sum_scope` placed there), leaving each average leaf's own
# scope intact. Identity comparison (metadata is invisible to `isequal`).
function _strip_mul_sum_scope(x)
    x isa SymbolicUtils.BasicSymbolic || return x
    _is_avg_leaf(x) && return x
    SymbolicUtils.iscall(x) || return x
    op = SymbolicUtils.operation(x)
    args = SymbolicUtils.arguments(x)
    new_args = Any[_strip_mul_sum_scope(a) for a in args]
    changed = any(((a, b),) -> a !== b, zip(args, new_args))
    y = changed ? TermInterface.maketerm(typeof(x), op, new_args, TermInterface.metadata(x)) : x
    if SymbolicUtils.hasmetadata(y, SQA.SumIndices)
        y = SymbolicUtils.setmetadata(y, SQA.SumIndices, SQA.Index[])
        SymbolicUtils.hasmetadata(y, SQA.SumNonEqual) &&
            (y = SymbolicUtils.setmetadata(y, SQA.SumNonEqual, SQA.NonEqualPair[]))
    end
    return y
end

function _scale_leaf(avg, ctx, selected)
    op = undo_average(avg)
    op isa QAdd || return avg
    pref = _sum_scope_prefactor(op, selected)
    # Coordinate isolation (spec Task 7a): fold ONLY the selected subspaces; every
    # OTHER subspace is left verbatim, including its sum scope. Key with a
    # coordinate where selected subspaces are Scaled (alpha-rename + symmetric_min
    # fold) and all other symmetric subspaces are Concrete (keep names). `_coord_key`
    # empties `.indices`, so re-attach the non-selected bound indices afterwards to
    # preserve their `Σ` (e.g. the filter `Σ_i` in the `⟨a†a⟩` drift, which
    # `orbit_key`/`canon_key` would otherwise strip).
    keep_idx = SQA.Index[b for b in op.indices if !(b.space_index in selected)]
    coords = Dict{Int, Coordinate}()
    for sp in ctx.symmetric
        coords[sp] = sp in selected ? Scaled : Concrete
    end
    folded = _coord_key(op, ctx, coords)
    # `_coord_key` runs `_drop_all_ne`, correct for the SELECTED (scaled) subspace
    # (its NE is part of the orbit identity) but WRONG for a non-selected one: an
    # off-diagonal sum like the filter's `Σ_{i≠i_2} ⟨b_i b_i_2†⟩` carries a physical
    # `i≠i_2` constraint that, if dropped, lets the diagonal `i=i_2` leak in on
    # `evaluate` (the spurious `b_i_2 b_i_2† = 1 + b_i_2† b_i_2` terms, a wrong
    # filter population). Non-selected subspaces are Concrete here, so `_coord_key`
    # keeps their index names and the original NE pairs still reference valid
    # indices: re-attach the pairs that involve any non-selected index.
    if folded isa QAdd && !isempty(keep_idx)
        kept_ne = SQA.NonEqualPair[]
        for (term, _) in op.arguments, p in term.ne
            (p[1].space_index in selected && p[2].space_index in selected) && continue
            p in kept_ne || push!(kept_ne, p)
        end
        reduced_op = isempty(kept_ne) ? SQA.QAdd(folded.arguments, keep_idx) :
            _reattach_ne(folded, kept_ne, keep_idx)
    else
        reduced_op = folded
    end
    reduced = average(reduced_op)
    pref === 1 && return reduced
    return SymbolicUtils.unwrap(pref) * reduced
end

# Re-attach NE pairs to a folded single-leaf QAdd, restricting each pair to the
# terms whose ops actually carry both of its indices, and set the kept bound
# indices as the new sum scope.
function _reattach_ne(folded::QAdd, kept_ne, keep_idx)
    out = SQA.QTermDict()
    for (term, c) in folded.arguments
        present = Set{SQA.Index}()
        for o in term.ops
            SQA.has_index(o.index) && push!(present, o.index)
        end
        ne = SQA.NonEqualPair[p for p in kept_ne if p[1] in present && p[2] in present]
        out[SQA.QTerm(copy(term.ops), vcat(term.ne, ne))] = c
    end
    return SQA.QAdd(out, keep_idx)
end

# Sum-scope prefactor: for each bound index `b` in `op.indices` that actually
# appears in some `term.ops` and lives on a selected subspace, contribute
# `(b.range - count_NE_involving_b)`. Bound indices not used in any op, or on an
# unselected subspace, contribute 1 (their sum was either eagerly scalarised or
# stays symbolic). Ported from the old `_sum_scope_prefactor`, with the
# `_in_h`/empty-means-all check replaced by direct membership in the resolved
# `selected` set.
function _sum_scope_prefactor(op::QAdd, selected::Set{Int})
    isempty(op.indices) && return 1
    op_indices = Set{SQA.Index}()
    for (term, _) in op.arguments, o in term.ops
        SQA.has_index(o.index) && push!(op_indices, o.index)
    end
    prefactor = 1
    used_ne = Set{Int}()
    for b in op.indices
        b in op_indices || continue
        b.space_index in selected || continue
        count_b = 0
        for (term, _) in op.arguments, (k, pair) in enumerate(term.ne)
            k in used_ne && continue
            (pair[1] == b || pair[2] == b) || continue
            count_b += 1
            push!(used_ne, k)
        end
        factor = b.range - count_b
        prefactor = (prefactor isa Number && prefactor == 1) ? factor : prefactor * factor
    end
    return prefactor
end

# IndexedVariable detection and flattening (`g(i)` -> scalar `g`). SQA
# materialises `IndexedVariable(:g, i)` as a call whose operation is a
# `Sym{SymReal}` of `FnType{Tuple{Int}, Real}`.
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

_flatten_indexed_var(x::SymbolicUtils.BasicSymbolic) =
    SymbolicUtils.Sym{SymbolicUtils.SymReal}(nameof(SymbolicUtils.operation(x)); type = Real)

# Flatten the IndexedVariable only when its source index sits on a selected
# subspace. Unrecognised arg shapes fall back to "flatten" so the default
# (all symmetric subspaces selected) collapses every coupling uniformly.
function _indexed_var_in_h(x::SymbolicUtils.BasicSymbolic, selected::Set{Int}, sym_to_space::Dict{Symbol, Int})
    args = SymbolicUtils.arguments(x)
    length(args) == 1 || return true
    a = args[1]
    a isa SymbolicUtils.BasicSymbolic || return true
    SymbolicUtils.iscall(a) && return true
    sp = get(sym_to_space, Base.nameof(a), 0)
    sp == 0 && return true
    return sp in selected
end

# Walk the coefficient tree flattening IndexedVariables; do NOT descend into
# `sym_average` leaves (those are the operator moments, handled by mapleaves).
function _flatten_indexed_vars_in_tree(x, selected::Set{Int}, sym_to_space::Dict{Symbol, Int})
    x isa SymbolicUtils.BasicSymbolic || return x
    if _is_indexed_var(x) && _indexed_var_in_h(x, selected, sym_to_space)
        return _flatten_indexed_var(x)
    end
    SymbolicUtils.iscall(x) || return x
    op = SymbolicUtils.operation(x)
    op === SQA.sym_average && return x
    args = SymbolicUtils.arguments(x)
    new_args = Any[_flatten_indexed_vars_in_tree(a, selected, sym_to_space) for a in args]
    all(((a, b),) -> isequal(a, b), zip(args, new_args)) && return x
    try
        return op(new_args...)
    catch err
        err isa MethodError || err isa ArgumentError || rethrow()
        return TermInterface.maketerm(typeof(x), op, new_args, TermInterface.metadata(x))
    end
end

"""
    scale(eqs::AbstractMeanFieldEquations; h=Int[])
    scale!(eqs::AbstractMeanFieldEquations; h=Int[])

Reduce a permutation-symmetric indexed system to one representative per symmetry
orbit. `h` selects which Hilbert subspaces participate (empty means all symmetric
subspaces). `scale!` mutates in place.
"""
# Build the graph from STORED drifts (`_graph_from_stored`), never re-derive.
# `complete(...; filter_func)` records filter substitutions (e.g. zeroing phase-
# variant averages) directly in `eqs.equations`; re-deriving from the operators
# would discard them and reintroduce untracked leaves that crash codegen. The
# quotient then re-keys those stored drifts onto orbit representatives.
scale(eqs::AbstractMeanFieldEquations; h::Vector{Int} = Int[]) =
    lower_to_eqs(quotient(_graph_from_stored(eqs); h))

scale!(eqs::AbstractMeanFieldEquations; h::Vector{Int} = Int[]) =
    _replace_contents!(eqs, scale(eqs; h))
