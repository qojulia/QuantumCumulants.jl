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
    nodes = OrderedCollections.OrderedDict{NodeKey, NodeData}()
    for (k, nd) in g.nodes
        ok = orbit_key(k, ctx; selected)
        haskey(nodes, ok) && continue   # symmetric image / folded conjugate already kept
        drift = _scale_expr(nd.drift, ctx, selected, sym_to_space)
        noise = nd.noise === nothing ? nothing : _scale_expr(nd.noise, ctx, selected, sym_to_space)
        nodes[ok] = NodeData(drift, nd.op_drift, noise, nd.op_noise, nd.order, nd.aon)
    end
    return MomentGraph(nodes, g.sys, ctx, true)
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
    reduced = mapleaves(l -> _scale_leaf(l, ctx, selected), SymbolicUtils.unwrap(x))
    return Symbolics.Num(_flatten_indexed_vars_in_tree(reduced, selected, sym_to_space))
end

function _scale_leaf(avg, ctx, selected)
    op = undo_average(avg)
    op isa QAdd || return avg
    pref = _sum_scope_prefactor(op, selected)
    reduced = average(orbit_key(op, ctx; selected))
    pref === 1 && return reduced
    return SymbolicUtils.unwrap(pref) * reduced
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
