"""Subspaces participating in the quotient: symmetric ones, restricted by `h` (empty = all)."""
_scale_selected(ctx::CanonCtx, h::Vector{Int}) =
    isempty(h) ? ctx.symmetric : intersect(ctx.symmetric, Set(h))

function quotient(g::MomentGraph; h::Vector{Int} = Int[])
    ctx = g.ctx
    selected = _scale_selected(ctx, h)
    sym_to_space = _build_sym_to_space(g)
    # Key nodes in the graph's treatment with selected subspaces now Scaled (not bare
    # `scaled_key`, which would alpha-rename already-Concrete subspaces and collapse them).
    treatments = copy(g.treatments)
    for sp in selected
        treatments[sp] = Scaled
    end
    nodes = OrderedCollections.OrderedDict{NodeKey, NodeData}()
    for (k, nd) in g.nodes
        ok = _materialised_key(k, ctx, treatments)
        haskey(nodes, ok) && continue   # symmetric image / folded conjugate already kept
        drift = _scale_expr(nd.drift, ctx, selected, sym_to_space)
        noise = nd.noise === nothing ? nothing : _scale_expr(nd.noise, ctx, selected, sym_to_space)
        nodes[ok] = NodeData(drift, nd.op_drift, noise, nd.op_noise, nd.order, nd.aon)
    end
    return MomentGraph(nodes, g.sys, ctx, treatments)
end

"""
Symbol -> space_index map so an `IndexedVariable(:g, i)` coefficient can be traced to
its source subspace. Built from the ctx vocabulary AND the indices actually present on the
graph's moments: closure mints extra per-subspace indices (e.g. a second `i_2_2`) for
two-body moments that the frozen ctx vocab does not list, and an unmapped coefficient index
would otherwise be wrongly flattened (`_indexed_var_in_h` treats unknown as in-`h`).
"""
function _build_sym_to_space(g::MomentGraph)
    out = Dict{Symbol, Int}()
    for (sp, indices) in g.ctx.vocab, idx in indices
        out[idx.name] = sp
    end
    for k in keys(g.nodes), idx in SQA.get_indices(k)
        out[idx.name] = idx.space_index
    end
    return out
end

"""
Rewrite an averaged RHS for scale: each leaf ⟨X⟩ becomes `prefactor * ⟨X_rep⟩`, where
`X_rep` is the permutation-symmetry representative of `X` (via `_scale_leaf`), then
flatten every `IndexedVariable(:g, i)` coefficient to its scalar `g` (all atoms share
the same coupling under the scale symmetry).
"""
function _scale_expr(x, ctx, selected, sym_to_space)
    # Graph drifts keep sum scope on each leaf, not on the enclosing node, so this strip of
    # any `*`-node scope is a no-op safeguard; `_scale_leaf` owns the leaf scope and the
    # prefactor it implies.
    stripped = _strip_mul_sum_scope(SymbolicUtils.unwrap(x))
    treatments = Dict{Int, SubspaceTreatment}()
    for sp in ctx.symmetric
        treatments[sp] = sp in selected ? Scaled : Concrete
    end
    reduced = mapleaves(l -> _scale_leaf(l, ctx, selected, treatments), stripped)
    return Symbolics.Num(_flatten_indexed_vars_in_tree(reduced, selected, sym_to_space))
end

"""
Remove the `SumIndices`/`SumNonEqual` metadata `_propagate_sum_scope` placed on non-leaf
`*`/`+` nodes, leaving each average leaf's own scope intact.
"""
function _strip_scope_node(y)
    if SymbolicUtils.hasmetadata(y, SQA.SumIndices)
        y = SymbolicUtils.setmetadata(y, SQA.SumIndices, SQA.Index[])
        SymbolicUtils.hasmetadata(y, SQA.SumNonEqual) &&
            (y = SymbolicUtils.setmetadata(y, SQA.SumNonEqual, SQA.NonEqualPair[]))
    end
    return y
end
_strip_mul_sum_scope(x) =
    rewrite(Returns(nothing), x; descend = _descend, post = _strip_scope_node, maketerm = _structural_maketerm)

function _scale_leaf(avg, ctx, selected, treatments)
    op = undo_average(avg)
    op isa QAdd || return avg
    pref = _sum_scope_prefactor(op, selected)
    # Fold only the selected subspaces (Scaled); keep the rest verbatim (Concrete).
    # `_treatment_key` empties `.indices`, so re-attach the non-selected bound indices
    # afterwards to preserve their `Σ`.
    keep_idx = SQA.Index[b for b in op.indices if !(b.space_index in selected)]
    folded = _treatment_key(op, ctx, treatments)
    # `_treatment_key` drops all non-equal index constraints, correct for the selected
    # subspace but wrong for a non-selected one whose off-diagonal constraint is
    # physical; re-attach the non-equal index pairs that involve any non-selected index.
    if folded isa QAdd && !isempty(keep_idx)
        kept_non_equal = SQA.NonEqualPair[]
        for (term, _) in op.arguments, p in term.ne
            (p[1].space_index in selected && p[2].space_index in selected) && continue
            p in kept_non_equal || push!(kept_non_equal, p)
        end
        reduced_op = isempty(kept_non_equal) ? SQA.QAdd(folded.arguments, keep_idx) :
            _reattach_non_equal(folded, kept_non_equal, keep_idx)
    else
        reduced_op = folded
    end
    reduced = average(reduced_op)
    pref === 1 && return reduced
    return SymbolicUtils.unwrap(pref) * reduced
end

"""
Re-attach non-equal index pairs to a folded QAdd, restricting each to terms whose ops carry both
its indices, and set the kept bound indices as the new sum scope.
"""
function _reattach_non_equal(folded::QAdd, kept_non_equal, keep_idx)
    out = SQA.QTermDict()
    for (term, c) in folded.arguments
        present = Set{SQA.Index}()
        for o in term.ops
            SQA.has_index(o.index) && push!(present, o.index)
        end
        non_equal = SQA.NonEqualPair[p for p in kept_non_equal if p[1] in present && p[2] in present]
        out[SQA.QTerm(copy(term.ops), vcat(term.ne, non_equal))] = c
    end
    return SQA.QAdd(out, keep_idx)
end

"""
Sum-scope prefactor: each bound index on a selected subspace that appears in some
`term.ops` contributes `(b.range - count_non_equal_involving_b)`; others contribute 1.
"""
function _sum_scope_prefactor(op::QAdd, selected::Set{Int})
    isempty(op.indices) && return 1
    op_indices = Set{SQA.Index}()
    for (term, _) in op.arguments, o in term.ops
        SQA.has_index(o.index) && push!(op_indices, o.index)
    end
    prefactor = 1
    for b in op.indices
        b in op_indices || continue
        b.space_index in selected || continue
        count_b = 0
        excluded = Set{SQA.Index}()
        for (term, _) in op.arguments, pair in term.ne
            if pair[1] == b
                other = pair[2]
            elseif pair[2] == b
                other = pair[1]
            else
                continue
            end
            other in excluded && continue
            count_b += 1
            push!(excluded, other)
        end
        factor = b.range - count_b
        prefactor = (prefactor isa Number && prefactor == 1) ? factor : prefactor * factor
    end
    return prefactor
end

"""
True for an `IndexedVariable` call (`g(i)`): its operation is a `Sym` of
`FnType{…, Real}`.
"""
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

"""
True when the `IndexedVariable`'s source index sits on a selected subspace;
unrecognised argument shapes default to `true`, flattening every coupling.
"""
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

"""
Walk the coefficient tree flattening `IndexedVariable`s, without descending into
`sym_average` leaves (those are the operator moments, handled by `mapleaves`).
"""
function _flatten_indexed_vars_in_tree(x, selected::Set{Int}, sym_to_space::Dict{Symbol, Int})
    return rewrite(x; descend = _descend) do y
        (_is_indexed_var(y) && _indexed_var_in_h(y, selected, sym_to_space)) ?
            _flatten_indexed_var(y) : nothing
    end
end

"""
    scale(eqs::AbstractMeanfieldEquations; h=Int[])
    scale!(eqs::AbstractMeanfieldEquations; h=Int[])

Reduce a permutation-symmetric indexed system to one representative per symmetry class,
exploiting that every atom of an indexed family acts on the system identically. `scale!`
mutates `eqs` in place; `scale` returns a new system.

# Keyword arguments
* `h=Int[]`: the Hilbert subspaces to scale, given by their `acts_on` indices. Empty
  scales every symmetric (indexed) subspace and leaves the others untouched.
"""
scale(eqs::MeanfieldEquations; h::Vector{Int} = Int[]) =
    MeanfieldEquations(quotient(eqs.graph; h))
scale(eqs::NoiseMeanfieldEquations; h::Vector{Int} = Int[]) =
    NoiseMeanfieldEquations(quotient(eqs.graph; h))

function scale!(eqs::AbstractMeanfieldEquations; h::Vector{Int} = Int[])
    eqs.graph = quotient(eqs.graph; h)
    resync!(eqs)
    return eqs
end

@doc (@doc scale) scale!
