# Jump-channel classification and index collection (Layer 1).
#
# Small structural predicates over the dissipator/Hamiltonian that the identity
# layer (`build_ctx`) and the closure policy depend on. Kept free of any
# orchestration so `canonical.jl` can use them without an upward dependency.

# Collective-decay `jumps`/`jumps_dagger` may be nested (a vector of mode
# vectors); flatten one level so each source is a single `QField`.
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

# Indices the Liouvillian treats as bound: sum-scope `.indices` and any free
# index carried by an atom inside the source (the Liouvillian sums collective
# jumps over their carried index). `build_ctx` excludes these names from the
# canonical free-index vocabulary so `derive` does not re-clash.
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

# A dissipator carries a dephasing channel when any (flattened) jump is an
# indexed diagonal atomic transition `σ^{αα}`. This selects the closure policy
# (`ctx.population`): population/dephasing systems fold cross-atom decay via
# `σ^gg = 1 - Σ σ^kk`; concrete-site systems keep the cross-decay correction.
_is_diag_atom_jump(j::SQA.QSym) =
    (j isa SQA.Transition) && SQA.has_index(j.index) && (j.i == j.j)
function _is_diag_atom_jump(j::SQA.QAdd)
    for (term, _) in j.arguments, o in term.ops
        _is_diag_atom_jump(o) && return true
    end
    return false
end
_is_diag_atom_jump(::Any) = false

function _has_dephasing_channel(jumps)
    for j in _flatten_jumps(jumps)
        _is_diag_atom_jump(j) && return true
    end
    return false
end

# Distinct free operator indices of a QAdd, in first-encounter order. Used by the
# identity layer's alpha-rename (`canonical.jl`).
function _free_op_indices(op::SQA.QAdd)
    out = SQA.Index[]
    for (term, _) in op.arguments, o in term.ops
        SQA.has_index(o.index) || continue
        o.index in out || push!(out, o.index)
    end
    return out
end

# A *leaf* average is one produced by `average(op)` directly: head `sym_average`.
# Distinguishes a single moment ⟨X⟩ from a product/expression of averages (which
# also carries the `AvgSym` symtype but is not a leaf). Used by the cumulant and
# closure layers.
function _is_leaf_average(x::SymbolicUtils.BasicSymbolic)
    SQA.is_average(x) || return false
    SymbolicUtils.iscall(x) || return false
    return SymbolicUtils.operation(x) === SQA.sym_average
end
_is_leaf_average(::Any) = false
