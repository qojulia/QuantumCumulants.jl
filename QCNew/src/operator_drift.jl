# Operator-level drift algebra (Layer 2). Pure SQA operator manipulation shared
# by `derive` (moments.jl) and the noise builders (noise.jl): coherent
# commutator, Lindblad recycling (forward and backward/retrodiction), per-jump
# index summation, the dY measurement-record correction, and the cross-atom
# distinctness fold. No equation containers, no orchestration.

# Build the operator-equation RHS for either evolution direction. Forward uses
# `+i[H,·]` and the standard Lindblad recycling; Backward uses `−i[H,·]` with the
# adjoint Lindblad recycling (`J ↔ J†`) and the trace-preserving term.
_operator_rhs(::Forward, op, imH, J, Jdagger, rates) =
    commutator(imH, op) + _lindblad_rhs(op, J, Jdagger, rates)

function _operator_rhs(::Backward, op, imH, J, Jdagger, rates)
    isempty(J) && return commutator(-imH, op)
    trace_term = Symbolics.Num(
        sum(
            rates[i] * (average(Jdagger[i] * J[i]) - average(J[i] * Jdagger[i]))
                for i in eachindex(J)
        )
    )
    return commutator(-imH, op)
    +_master_lindblad_backward(op, J, Jdagger, rates)
    +op * trace_term
end

function _lindblad_rhs(op, J, Jdagger, rates)
    isempty(J) && return zero(op)
    op_idx = _op_free_indices(op)
    acc = zero(op)
    @inbounds for k in eachindex(J)
        rk = rates[k]
        if rk isa AbstractMatrix
            # Collective decay: `J[k]`/`Jdagger[k]` are vectors of mode operators;
            # `rk[i,j]` is the cross-rate between modes `i` and `j`.
            Jv, Jdv = J[k], Jdagger[k]
            for i in eachindex(Jv), j in eachindex(Jv)
                acc += (rk[i, j] / 2) * (
                    Jdv[i] * commutator(op, Jv[j]) +
                        commutator(Jdv[i], op) * Jv[j]
                )
            end
        else
            term = (rk / 2) * (
                Jdagger[k] * commutator(op, J[k]) +
                    commutator(Jdagger[k], op) * J[k]
            )
            acc += _sum_over_jump_indices(term, J[k], op_idx)
        end
    end
    return acc
end

# Free indices of an operator not bound by a sum scope. An indexed jump σ_i^{21}
# carries a free `Index` i; wrapping the per-jump dissipator in Σ_i produces
# independent-decay (rather than collective) semantics and lets SQA's diagonal
# split fire.
_op_free_indices(op::SQA.QSym) =
    SQA.has_index(op.index) ? SQA.Index[op.index] : SQA.Index[]
function _op_free_indices(op::QAdd)
    bound = Set(op.indices)
    free = SQA.Index[]
    for term in keys(op.arguments), o in term.ops
        SQA.has_index(o.index) || continue
        o.index in bound && continue
        o.index in free && continue
        push!(free, o.index)
    end
    return free
end
_op_free_indices(_) = SQA.Index[]

# Wrap `term` in Σ over each free index that originates in the jump and is not
# already an LHS observable index. Returns `term` unchanged when there are no
# jump-specific free indices.
function _sum_over_jump_indices(term, jump, op_idx)
    jump_idx = _op_free_indices(jump)
    isempty(jump_idx) && return term
    free = SQA.Index[i for i in jump_idx if !(i in op_idx)]
    isempty(free) && return term
    return SQA.Σ(term, free[1], free[2:end]...)
end

# Adjoint-action Lindblad recycling for the backward Heisenberg/Kalman
# retrodiction picture.
function _master_lindblad_backward(op, J, Jdagger, rates)
    isempty(J) && return zero(op)
    op_idx = _op_free_indices(op)
    acc = zero(op)
    @inbounds for k in eachindex(J)
        rk = rates[k]
        rk isa AbstractMatrix && error(
            "Nondiagonal measurements not supported in backward retrodiction",
        )
        term = (-rk / 2) * op * Jdagger[k] * J[k] +
            (-rk / 2) * Jdagger[k] * J[k] * op +
            rk * J[k] * op * Jdagger[k]
        acc += _sum_over_jump_indices(term, J[k], op_idx)
    end
    return acc
end

# Extra deterministic term arising when the SDE is recast in terms of the
# measurement record `dY` instead of the Wiener increment `dW`.
function _dY_dS_extra_term(op, J, Jdagger, rates)
    out = 0
    @inbounds for k in eachindex(J)
        iszero(rates[k]) && continue
        rk = rates[k]
        rk isa AbstractMatrix && error(
            "Nondiagonal measurements not supported in meanfield_backward",
        )
        c1 = rk * (average(J[k] * op) - average(J[k]) * average(op))
        c2 = rk * (average(op * Jdagger[k]) - average(op) * average(Jdagger[k]))
        out = out + (-(c1 + c2)) * average(Jdagger[k] + J[k])
    end
    return out
end

# Assert NE between every member of `distinct` and every other distinct
# atom-space index in `q` on the same Hilbert subspace, then route through
# `SQA.assume_distinct_index` (which propagates NE, canonicalises, expands
# completeness). Returns `q` unchanged when `distinct` is empty. Drives the
# `σ^gg = 1 - Σ σ^kk` fold for population/dephasing systems.
function _assume_distinct_atom_indices(q::SQA.QAdd, distinct::Vector{SQA.Index})
    isempty(distinct) && return q
    by_space = Dict{Int, Set{SQA.Index}}()
    for (term, _) in q.arguments, o in term.ops
        SQA.has_index(o.index) || continue
        push!(get!(by_space, o.index.space_index, Set{SQA.Index}()), o.index)
    end
    for b in q.indices
        push!(get!(by_space, b.space_index, Set{SQA.Index}()), b)
    end
    pairs = Tuple{SQA.Index, SQA.Index}[]
    seen_pairs = Set{Tuple{SQA.Index, SQA.Index}}()
    for d in distinct
        idxs = get(by_space, d.space_index, Set{SQA.Index}())
        for other in idxs
            other == d && continue
            key = d < other ? (d, other) : (other, d)
            key in seen_pairs && continue
            push!(seen_pairs, key)
            push!(pairs, (d, other))
        end
    end
    return isempty(pairs) ? q : SQA.assume_distinct_index(q, pairs)
end
_assume_distinct_atom_indices(q, _distinct, _concretes = nothing) = q

# Free LHS indices on `new_ops` living on an N-level (atom) subspace, deduped, but
# only when a subspace carries 2+ such indices (the configuration that benefits
# from same-site collapse). Restricted to Transition-carrying indices: Fock-space
# (filter/mode) indices refer to physically distinct sites by construction and
# must not be folded.
function _distinct_atom_indices(new_ops)
    by_space = Dict{Int, Vector{SQA.Index}}()
    for op in new_ops
        _collect_atom_space_indices_by_space!(by_space, op)
    end
    out = SQA.Index[]
    seen = Set{SQA.Index}()
    for (_, idxs) in by_space
        length(idxs) >= 2 || continue
        for idx in idxs
            idx in seen && continue
            push!(seen, idx)
            push!(out, idx)
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
function _collect_atom_space_indices_by_space!(by_space, op::QAdd)
    for (term, _) in op.arguments, o in term.ops
        _collect_atom_space_indices_by_space!(by_space, o)
    end
    return nothing
end

# Noise-builder dispatch (the builders live in noise.jl) and the backward dY
# averaged-drift correction (used by the noise meanfield path).
_noise_builder(::Forward) = _build_noise_equations_forward
_noise_builder(::Backward) = _build_noise_equations_backward

_avg_extra_term(::Forward, _, _, _) = nothing
function _avg_extra_term(::Backward, J, Jdagger, eff_rates)
    return op -> _dY_dS_extra_term(op, J, Jdagger, eff_rates)
end
