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
        elseif _is_double_indexed_var(rk) && !isempty(_op_free_indices(J[k]))
            # Collective indexed dissipation: a `DoubleIndexedVariable` rate Γ(i,j)
            # on a singly-indexed jump gives the cross-jump dissipator Σ_{i,j} Γ D[J_i,J_j].
            acc += _collective_indexed_lindblad(op, J[k], Jdagger[k], rk)
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

# A `DoubleIndexedVariable` rate `Γ(i,j)` materialises as a callable `Sym` whose
# `FnType` domain is a 2-tuple. (A single `IndexedVariable` has a 1-tuple domain;
# a plain scalar rate is not a call.) Used to switch the dissipator onto the
# collective cross-jump form.
_fn_domain_len(::Type{<:SymbolicUtils.FnType{A}}) where {A <: Tuple} = length(A.parameters)
_fn_domain_len(::Type) = -1
function _is_double_indexed_var(x)
    u = SymbolicUtils.unwrap(x)
    u isa SymbolicUtils.BasicSymbolic || return false
    SymbolicUtils.iscall(u) || return false
    f = SymbolicUtils.operation(u)
    (f isa SymbolicUtils.BasicSymbolic && !SymbolicUtils.iscall(f)) || return false
    return _fn_domain_len(SymbolicUtils.symtype(f)) == 2 &&
        length(SymbolicUtils.arguments(u)) == 2
end

# Partner `Index` from the rate's second argument, inheriting the jump index's
# range and subspace (derive the name from the user's vocabulary, never invent one).
function _partner_index(rk, i_jump::SQA.Index)
    args = SymbolicUtils.arguments(SymbolicUtils.unwrap(rk))
    sym_j = args[2]
    return SQA.Index(Base.nameof(sym_j), i_jump.range, i_jump.space_index, Symbolics.Num(sym_j))
end

# Σ_{i,j} Γ(i,j) D[J_i, J_j] for a singly-indexed jump `Jk`, rate `rk=Γ(i,j)`.
# Standard collective-rate-matrix split: explicit diagonal self-decay
# `Σ_i Γ(i,i) D[J_i]` + off-diagonal cross-recycling `Σ_{i≠j} Γ(i,j) D[J_i,J_j]`.
# The diagonal must be explicit — completeness does NOT recover it (the off-diagonal
# emission has no ground projector to fold), so off-diagonal-only loses the self-decay.
function _collective_indexed_lindblad(op, Jk, Jdk, rk)
    i_jump = first(_op_free_indices(Jk))
    jdx = _partner_index(rk, i_jump)
    Jj = SQA.change_index(Jk, i_jump, jdx)
    off_term = (rk / 2) * (Jdk * commutator(op, Jj) + commutator(Jdk, op) * Jj)
    off = SQA.Σ(SQA.Σ(off_term, jdx, SQA.Index[i_jump]), i_jump)
    rkk = _diag_rate(rk)
    dia_term = (rkk / 2) * (Jdk * commutator(op, Jk) + commutator(Jdk, op) * Jk)
    dia = SQA.Σ(dia_term, i_jump)
    return off + dia
end

# Diagonal rate `Γ(i,i)` from `Γ(i,j)` (collapse the partner onto the jump index).
function _diag_rate(rk)
    u = SymbolicUtils.unwrap(rk)
    f = SymbolicUtils.operation(u)
    a = SymbolicUtils.arguments(u)
    return f(a[1], a[1])
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

# Free LHS indices on `new_ops` that share a Hilbert subspace, deduped, but only
# when a subspace carries 2+ such indices. A multi-index moment ⟨X_i X_j⟩ is by
# construction a distinct-slot cumulant (i≠j); its diagonal (i=j) is a separate,
# lower node. Asserting the distinctness here lets SQA's diagonal split collapse
# the dissipator's same-index contribution (the `k=i,j` terms of `Σ_k D[c_k]`)
# instead of leaking spurious higher-order cumulants. This applies to BOTH atom
# (Transition) AND Fock (Destroy/Create filter/mode) subspaces: without it a
# cross-mode moment like ⟨b_i b_j⟩ grows a bogus nonlinear `κf` dissipator
# instead of the exact `-κf⟨b_i b_j⟩`. Distinctness assertion is independent of
# the scaling permutation fold (which lives in scaling.jl and never folds Fock
# modes); it only fixes operator order and feeds the diagonal split.
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
