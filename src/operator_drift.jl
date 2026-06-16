"""
Build the operator-equation RHS for the given evolution direction. Forward uses
`+i[H,·]` and the standard Lindblad recycling; Backward uses `−i[H,·]` with the
adjoint Lindblad recycling (`J ↔ J†`) and the trace-preserving term.
"""
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
    return commutator(-imH, op) +
        _master_lindblad_backward(op, J, Jdagger, rates) +
        op * trace_term
end

"""
Standard (forward) Lindblad recycling `Σ_k (rₖ/2)(Jₖ† [op, Jₖ] + [Jₖ†, op] Jₖ)`.
A matrix rate selects collective decay across mode-operator vectors; a
`DoubleIndexedVariable` rate on a singly-indexed jump selects the collective
cross-jump dissipator; otherwise the per-jump term is summed over its free
indices.
"""
function _lindblad_rhs(op, J, Jdagger, rates)
    isempty(J) && return zero(op)
    op_idx = _op_free_indices(op)
    acc = zero(op)
    @inbounds for k in eachindex(J)
        rk = rates[k]
        if rk isa AbstractMatrix
            # `rk[i,j]` is the cross-rate between mode operators `J[k][i]` and `J[k][j]`.
            Jv, Jdv = J[k], Jdagger[k]
            for i in eachindex(Jv), j in eachindex(Jv)
                acc += (rk[i, j] / 2) * (
                    Jdv[i] * commutator(op, Jv[j]) +
                        commutator(Jdv[i], op) * Jv[j]
                )
            end
        elseif _is_double_indexed_var(rk) && !isempty(_op_free_indices(J[k]))
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

"""
Free indices of an operator that are not bound by a sum scope. An indexed jump
`σ_i^{21}` carries a free `Index` `i`; wrapping the per-jump dissipator in `Σ_i`
gives independent-decay (rather than collective) semantics and lets SQA's diagonal
split fire.
"""
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

"""
Wrap `term` in `Σ` over each free index that originates in the jump and is not
already an LHS observable index. Returns `term` unchanged when there are no
jump-specific free indices.
"""
function _sum_over_jump_indices(term, jump, op_idx)
    jump_idx = _op_free_indices(jump)
    isempty(jump_idx) && return term
    free = SQA.Index[i for i in jump_idx if !(i in op_idx)]
    isempty(free) && return term
    return SQA.Σ(term, free[1], free[2:end]...)
end

_fn_domain_len(::Type{<:SymbolicUtils.FnType{A}}) where {A <: Tuple} = length(A.parameters)
_fn_domain_len(::Type) = -1

"""
True when `x` is a `DoubleIndexedVariable` rate `Γ(i,j)`: a callable `Sym` whose
`FnType` domain is a 2-tuple. (A single `IndexedVariable` has a 1-tuple domain; a
plain scalar rate is not a call.) Used to switch the dissipator onto the
collective cross-jump form.
"""
function _is_double_indexed_var(x)
    u = SymbolicUtils.unwrap(x)
    u isa SymbolicUtils.BasicSymbolic || return false
    SymbolicUtils.iscall(u) || return false
    f = SymbolicUtils.operation(u)
    (f isa SymbolicUtils.BasicSymbolic && !SymbolicUtils.iscall(f)) || return false
    return _fn_domain_len(SymbolicUtils.symtype(f)) == 2 &&
        length(SymbolicUtils.arguments(u)) == 2
end

"""
Partner `Index` for a collective cross-jump dissipator, taken from the rate's
second argument and inheriting the jump index's range and subspace. The name comes
from the user's existing vocabulary rather than being invented.
"""
function _partner_index(rk, i_jump::SQA.Index)
    args = SymbolicUtils.arguments(SymbolicUtils.unwrap(rk))
    sym_j = args[2]
    return SQA.Index(Base.nameof(sym_j), i_jump.range, i_jump.space_index, Symbolics.Num(sym_j))
end

"""
Collective indexed dissipator `Σ_{i,j} Γ(i,j) D[J_i, J_j]` for a singly-indexed
jump `Jk` with rate `rk = Γ(i,j)`, split as explicit diagonal self-decay
`Σ_i Γ(i,i) D[J_i]` plus off-diagonal cross-recycling `Σ_{i≠j} Γ(i,j) D[J_i,J_j]`.
The diagonal must be explicit: completeness does not recover it (the off-diagonal
emission has no ground projector to fold), so an off-diagonal-only form would lose
the self-decay.
"""
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

"""
Diagonal rate `Γ(i,i)` from `Γ(i,j)`, collapsing the partner argument onto the
jump index.
"""
function _diag_rate(rk)
    u = SymbolicUtils.unwrap(rk)
    f = SymbolicUtils.operation(u)
    a = SymbolicUtils.arguments(u)
    return f(a[1], a[1])
end

"""
Adjoint-action Lindblad recycling for the backward Heisenberg/Kalman retrodiction
picture. Matrix (nondiagonal) measurement rates are not supported here.
"""
function _master_lindblad_backward(op, J, Jdagger, rates)
    isempty(J) && return zero(op)
    op_idx = _op_free_indices(op)
    acc = zero(op)
    @inbounds for k in eachindex(J)
        rk = rates[k]
        rk isa AbstractMatrix && throw(
            ArgumentError(
                "backward retrodiction does not support nondiagonal (matrix-valued) measurement \
            rates; jump channel $k carries a matrix rate. The backward Heisenberg/Kalman \
            picture requires each monitored channel to be a single diagonal observable. \
            Pass scalar rates, or build the forward equations instead.",
            )
        )
        term = (-rk / 2) * op * Jdagger[k] * J[k] +
            (-rk / 2) * Jdagger[k] * J[k] * op +
            rk * J[k] * op * Jdagger[k]
        acc += _sum_over_jump_indices(term, J[k], op_idx)
    end
    return acc
end

"""
Extra deterministic term arising when the SDE is recast in terms of the
measurement record `dY` instead of the Wiener increment `dW`.
"""
function _dY_dS_extra_term(op, J, Jdagger, rates)
    out = 0
    @inbounds for k in eachindex(J)
        iszero(rates[k]) && continue
        rk = rates[k]
        rk isa AbstractMatrix && throw(
            ArgumentError(
                "backward retrodiction does not support nondiagonal (matrix-valued) measurement \
            rates; jump channel $k carries a matrix rate. The `dY`-record drift term of the \
            backward picture requires each monitored channel to be a single diagonal \
            observable. Pass scalar rates, or build the forward equations instead.",
            )
        )
        c1 = rk * (average(J[k] * op) - average(J[k]) * average(op))
        c2 = rk * (average(op * Jdagger[k]) - average(op) * average(Jdagger[k]))
        out = out + (-(c1 + c2)) * average(Jdagger[k] + J[k])
    end
    return out
end

"""
Assert non-equality between every member of `distinct` and every other atom-space index in
`q` on the same Hilbert subspace, then route through `SQA.assume_distinct_index`
(which propagates the non-equal constraints, canonicalises, and expands completeness). Returns `q`
unchanged when `distinct` is empty. Drives the `σ^gg = 1 - Σ σ^kk` fold for
population/dephasing systems.
"""
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

"""
Free LHS indices on `new_ops` that share a Hilbert subspace, deduplicated, returned only
for subspaces carrying two or more such indices. A multi-index moment `⟨X_i X_j⟩`
is by construction a distinct-slot cumulant (`i≠j`), its diagonal (`i=j`) being a
separate lower node. Asserting that distinctness lets SQA's diagonal split collapse
the dissipator's same-index contribution (the `k=i,j` terms of `Σ_k D[c_k]`) rather
than leaking spurious higher-order cumulants.

This applies to both atom (Transition) and Fock (Destroy/Create filter/mode)
subspaces: without it a cross-mode moment like `⟨b_i b_j⟩` grows a bogus nonlinear
`κf` dissipator instead of the exact `-κf⟨b_i b_j⟩`. The assertion is independent of
the scaling permutation fold (in scaling.jl, which never folds Fock modes); it only
fixes operator order and feeds the diagonal split.
"""
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

# ---- noise builders + dispatch ----

"""
Operator form of the measurement-backaction stochastic increment for `⟨op⟩` under jump
`J_k` with detector efficiency `eff`, cumulant-truncated to second order:

    √(η r) [ ⟨J†·op⟩ + ⟨op·J⟩ - (⟨J†⟩ + ⟨J⟩) ⟨op⟩ ]

The leading average is taken downstream; this returns the operator whose average
reproduces the expression above.
"""
function _noise_drift_one(op, J_k, Jd_k, rate, eff)
    coeff = sqrt(eff * rate)
    avg_term = Symbolics.Num(average(Jd_k + J_k))
    return coeff * (Jd_k * op + op * J_k - avg_term * op)
end

function _build_noise_equations_forward(ops, J, Jdagger, rates, efficiencies)
    n_ops = length(ops)
    n_J = length(J)
    operator_noise_eqs = Vector{Symbolics.Equation}(undef, n_ops)
    noise_eqs = Vector{Symbolics.Equation}(undef, n_ops)
    @inbounds for (i, op) in enumerate(ops)
        op_drift = 0 * op
        avg_drift = 0
        for k in 1:n_J
            iszero(efficiencies[k]) && continue
            d = SQA.expand_completeness(
                _noise_drift_one(op, J[k], Jdagger[k], rates[k], efficiencies[k])
            )
            op_drift = op_drift + d
            avg_drift = avg_drift + average(d)
        end
        operator_noise_eqs[i] = op ~ op_drift
        noise_eqs[i] = average(op) ~ avg_drift
    end
    return operator_noise_eqs, noise_eqs
end

"""
Backward / retrodiction noise drift: the forward formula with `J ↔ J†` swapped.
"""
function _build_noise_equations_backward(ops, J, Jdagger, rates, efficiencies)
    return _build_noise_equations_forward(ops, Jdagger, J, rates, efficiencies)
end

# `_avg_extra_term` supplies the backward dY averaged-drift correction used by the
# noise meanfield path.
_noise_builder(::Forward) = _build_noise_equations_forward
_noise_builder(::Backward) = _build_noise_equations_backward

_avg_extra_term(::Forward, _, _, _) = nothing
function _avg_extra_term(::Backward, J, Jdagger, eff_rates)
    return op -> _dY_dS_extra_term(op, J, Jdagger, eff_rates)
end
