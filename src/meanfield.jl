_make_iv() = MTK.t_nounits

"""
    meanfield(ops, H, J=QField[]; Jdagger=adjoint.(J), rates=ones(length(J)),
              efficiencies=nothing, direction=Forward(),
              order=nothing, mix_choice=maximum, iv=ModelingToolkitBase.t_nounits)

Compute equations of motion for the averages of `ops` under Hamiltonian `H`
and collapse operators `J` (with rates `rates`). Returns a `MeanFieldEquations`
unless `efficiencies` is given (then a `NoiseMeanFieldEquations`).

The returned RHS is left in its raw, unsimplified form. Apply
`SymbolicUtils.simplify` (e.g. via `Symbolics.simplify(eq.rhs; expand=true)`)
to the equations you want to inspect.
"""
function meanfield(
        ops::AbstractVector,
        H::QField,
        J::AbstractVector = QField[];
        Jdagger = nothing,
        rates::Union{AbstractVector, AbstractMatrix} = ones(length(J)),
        efficiencies = nothing,
        direction::EvolutionDirection = Forward(),
        order = nothing,
        mix_choice::Function = maximum,
        iv::Symbolics.Num = _make_iv(),
    )
    # Collective dissipation: wrap a flat `(J, rates::Matrix)` into the
    # nested-cluster layout `_lindblad_rhs` expects.
    if rates isa AbstractMatrix
        size(rates, 1) == size(rates, 2) == length(J) || throw(ArgumentError(
            "rates::Matrix must be square with side length(J) for collective dissipation",
        ))
        J = [collect(J)]
        Jdagger = Jdagger === nothing ? nothing : [collect(Jdagger)]
        rates = [rates]
    end
    if Jdagger === nothing
        Jdagger = _default_jdagger(J)
    end
    Jn, Jdn = _normalize_jumps(J, Jdagger)
    rn = _normalize_rates(rates, length(Jn))
    if efficiencies === nothing
        return _meanfield_deterministic(
            direction, ops, H, Jn, Jdn, rn,
            order, mix_choice, iv,
        )
    else
        en = _normalize_rates(efficiencies, length(Jn))
        return _meanfield_noise(
            direction, ops, H, Jn, Jdn, rn, en,
            order, mix_choice, iv,
        )
    end
end

function _normalize_jumps(J, Jdagger)
    if isempty(J)
        return QField[], QField[]
    end
    return collect(J), collect(Jdagger)
end

# Build a default `Jdagger`. For the standard flat form `J::Vector{QField}`
# this is `adjoint.(J)`; for collective decay `J::Vector{Vector{QField}}` we
# adjoint each inner mode-vector entry-wise (the outer-level `adjoint.` would
# transpose the inner vector, which is the wrong layout for `rates::Matrix`).
function _default_jdagger(J)
    isempty(J) && return QField[]
    if eltype(J) <: AbstractVector
        return [adjoint.(jk) for jk in J]
    end
    return adjoint.(J)
end

function _normalize_rates(rates, n::Int)
    if isempty(rates) && n == 0
        return Symbolics.Num[]
    end
    return collect(rates)
end

meanfield(op::QField, H::QField, args...; kw...) = meanfield([op], H, args...; kw...)

# Build the operator-equation RHS for either direction. Forward uses `+i[H,·]`
# and the standard Lindblad recycling; Backward uses `−i[H,·]` together with
# the adjoint Lindblad recycling (`J ↔ J†`) and the trace-preserving term.
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
    + _master_lindblad_backward(op, J, Jdagger, rates)
    + op * trace_term
end

# Operator-level drift derivation shared by every derivation path: build a
# `op ~ rhs` equation for each requested observable. The direction selects
# the coherent + Lindblad recycling formula via `_operator_rhs`.
function _build_op_drift(direction, ops_qa, H, J, Jdagger, rates)
    imH = im * H
    op_rhs = Vector{QAdd}(undef, length(ops_qa))
    operator_eqs = Vector{Symbolics.Equation}(undef, length(ops_qa))
    @inbounds for (i, op) in enumerate(ops_qa)
        rhs = SQA.expand_completeness(
            _operator_rhs(direction, op, imH, J, Jdagger, rates),
        )
        op_rhs[i] = rhs
        operator_eqs[i] = op ~ rhs
    end
    return op_rhs, operator_eqs
end

# Lift operator drift to the averaged level, optionally adding a per-op
# correction (e.g. `_dY_dS_extra_term` for backward noise) and applying
# cumulant truncation. The output is left unsimplified; callers apply
# `SymbolicUtils.simplify` themselves if they want a canonical form.
function _build_avg_drift_eqs(
        op_rhs, states, ops_qa, extra_term, order_vec, mix_choice,
    )
    avg_eqs = Vector{Symbolics.Equation}(undef, length(op_rhs))
    @inbounds for i in eachindex(op_rhs)
        rhs = average(op_rhs[i])
        extra_term === nothing || (rhs = rhs + extra_term(ops_qa[i]))
        if order_vec !== nothing
            rhs = cumulant_expansion(rhs, order_vec; mix_choice)
        end
        avg_eqs[i] = states[i] ~ rhs
    end
    return avg_eqs
end

# Apply cumulant truncation to pre-built noise drifts; no symbolic simplify.
function _finalize_noise_eqs(avg_noise, order_vec, mix_choice)
    if order_vec !== nothing
        avg_noise = [
            eq.lhs ~ cumulant_expansion(eq.rhs, order_vec; mix_choice)
                for eq in avg_noise
        ]
    end
    return avg_noise
end

function _meanfield_deterministic(
        direction::EvolutionDirection,
        ops, H, J, Jdagger, rates, order,
        mix_choice, iv,
    )
    ops_qa = QAdd[op * 1 for op in ops]
    op_rhs, operator_eqs = _build_op_drift(direction, ops_qa, H, J, Jdagger, rates)
    states = SymbolicUtils.BasicSymbolic[average(op) for op in ops_qa]
    order_vec = order === nothing ? nothing :
        _normalize_order(order, (; hamiltonian = H))
    avg_eqs = _build_avg_drift_eqs(
        op_rhs, states, ops_qa, nothing, order_vec, mix_choice,
    )
    return MeanFieldEquations(
        avg_eqs, operator_eqs, states, ops_qa,
        H, collect(J), collect(Jdagger), collect(rates),
        iv, order_vec, direction,
    )
end

function _lindblad_rhs(op, J, Jdagger, rates)
    isempty(J) && return zero(op)
    op_idx = _op_free_indices(op)
    acc = zero(op)
    @inbounds for k in eachindex(J)
        rk = rates[k]
        if rk isa AbstractMatrix
            # Collective decay: `J[k]` and `Jdagger[k]` are vectors of mode
            # operators; `rk[i,j]` is the cross-rate between mode `i` and `j`.
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

# Free indices of an operator that are not bound by a sum scope. An indexed
# jump σ_i^{21} carries a free `Index` i; wrapping the per-jump dissipator in
# Σ_i is what produces independent-decay (rather than collective) semantics
# and lets SQA's diagonal split fire (same-site i = j contributions reduce
# algebraically, different-site i ≠ j contributions vanish via commutation).
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

# Wrap `term` in Σ over each free index that originates in the jump and is
# not already an LHS observable index. Returns `term` unchanged when there
# are no jump-specific free indices (non-indexed jumps, or user-supplied
# explicitly Σ-wrapped jumps).
function _sum_over_jump_indices(term, jump, op_idx)
    jump_idx = _op_free_indices(jump)
    isempty(jump_idx) && return term
    free = SQA.Index[i for i in jump_idx if !(i in op_idx)]
    isempty(free) && return term
    return SQA.Σ(term, free[1], free[2:end]...)
end

# Adjoint-action Lindblad recycling for the backward Heisenberg/Kalman
# retrodiction picture. Matches master's `_master_lindblad_backward`.
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

# Extra deterministic term that arises when the SDE is recast in terms of
# the measurement record `dY` instead of the Wiener increment `dW`. Matches
# master's `_dY_dS_extra_term`.
function _dY_dS_extra_term(op, J, Jdagger, rates)
    out = 0
    @inbounds for k in eachindex(J)
        iszero(rates[k]) && continue
        rk = rates[k]
        rk isa AbstractMatrix && error(
            "Nondiagonal measurements not supported in meanfield_backward",
        )
        # Form `average(J·op) − ⟨J⟩·⟨op⟩` in average-space (subtracting QAdds
        # and BasicSymbolics directly is not supported in v1).
        c1 = rk * (average(J[k] * op) - average(J[k]) * average(op))
        c2 = rk * (average(op * Jdagger[k]) - average(op) * average(Jdagger[k]))
        out = out + (-(c1 + c2)) * average(Jdagger[k] + J[k])
    end
    return out
end

_noise_builder(::Forward) = _build_noise_equations_forward
_noise_builder(::Backward) = _build_noise_equations_backward

# Backward + noise picks up an extra deterministic correction in the dY-form
# of the SDE (master `_master_noise_dY` / `_dY_dS_extra_term`). Forward noise
# and pure-deterministic paths have no such correction.
_avg_extra_term(::Forward, _, _, _) = nothing
function _avg_extra_term(::Backward, J, Jdagger, eff_rates)
    return op -> _dY_dS_extra_term(op, J, Jdagger, eff_rates)
end

function _meanfield_noise(
        direction::EvolutionDirection, ops, H, J, Jdagger, rates,
        efficiencies, order, mix_choice, iv,
    )
    ops_qa = QAdd[op * 1 for op in ops]
    op_rhs, operator_eqs = _build_op_drift(direction, ops_qa, H, J, Jdagger, rates)
    states = SymbolicUtils.BasicSymbolic[average(op) for op in ops_qa]
    order_vec = order === nothing ? nothing :
        _normalize_order(order, (; hamiltonian = H))

    # The dY-form average correction only fires for `Backward` + noise;
    # `_avg_extra_term` returns `nothing` for `Forward`, so the averaged
    # drift comes out identical to the deterministic forward path.
    extra = _avg_extra_term(direction, J, Jdagger, efficiencies .* rates)
    avg_eqs = _build_avg_drift_eqs(
        op_rhs, states, ops_qa, extra, order_vec, mix_choice,
    )

    op_noise, avg_noise = _noise_builder(direction)(
        ops_qa, J, Jdagger, rates, efficiencies,
    )
    avg_noise = _finalize_noise_eqs(avg_noise, order_vec, mix_choice)
    return NoiseMeanFieldEquations(
        avg_eqs, avg_noise, operator_eqs, op_noise,
        states, ops_qa, H,
        collect(J), collect(Jdagger),
        collect(rates), collect(efficiencies),
        iv, order_vec, direction,
    )
end

"""
    simplify!(eqs::AbstractMeanFieldEquations; kwargs...)
    simplify(eqs::AbstractMeanFieldEquations; kwargs...)

Run `SymbolicUtils.simplify` on every RHS in `eqs` (and on the noise drift
RHSs of a `NoiseMeanFieldEquations`). `simplify!` mutates in place; `simplify`
returns a fresh struct. Extra `kwargs` are forwarded to `SymbolicUtils.simplify`
(e.g. `expand=true`).

Use this to opt in to canonical-form RHS for inspection or LaTeX rendering;
the derivation pipeline (`meanfield`, `complete!`, `CorrelationFunction`,
`evaluate`, …) leaves expressions in their raw form so the heavy cost of
`SymbolicUtils.simplify` is only paid when explicitly requested.
"""
function simplify! end

function simplify!(eqs::MeanFieldEquations; kwargs...)
    for (i, eq) in enumerate(eqs.equations)
        eqs.equations[i] = eq.lhs ~ SymbolicUtils.simplify(eq.rhs; kwargs...)
    end
    return eqs
end

function simplify!(eqs::NoiseMeanFieldEquations; kwargs...)
    for (i, eq) in enumerate(eqs.equations)
        eqs.equations[i] = eq.lhs ~ SymbolicUtils.simplify(eq.rhs; kwargs...)
    end
    for (i, eq) in enumerate(eqs.noise_equations)
        eqs.noise_equations[i] = eq.lhs ~ SymbolicUtils.simplify(eq.rhs; kwargs...)
    end
    return eqs
end

SymbolicUtils.simplify(eqs::AbstractMeanFieldEquations; kwargs...) =
    simplify!(_copy(eqs); kwargs...)
