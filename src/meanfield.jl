"""
    meanfield(ops, H, J=QField[]; Jdagger=adjoint.(J), rates=ones(length(J)),
              efficiencies=nothing, direction=Forward(),
              order=nothing, simplify=true,
              mix_choice=maximum, iv=Symbolics.variable(:t))

Compute equations of motion for the averages of `ops` under Hamiltonian `H`
and collapse operators `J` (with rates `rates`). Returns a `MeanFieldEquations`
unless `efficiencies` is given (then a `NoiseMeanFieldEquations`).
"""
_make_iv() = first(MTK.@independent_variables t)

function meanfield(
        ops::AbstractVector,
        H::QField,
        J::AbstractVector = QField[];
        Jdagger = nothing,
        rates::AbstractVector = ones(length(J)),
        efficiencies = nothing,
        direction::EvolutionDirection = Forward(),
        order = nothing,
        simplify::Bool = true,
        mix_choice::Function = maximum,
        iv::Symbolics.Num = _make_iv(),
    )
    if Jdagger === nothing
        Jdagger = _default_jdagger(J)
    end
    Jn, Jdn = _normalize_jumps(J, Jdagger)
    rn = _normalize_rates(rates, length(Jn))
    if efficiencies === nothing
        return _meanfield_deterministic(
            direction, ops, H, Jn, Jdn, rn,
            order, simplify, mix_choice, iv,
        )
    else
        en = _normalize_rates(efficiencies, length(Jn))
        return _meanfield_noise(
            direction, ops, H, Jn, Jdn, rn, en,
            order, simplify, mix_choice, iv,
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
    trace_term = Symbolics.Num(sum(
        rates[i] * (average(Jdagger[i] * J[i]) - average(J[i] * Jdagger[i]))
        for i in eachindex(J)
    ))
    return commutator(-imH, op)
        + _master_lindblad_backward(op, J, Jdagger, rates)
        + op * trace_term
end

function _meanfield_deterministic(
        direction::EvolutionDirection,
        ops, H, J, Jdagger, rates, order,
        simplify, mix_choice, iv,
    )
    imH = im * H
    ops_qa = QAdd[op * 1 for op in ops]
    op_rhs = Vector{QAdd}(undef, length(ops_qa))
    operator_eqs = Vector{Symbolics.Equation}(undef, length(ops_qa))
    @inbounds for (i, op) in enumerate(ops_qa)
        rhs = SQA.expand_completeness(
            _operator_rhs(direction, op, imH, J, Jdagger, rates),
        )
        op_rhs[i] = rhs
        operator_eqs[i] = op ~ rhs
    end
    states = SymbolicUtils.BasicSymbolic[average(op) for op in ops_qa]
    order_vec = order === nothing ? nothing :
        _normalize_order(order, (; hamiltonian = H))
    avg_eqs = Vector{Symbolics.Equation}(undef, length(ops_qa))
    @inbounds for i in eachindex(op_rhs)
        rhs = average(op_rhs[i])
        if order_vec !== nothing
            rhs = cumulant_expansion(rhs, order_vec; simplify = false, mix_choice)
        end
        simplify && (rhs = SymbolicUtils.simplify(rhs))
        avg_eqs[i] = states[i] ~ rhs
    end
    return MeanFieldEquations(
        avg_eqs, operator_eqs, states, ops_qa,
        H, collect(J), collect(Jdagger), collect(rates),
        iv, order_vec, direction,
    )
end

function _lindblad_rhs(op, J, Jdagger, rates)
    isempty(J) && return zero(op)
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
            acc += (rk / 2) * (
                Jdagger[k] * commutator(op, J[k]) +
                    commutator(Jdagger[k], op) * J[k]
            )
        end
    end
    return acc
end

# Adjoint-action Lindblad recycling for the backward Heisenberg/Kalman
# retrodiction picture. Matches master's `_master_lindblad_backward`.
function _master_lindblad_backward(op, J, Jdagger, rates)
    isempty(J) && return zero(op)
    acc = zero(op)
    @inbounds for k in eachindex(J)
        rk = rates[k]
        rk isa AbstractMatrix && error(
            "Nondiagonal measurements not supported in backward retrodiction",
        )
        acc += (-rk / 2) * op * Jdagger[k] * J[k]
        acc += (-rk / 2) * Jdagger[k] * J[k] * op
        acc += rk * J[k] * op * Jdagger[k]
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

function _meanfield_noise(
        direction::Forward, ops, H, J, Jdagger, rates,
        efficiencies, order, simplify, mix_choice, iv
    )
    imH = im * H
    ops_qa = QAdd[op * 1 for op in ops]
    op_rhs = Vector{QAdd}(undef, length(ops_qa))
    operator_eqs = Vector{Symbolics.Equation}(undef, length(ops_qa))
    @inbounds for (i, op) in enumerate(ops_qa)
        rhs = SQA.expand_completeness(
            _operator_rhs(direction, op, imH, J, Jdagger, rates),
        )
        op_rhs[i] = rhs
        operator_eqs[i] = op ~ rhs
    end
    states = SymbolicUtils.BasicSymbolic[average(op) for op in ops_qa]
    order_vec = order === nothing ? nothing :
        _normalize_order(order, (; hamiltonian = H))
    avg_eqs = Vector{Symbolics.Equation}(undef, length(ops_qa))
    @inbounds for i in eachindex(op_rhs)
        rhs = average(op_rhs[i])
        if order_vec !== nothing
            rhs = cumulant_expansion(rhs, order_vec; simplify = false, mix_choice)
        end
        simplify && (rhs = SymbolicUtils.simplify(rhs))
        avg_eqs[i] = states[i] ~ rhs
    end
    op_noise, avg_noise = _build_noise_equations_forward(
        ops_qa, J, Jdagger, rates, efficiencies, simplify
    )
    if order_vec !== nothing
        avg_noise = [
            eq.lhs ~ cumulant_expansion(
                    eq.rhs, order_vec;
                    simplify = false, mix_choice
                )
                for eq in avg_noise
        ]
    end
    if simplify
        avg_noise = [eq.lhs ~ SymbolicUtils.simplify(eq.rhs) for eq in avg_noise]
    end
    return NoiseMeanFieldEquations(
        avg_eqs, avg_noise, operator_eqs, op_noise,
        states, ops_qa, H,
        collect(J), collect(Jdagger),
        collect(rates), collect(efficiencies),
        iv, order_vec, Forward()
    )
end

function _meanfield_noise(
        direction::Backward, ops, H, J, Jdagger, rates,
        efficiencies, order, simplify, mix_choice, iv
    )
    imH = im * H
    ops_qa = QAdd[op * 1 for op in ops]
    op_rhs = Vector{QAdd}(undef, length(ops_qa))
    operator_eqs = Vector{Symbolics.Equation}(undef, length(ops_qa))
    eff_rates = efficiencies .* rates
    @inbounds for (i, op) in enumerate(ops_qa)
        # Operator equation: time-reversed coherent drift + adjoint Lindblad
        # recycling + trace-preserving term. The `_dY_dS_extra_term` lives in
        # average-space (it contains `⟨J⟩·⟨op⟩` cumulant cross-products), so
        # it's added only to the averaged equation below — matching master's
        # `meanfield_backward` semantics where `rhs[i] = rhs_ + rhs_diss +
        # rhs_trace + rhs_dY_dS` is then `average`d.
        rhs = SQA.expand_completeness(
            _operator_rhs(direction, op, imH, J, Jdagger, rates),
        )
        op_rhs[i] = rhs
        operator_eqs[i] = op ~ rhs
    end
    states = SymbolicUtils.BasicSymbolic[average(op) for op in ops_qa]
    order_vec = order === nothing ? nothing :
        _normalize_order(order, (; hamiltonian = H))
    avg_eqs = Vector{Symbolics.Equation}(undef, length(ops_qa))
    @inbounds for i in eachindex(op_rhs)
        rhs = average(op_rhs[i]) +
            _dY_dS_extra_term(ops_qa[i], J, Jdagger, eff_rates)
        if order_vec !== nothing
            rhs = cumulant_expansion(rhs, order_vec; simplify = false, mix_choice)
        end
        simplify && (rhs = SymbolicUtils.simplify(rhs))
        avg_eqs[i] = states[i] ~ rhs
    end
    op_noise, avg_noise = _build_noise_equations_backward(
        ops_qa, J, Jdagger, rates, efficiencies, simplify
    )
    if order_vec !== nothing
        avg_noise = [
            eq.lhs ~ cumulant_expansion(
                    eq.rhs, order_vec;
                    simplify = false, mix_choice
                )
                for eq in avg_noise
        ]
    end
    if simplify
        avg_noise = [eq.lhs ~ SymbolicUtils.simplify(eq.rhs) for eq in avg_noise]
    end
    return NoiseMeanFieldEquations(
        avg_eqs, avg_noise, operator_eqs, op_noise,
        states, ops_qa, H,
        collect(J), collect(Jdagger),
        collect(rates), collect(efficiencies),
        iv, order_vec, Backward()
    )
end
