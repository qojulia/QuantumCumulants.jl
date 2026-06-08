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
