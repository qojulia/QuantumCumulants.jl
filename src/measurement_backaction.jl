"""
    translate_W_to_Y(eqs::NoiseMeanFieldEquations; simplify=true, mix_choice=maximum)

Rewrite an SDE whose noise drift is parametrised by the underlying Wiener
process `dW` into one parametrised by the measurement record `dY` instead.
The substitution `dW = dY − √(2η)·⟨J + J†⟩·dt` adds a deterministic
correction to the drift; the noise drift itself is unchanged.

Concretely, for each equation the RHS is augmented by

    cumulant_expansion(-_dY_dS_extra_term(lhs_op, J, Jdagger, rates .* efficiencies))

cumulant-expanded to `eqs.order` (when set). Returns a fresh
`NoiseMeanFieldEquations` of the same direction.

Matches master `src/measurement_backaction.jl::translate_W_to_Y`.
"""
function translate_W_to_Y(
        eqs::NoiseMeanFieldEquations;
        simplify::Bool = true,
        mix_choice::Function = maximum,
    )
    out = _copy(eqs)
    J, Jd = out.jumps, out.jumps_dagger
    rates = out.rates
    eff = out.efficiencies
    rates_eff = rates .* eff
    for i in eachindex(out.equations)
        eq_i = out.equations[i]
        lhs_op = SQA.undo_average(eq_i.lhs)
        # `_dY_dS_extra_term` already lives in average-space and matches
        # master's `_master_noise_dY` algebra (they are the same expression).
        # Master then wraps the result in `average(-...)`; ours is already
        # in the right form, just negate.
        term = -_dY_dS_extra_term(lhs_op, J, Jd, rates_eff)
        if out.order !== nothing
            term = cumulant_expansion(term, out.order; simplify, mix_choice)
        end
        if simplify
            term = SymbolicUtils.simplify(term)
        end
        out.equations[i] = eq_i.lhs ~ eq_i.rhs + term
    end
    return out
end
