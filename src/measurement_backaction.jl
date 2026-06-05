# Measurement backaction (Layer 5). Recast a homodyne SDE from the Wiener
# increment `dW` to the measurement record `dY`. The substitution
# `dW = dY - sqrt(2η)·⟨J + J†⟩·dt` leaves the noise drift untouched and adds a
# deterministic correction term to the drift. The correction is exactly
# `-_dY_dS_extra_term` (the same average-space algebra the backward recycling
# uses), cumulant-expanded to the equation set's order.

"""
    translate_W_to_Y(eqs::NoiseMeanFieldEquations; mix_choice=maximum)

Rewrite an SDE whose noise drift is parametrised by the underlying Wiener
process `dW` into one parametrised by the measurement record `dY` instead.
The substitution `dW = dY - sqrt(2η)·⟨J + J†⟩·dt` adds a deterministic
correction to the drift; the noise drift itself is unchanged.

For each equation the RHS is augmented by
`cumulant_expansion(-_dY_dS_extra_term(lhs_op, J, Jdagger, rates .* efficiencies))`,
cumulant-expanded to `eqs.order` when set. Returns a fresh
`NoiseMeanFieldEquations` of the same direction. The augmented RHS is left
unsimplified; apply `simplify` yourself for a canonical form.
"""
function translate_W_to_Y(
        eqs::NoiseMeanFieldEquations;
        mix_choice::Function = maximum,
    )
    out = _copy(eqs)
    J, Jd = out.jumps, out.jumps_dagger
    rates_eff = out.rates .* out.efficiencies
    for i in eachindex(out.equations)
        eq_i = out.equations[i]
        lhs_op = undo_average(eq_i.lhs)
        term = -_dY_dS_extra_term(lhs_op, J, Jd, rates_eff)
        out.order !== nothing && (term = cumulant_expansion(term, out.order; mix_choice))
        out.equations[i] = eq_i.lhs ~ eq_i.rhs + term
    end
    return out
end
