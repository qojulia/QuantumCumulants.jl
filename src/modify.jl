# Equation post-processing: apply a user hook to every RHS of a pre-derived
# equation set without re-deriving the whole system.

"""
    modify_equations(eqs::AbstractMeanFieldEquations, f::Function)

Return a copy of `eqs` whose RHS for each equation has been rewritten by `f`.
The function is called as `f(lhs_op, rhs)`, receiving the *operator* form of
the LHS (`undo_average(eq.lhs)`) and the symbolic RHS, and returning a new RHS.

```julia
f(lhs, rhs) = rhs + cumulant_expansion(average(commutator(1im * Hadd, lhs)), 2)
eqs_mod = modify_equations(eqs, f)
```

See also [`modify_equations!`](@ref).
"""
modify_equations(eqs::AbstractMeanFieldEquations, f::Function) =
    modify_equations!(_copy(eqs), f)

"""
    modify_equations!(eqs::AbstractMeanFieldEquations, f::Function)

In-place version of [`modify_equations`](@ref). Walks `eqs.equations` and replaces each
RHS with `f(undo_average(lhs), rhs)`.
"""
function modify_equations!(eqs::AbstractMeanFieldEquations, f::Function)
    for i in eachindex(eqs.equations)
        lhs = eqs.equations[i].lhs
        rhs = eqs.equations[i].rhs
        eqs.equations[i] = lhs ~ f(undo_average(lhs), rhs)
    end
    return eqs
end
