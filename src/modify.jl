"""
    modify_equations(eqs::AbstractMeanFieldEquations, f::Function)

Return a copy of `eqs` whose RHS for each equation has been rewritten by `f`.
The function is called as `f(lhs_op, rhs)`, receiving the *operator* form of
the LHS (i.e. `undo_average(eq.lhs)`) and the symbolic RHS expression, and
returning a new RHS expression. Useful for adding measurement-record
correction terms or other deterministic adjustments to a pre-derived set of
equations.

```julia
function f(lhs, rhs)
    term = cumulant_expansion(average(commutator(1im*Hadd, lhs)), 2)
    return rhs + term
end
eqs_mod = modify_equations(eqs, f)
```

See also [`modify_equations!`](@ref).
"""
modify_equations(eqs::AbstractMeanFieldEquations, f::Function) =
    modify_equations!(_copy(eqs), f)

"""
    modify_equations!(eqs::AbstractMeanFieldEquations, f::Function)

In-place version of [`modify_equations`](@ref). Walks `eqs.equations` and
replaces each RHS with `f(undo_average(lhs), rhs)`.
"""
function modify_equations!(eqs::AbstractMeanFieldEquations, f::Function)
    for i in eachindex(eqs.equations)
        lhs = eqs.equations[i].lhs
        rhs = eqs.equations[i].rhs
        eqs.equations[i] = lhs ~ f(SQA.undo_average(lhs), rhs)
    end
    return eqs
end
