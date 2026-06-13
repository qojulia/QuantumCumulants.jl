"""
    complete!(eqs::AbstractMeanfieldEquations; max_iter=100_000, filter_func=nothing,
              mix_choice=maximum, get_adjoints=true)

Close `eqs` in place by repeatedly deriving equations of motion for every average that
appears on a right-hand side but is not yet a state, until the set is self-contained.
Right-hand sides are left unsimplified.

# Keyword arguments
* `max_iter=100_000`: runaway backstop; an error is raised if closure has not converged
  within this many iterations.
* `filter_func=nothing`: a predicate `filter_func(avg)`; rejected averages are dropped
  (substituted by 0), e.g. to discard phase-invariant terms.
* `mix_choice=maximum`: passed to [`cumulant_expansion`](@ref) for mixed-order truncation.
* `get_adjoints=true`: when `false`, track one representative per conjugate pair, with the
  partner recovered by `conj` at code generation.

See also: [`complete`](@ref), [`find_missing`](@ref), [`meanfield`](@ref).
"""
function complete!(
        eqs::AbstractMeanfieldEquations; max_iter::Int = 100_000,
        filter_func = nothing, mix_choice = maximum, get_adjoints::Bool = true,
    )
    g = _graph_from_eqs(eqs; mix_choice)
    keep = filter_func === nothing ? _alltrue : filter_func
    closure!(g; filter = keep, get_adjoints, max_iter)
    closed = assemble_equations(g)
    filter_func !== nothing && _filter_rhs!(closed, filter_func)
    _replace_contents!(eqs, closed)
    return eqs
end

"""
    complete(eqs::AbstractMeanfieldEquations; order=nothing, kw...)

Non-mutating [`complete!`](@ref): derive equations of motion for every average that
appears on a right-hand side but is not yet a state, returning a new closed system.

# Keyword arguments
* `order=nothing`: if given, the copy is first cumulant-expanded to this order
  (deterministic systems only) before closing.
* `kw...`: forwarded to [`complete!`](@ref) (`max_iter`, `filter_func`, `mix_choice`,
  `get_adjoints`).

See also: [`complete!`](@ref), [`find_missing`](@ref), [`meanfield`](@ref).
"""
function MTK.complete(
        eqs::MeanfieldEquations; order = nothing, mix_choice = maximum, kw...
    )
    eqs_copy = order === nothing ? _copy(eqs) :
        cumulant_expansion(_copy(eqs), order; mix_choice)
    return complete!(eqs_copy; mix_choice, kw...)
end

function MTK.complete(
        eqs::NoiseMeanfieldEquations; order = nothing, mix_choice = maximum, kw...
    )
    order === nothing || throw(
        ArgumentError(
            "`complete(::NoiseMeanfieldEquations; order=...)` is not supported; pass `order` to `meanfield` instead.",
        ),
    )
    return complete!(_copy(eqs); mix_choice, kw...)
end

"""
    find_missing(eqs::AbstractMeanfieldEquations; filter_func=nothing)

Return the averages that appear on a right-hand side of `eqs` but are not yet among
its tracked states: the equations still needed to close the system. A moment and
its conjugate count as one: if either is already a state, neither is reported.
`filter_func` mirrors [`complete!`](@ref).

The report is representation-invariant: conjugate pairs are always folded onto one key,
so the answer is the same whether `eqs` was closed with `get_adjoints` true or false.

See also: [`complete`](@ref), [`complete!`](@ref).
"""
function find_missing(
        eqs::AbstractMeanfieldEquations; filter_func = nothing,
    )
    # Key states and RHS leaves with the same `canonical_rep` + recorded treatment the
    # numerical-system resolver uses, so a leaf is "missing" exactly when closure would not track it.
    ctx = build_ctx(eqs)
    treatments = _treatments(eqs, ctx)
    seen = Set{QAdd}()
    for s in eqs.states
        op = undo_average(s)
        op isa QAdd && push!(seen, canonical_rep(op, ctx; treatments)[1])
    end
    eqs_list = eqs.equations
    if eqs isa NoiseMeanfieldEquations
        eqs_list = vcat(eqs.equations, eqs.noise_equations)
    end
    missing_states = SymbolicUtils.BasicSymbolic[]
    seen_missing = Set{QAdd}()
    for eq in eqs_list
        for leaf in eachleaf(eq.rhs)
            op = undo_average(leaf)
            op isa QAdd || continue
            # `canonical_rep` folds ⟨X⟩ and ⟨X†⟩ onto one key, so the leaf's own rep is
            # enough; no separate adjoint lookup needed.
            rep, _ = canonical_rep(op, ctx; treatments)
            (rep in seen || rep in seen_missing) && continue
            push!(seen_missing, rep)
            push!(missing_states, average(rep))
        end
    end
    filter_func !== nothing && filter!(filter_func, missing_states)
    return missing_states
end

# ---- RHS filtering -----------------------------------------------------------

function _filter_eq_vec!(v, filter_func)
    for (i, eq) in enumerate(v)
        v[i] = eq.lhs ~ _filter_expr(eq.rhs, filter_func)
    end
    return
end

"""
Zero every average leaf the user's `filter_func` rejects on each equation RHS of `eqs`
(and on the noise-equation RHSs for a noise system), in place.
"""
function _filter_rhs!(eqs::MeanfieldEquations, filter_func)
    _filter_eq_vec!(eqs.equations, filter_func)
    return eqs
end

function _filter_rhs!(eqs::NoiseMeanfieldEquations, filter_func)
    _filter_eq_vec!(eqs.equations, filter_func)
    _filter_eq_vec!(eqs.noise_equations, filter_func)
    return eqs
end

"""
Walk a c-number expression tree, replacing each average leaf `filter_func` rejects with
`0` and rebuilding the enclosing calls; non-average subexpressions pass through.
"""
function _filter_expr(x, filter_func)
    return rewrite(x; descend = y -> SymbolicUtils.iscall(y) && _has_average(y)) do y
        SQA.is_average(y) ? (filter_func(y) ? y : 0) : nothing
    end
end

# ---- struct copy / in-place content swap -------------------------------------

"""
Graft `src`'s contents into `eqs` in place: empty and refill each mutable `Vector` field
(`equations`, `states`, …) from `src`. Lets `complete!` keep its in-place contract while
the closed system is built functionally via the graph. Only the fields completion changes
are swapped; `hamiltonian`, `jumps`, `treatments` and the other shared scalars stay put.
"""
function _replace_contents!(eqs::MeanfieldEquations, src::MeanfieldEquations)
    for (dst, s) in (
            (eqs.equations, src.equations),
            (eqs.operator_equations, src.operator_equations),
            (eqs.states, src.states),
            (eqs.operators, src.operators),
        )
        empty!(dst)
        append!(dst, s)
    end
    return eqs
end

function _replace_contents!(eqs::NoiseMeanfieldEquations, src::NoiseMeanfieldEquations)
    for (dst, s) in (
            (eqs.equations, src.equations),
            (eqs.noise_equations, src.noise_equations),
            (eqs.operator_equations, src.operator_equations),
            (eqs.operator_noise_equations, src.operator_noise_equations),
            (eqs.states, src.states),
            (eqs.operators, src.operators),
        )
        empty!(dst)
        append!(dst, s)
    end
    return eqs
end

"""
Shallow copy of an equation set: fresh `Vector` fields (so an in-place `!` pass cannot
reach the caller's object) over the shared immutable `hamiltonian`, `iv`, `order` and
`direction`. Backs every non-mutating wrapper (`complete`, `simplify`, `modify_equations`,
`evaluate`) that runs a mutating engine on a throwaway. A deep copy would needlessly
clone every operator into a structurally-equal but non-identical object, which can break
SQA logic that compares operators by identity.
"""
function _copy(eqs::MeanfieldEquations)
    return MeanfieldEquations(
        copy(eqs.equations), copy(eqs.operator_equations),
        copy(eqs.states), copy(eqs.operators),
        eqs.hamiltonian, copy(eqs.jumps), copy(eqs.jumps_dagger),
        copy(eqs.rates), eqs.iv, eqs.order, eqs.direction;
        treatments = copy(eqs.treatments),
    )
end

function _copy(eqs::NoiseMeanfieldEquations)
    return NoiseMeanfieldEquations(
        copy(eqs.equations), copy(eqs.noise_equations),
        copy(eqs.operator_equations), copy(eqs.operator_noise_equations),
        copy(eqs.states), copy(eqs.operators),
        eqs.hamiltonian, copy(eqs.jumps), copy(eqs.jumps_dagger),
        copy(eqs.rates), copy(eqs.efficiencies),
        eqs.iv, eqs.order, eqs.direction;
        treatments = copy(eqs.treatments),
    )
end
