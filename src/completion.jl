"""
    complete!(eqs::AbstractMeanFieldEquations; max_iter=100_000, filter_func=nothing,
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
        eqs::MeanFieldEquations; max_iter::Int = 100_000,
        filter_func = nothing, mix_choice = maximum, get_adjoints::Bool = true,
    )
    g = _graph_from_eqs(eqs; mix_choice)
    flt = filter_func === nothing ? _alltrue : avg -> filter_func(avg)
    closure!(g; filter = flt, get_adjoints, max_iter)
    closed = lower_to_eqs(g)
    filter_func !== nothing && _filter_rhs!(closed, filter_func)
    _replace_contents!(eqs, closed)
    return eqs
end

function complete!(
        eqs::NoiseMeanFieldEquations; max_iter::Int = 100_000,
        filter_func = nothing, mix_choice = maximum, get_adjoints::Bool = true,
    )
    g = _graph_from_eqs(eqs; mix_choice)
    flt = filter_func === nothing ? _alltrue : avg -> filter_func(avg)
    closure!(g; filter = flt, get_adjoints, max_iter)
    closed = lower_to_eqs(g)
    filter_func !== nothing && _filter_rhs!(closed, filter_func)
    _replace_contents!(eqs, closed)
    return eqs
end

"""
    complete(eqs::AbstractMeanFieldEquations; order=nothing, kw...)

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
        eqs::MeanFieldEquations; order = nothing, mix_choice = maximum, kw...
    )
    eqs_copy = order === nothing ? _copy(eqs) :
        cumulant_expansion(_copy(eqs), order; mix_choice)
    return complete!(eqs_copy; mix_choice, kw...)
end

function MTK.complete(
        eqs::NoiseMeanFieldEquations; order = nothing, mix_choice = maximum, kw...
    )
    order === nothing || throw(
        ArgumentError(
            "`complete(::NoiseMeanFieldEquations; order=...)` is not supported; pass `order` to `meanfield` instead.",
        ),
    )
    return complete!(_copy(eqs); mix_choice, kw...)
end

"""
    find_missing(eqs::AbstractMeanFieldEquations; filter_func=nothing, get_adjoints=true, mix_choice=maximum)

Return the averages that appear on a right-hand side of `eqs` but are not yet among
its tracked states: the equations still needed to close the system. A moment and
its conjugate count as one: if either is already a state, neither is reported.
`get_adjoints`, `filter_func` and `mix_choice` mirror [`complete!`](@ref).

See also: [`complete`](@ref), [`complete!`](@ref).
"""
function find_missing(
        eqs::AbstractMeanFieldEquations; filter_func = nothing,
        get_adjoints::Bool = true, mix_choice = maximum,
    )
    # Fold every tracked state AND every RHS leaf through the SAME
    # `canonical_rep(·; treatments)` the codegen resolver uses, in the system's
    # recorded treatment. The seen-set is built from the faithful state operators
    # (not a canon_key-collapsed graph), so concrete per-site atoms stay distinct
    # under the Concrete treatment. A leaf is missing only if neither it nor its
    # conjugate folds to a tracked state.
    ctx = build_ctx(eqs.operators, eqs.hamiltonian, eqs.jumps, eqs.jumps_dagger)
    treatments = isempty(eqs.treatments) ? all_free_treatments(ctx) :
        Dict{Int, SubspaceTreatment}(sp => SubspaceTreatment(c) for (sp, c) in eqs.treatments)
    seen = Set{QAdd}()
    for s in eqs.states
        op = undo_average(s)
        op isa QAdd && push!(seen, canonical_rep(op, ctx; treatments)[1])
    end
    eqs_list = eqs.equations
    if eqs isa NoiseMeanFieldEquations
        eqs_list = vcat(eqs.equations, eqs.noise_equations)
    end
    missing_states = SymbolicUtils.BasicSymbolic[]
    seen_missing = Set{QAdd}()
    for eq in eqs_list
        for leaf in eachleaf(eq.rhs)
            op = undo_average(leaf)
            op isa QAdd || continue
            rep, _ = canonical_rep(op, ctx; treatments)
            (rep in seen || rep in seen_missing) && continue
            repc, _ = canonical_rep(adjoint(op), ctx; treatments)
            (repc in seen || repc in seen_missing) && continue
            push!(seen_missing, rep)
            push!(missing_states, average(rep))
            if get_adjoints && !isequal(repc, rep)
                push!(seen_missing, repc)
                push!(missing_states, average(repc))
            end
        end
    end
    filter_func !== nothing && filter!(filter_func, missing_states)
    return missing_states
end

# ---- RHS filtering (zero out leaves the user's filter rejects) ---------------

function _filter_rhs!(eqs::MeanFieldEquations, filter_func)
    for (i, eq) in enumerate(eqs.equations)
        eqs.equations[i] = eq.lhs ~ _filter_expr(eq.rhs, filter_func)
    end
    return eqs
end

function _filter_rhs!(eqs::NoiseMeanFieldEquations, filter_func)
    for (i, eq) in enumerate(eqs.equations)
        eqs.equations[i] = eq.lhs ~ _filter_expr(eq.rhs, filter_func)
    end
    for (i, eq) in enumerate(eqs.noise_equations)
        eqs.noise_equations[i] = eq.lhs ~ _filter_expr(eq.rhs, filter_func)
    end
    return eqs
end

function _filter_expr(x, filter_func)
    if x isa SymbolicUtils.BasicSymbolic
        if SQA.is_average(x)
            return filter_func(x) ? x : 0
        end
        if SymbolicUtils.iscall(x) && _has_average(x)
            op = SymbolicUtils.operation(x)
            args = SymbolicUtils.arguments(x)
            new_args = Any[_filter_expr(a, filter_func) for a in args]
            return op(new_args...)
        end
    end
    return x
end

# ---- struct copy / in-place content swap -------------------------------------

# Swap `eqs`'s array fields in place to match `src` (keeps `complete!`'s in-place
# contract while the closed system is built functionally via the graph).
function _replace_contents!(eqs::MeanFieldEquations, src::MeanFieldEquations)
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

function _replace_contents!(eqs::NoiseMeanFieldEquations, src::NoiseMeanFieldEquations)
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

function _copy(eqs::MeanFieldEquations)
    return MeanFieldEquations(
        copy(eqs.equations), copy(eqs.operator_equations),
        copy(eqs.states), copy(eqs.operators),
        eqs.hamiltonian, copy(eqs.jumps), copy(eqs.jumps_dagger),
        copy(eqs.rates), eqs.iv, eqs.order, eqs.direction;
        initial_operators = copy(eqs.initial_operators),
        treatments = copy(eqs.treatments),
    )
end

function _copy(eqs::NoiseMeanFieldEquations)
    return NoiseMeanFieldEquations(
        copy(eqs.equations), copy(eqs.noise_equations),
        copy(eqs.operator_equations), copy(eqs.operator_noise_equations),
        copy(eqs.states), copy(eqs.operators),
        eqs.hamiltonian, copy(eqs.jumps), copy(eqs.jumps_dagger),
        copy(eqs.rates), copy(eqs.efficiencies),
        eqs.iv, eqs.order, eqs.direction;
        initial_operators = copy(eqs.initial_operators),
        treatments = copy(eqs.treatments),
    )
end
