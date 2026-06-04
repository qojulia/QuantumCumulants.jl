# Closure (Layer 5). `complete`/`complete!`/`find_missing` are thin wrappers over
# the graph passes: build a graph from the seed equations, run `closure!` (BFS to
# fixpoint) or `frontier` (the missing edge targets), and lower back to the array
# container. All identity/dedup/conjugate logic lives in the kernel.

"""
    complete!(eqs::AbstractMeanFieldEquations; max_iter=100_000, filter_func=nothing,
              mix_choice=maximum, get_adjoints=true)

Close the system in place by deriving equations for every average reachable on a
RHS. `filter_func(avg)` keeps only the averages it accepts (e.g. phase-invariant
ones). `get_adjoints=false` tracks one representative per conjugate pair (the
partner resolved by `conj` at codegen). RHS expressions are left unsimplified.
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

Non-mutating [`complete!`](@ref). If `order` is given, the copy is first
cumulant-expanded to that order (deterministic systems only).
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
    find_missing(eqs::AbstractMeanFieldEquations; filter_func=nothing,
                 get_adjoints=true)

The averages appearing on a RHS whose identity is not yet a state (the graph
frontier). A state and its conjugate count as one: if either is present, neither
is missing.
"""
function find_missing(
        eqs::AbstractMeanFieldEquations; filter_func = nothing,
        get_adjoints::Bool = true, mix_choice = maximum,
    )
    # Coordinate-consistent matching (spec Task 2): fold every tracked state AND
    # every RHS leaf through the SAME `canonical_rep(·; coords)` the codegen
    # resolver uses, in the system's recorded coordinate. The seen-set is built
    # from the FAITHFUL state operators (not a canon_key-collapsed graph), so
    # concrete per-site atoms stay distinct under the Concrete coordinate. A leaf
    # is missing only if neither it nor its conjugate folds to a tracked state.
    ctx = build_ctx(eqs.operators, eqs.hamiltonian, eqs.jumps, eqs.jumps_dagger)
    coords = isempty(eqs.coords) ? all_free_coords(ctx) :
        Dict{Int, Coordinate}(sp => Coordinate(c) for (sp, c) in eqs.coords)
    seen = Set{QAdd}()
    for s in eqs.states
        op = undo_average(s)
        op isa QAdd && push!(seen, canonical_rep(op, ctx; coords)[1])
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
            rep, _ = canonical_rep(op, ctx; coords)
            (rep in seen || rep in seen_missing) && continue
            repc, _ = canonical_rep(adjoint(op), ctx; coords)
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
        coords = copy(eqs.coords),
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
        coords = copy(eqs.coords),
    )
end
