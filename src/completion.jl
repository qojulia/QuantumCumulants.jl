"""
    complete!(eqs::AbstractMeanfieldEquations; max_iter=100_000, filter_func=nothing,
              get_adjoints=true)

Close `eqs` in place by repeatedly deriving equations of motion for every average that
appears on a right-hand side but is not yet a state, until the set is self-contained.
Closure reads and rewrites `eqs.graph`, so user RHS modifications and filter substitutions
already recorded there survive. Right-hand sides are left unsimplified.

# Keyword arguments
* `max_iter=100_000`: runaway backstop; an error is raised if closure has not converged
  within this many iterations.
* `filter_func=nothing`: a predicate `filter_func(avg)`; rejected averages are dropped
  (substituted by 0), e.g. to discard phase-invariant terms.
* `get_adjoints=true`: when `false`, track one representative per conjugate pair, with the
  partner recovered by `conj` at code generation.

See also: [`complete`](@ref), [`find_missing`](@ref), [`meanfield`](@ref).
"""
function complete!(
        eqs::AbstractMeanfieldEquations; max_iter::Int = 100_000,
        filter_func = nothing, get_adjoints::Bool = true,
    )
    keep = filter_func === nothing ? _alltrue : filter_func
    g = closure(eqs.graph; filter = keep, get_adjoints, max_iter)
    filter_func === nothing || (g = map_drifts(g, (_, d) -> _filter_expr(d, filter_func)))
    eqs.graph = g
    resync!(eqs)
    return eqs
end

"""
    complete(eqs::AbstractMeanfieldEquations; order=nothing, kw...)

Non-mutating [`complete!`](@ref): derive equations of motion for every average that
appears on a right-hand side but is not yet a state, returning a new closed system.

# Keyword arguments
* `order=nothing`: if given, the copy is first cumulant-expanded to this order
  (deterministic systems only) before closing.
* `kw...`: forwarded to [`complete!`](@ref) (`max_iter`, `filter_func`, `get_adjoints`).

See also: [`complete!`](@ref), [`find_missing`](@ref), [`meanfield`](@ref).
"""
function MTK.complete(eqs::MeanfieldEquations; order = nothing, kw...)
    base = order === nothing ? _copy(eqs) : cumulant_expansion(eqs, order)
    return complete!(base; kw...)
end

function MTK.complete(eqs::NoiseMeanfieldEquations; order = nothing, kw...)
    order === nothing || throw(
        ArgumentError(
            "`complete(::NoiseMeanfieldEquations; order=...)` is not supported; pass `order` to `meanfield` instead.",
        ),
    )
    return complete!(_copy(eqs); kw...)
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
    g = eqs.graph
    ctx = g.ctx
    treatments = _treatments(eqs, ctx)
    seen = Set{QAdd}()
    for k in keys(g.nodes)
        push!(seen, canonical_rep(k, ctx; treatments)[1])
    end
    missing_states = SymbolicUtils.BasicSymbolic[]
    seen_missing = Set{QAdd}()
    for (_, nd) in g.nodes
        for leaf in _drift_leaves(nd)   # includes noise leaves for a noise system
            op = undo_average(leaf)
            op isa QAdd || continue
            rep, _ = canonical_rep(op, ctx; treatments)
            (rep in seen || rep in seen_missing) && continue
            push!(seen_missing, rep)
            push!(missing_states, average(rep))
        end
    end
    filter_func !== nothing && filter!(filter_func, missing_states)
    return missing_states
end

"""
    closure_report(eqs::AbstractMeanfieldEquations; filter_func=nothing)

Summarise the closure status of `eqs`:
* `n_states`: number of tracked moments.
* `by_order`: tracked-moment count per cumulant order.
* `missing`: moments still needed to close the system ([`find_missing`](@ref)).
* `closed`: whether nothing is missing.

`filter_func` mirrors [`find_missing`](@ref) / [`complete!`](@ref).
"""
function closure_report(eqs::AbstractMeanfieldEquations; filter_func = nothing)
    missing_moments = find_missing(eqs; filter_func)
    orders = get_order.(eqs.operators)
    by_order = OrderedCollections.OrderedDict{Int, Int}(
        ord => count(==(ord), orders) for ord in sort!(unique(orders))
    )
    return (;
        n_states = length(eqs.states),
        by_order,
        missing = missing_moments,
        closed = isempty(missing_moments),
    )
end

# ---- RHS filtering -----------------------------------------------------------

"""
Walk a c-number expression tree, replacing each average leaf `filter_func` rejects with
`0` and rebuilding the enclosing calls; non-average subexpressions pass through.
"""
function _filter_expr(x, filter_func)
    return rewrite(x; descend = y -> SymbolicUtils.iscall(y) && _has_average(y)) do y
        SQA.is_average(y) ? (filter_func(y) ? y : 0) : nothing
    end
end

# ---- struct copy -------------------------------------------------------------

"""
Shallow copy: fresh `graph` (copied `nodes`/`treatments`, shared `sys`/`ctx`) and a fresh
view, so an in-place `!` pass cannot reach the caller's object. Backs the non-mutating
wrappers (`complete`, `simplify`, `modify_equations`).
"""
_copy(eqs::AbstractMeanfieldEquations) = assemble_equations(copy_graph(eqs.graph))
