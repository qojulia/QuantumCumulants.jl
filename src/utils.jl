"""
    find_missing(me::MeanfieldEquations, vs_adj=nothing, get_adjoints=true)

Find all averages on the right-hand-side of in `me.equations` that are not
listed `me.states`. For a complete system this list is empty.

Optional arguments
=================

*`vs_adj`: List of the complex conjugates of `me.states`. If set to `nothing`
    the list is generated internally.
*`get_adjoints=true`: Specify whether a complex conjugate of an average should be
    explicitly listed as missing.

see also: [`complete`](@ref), [`complete!`](@ref)
"""
function find_missing(me::AbstractMeanfieldEquations; vs_adj = nothing, get_adjoints = true)
    vs = me.states
    vhash = map(hash, vs)
    vs′ = if vs_adj===nothing
        map(_conj, vs)
    else
        vs_adj
    end
    vs′hash = map(hash, vs′)
    filter!(!in(vhash), vs′hash)

    missed = []
    missed_hashes = UInt[]

    eqs = me.equations
    for i = 1:length(eqs)
        find_missing!(
            missed,
            missed_hashes,
            eqs[i].rhs,
            vhash,
            vs′hash;
            get_adjoints = get_adjoints,
        )
    end
    return missed
end

function find_missing!(
    missed,
    missed_hashes,
    r::SymbolicUtils.Symbolic,
    vhash,
    vs′hash;
    get_adjoints = true,
)
    if SymbolicUtils.iscall(r)
        for arg ∈ SymbolicUtils.arguments(r)
            find_missing!(
                missed,
                missed_hashes,
                arg,
                vhash,
                vs′hash;
                get_adjoints = get_adjoints,
            )
        end
    end
    return missed
end
function find_missing!(
    missed,
    missed_hashes,
    r::Average,
    vhash,
    vs′hash;
    get_adjoints = true,
)
    rhash = hash(r)
    if !(rhash ∈ vhash) && !(rhash ∈ vs′hash) && !(rhash ∈ missed_hashes)
        push!(missed, r)
        push!(missed_hashes, rhash)
        if !get_adjoints
            # To avoid collecting adjoints as missing variables,
            # collect the hash of the adjoint right away
            r′ = _conj(r)
            r′hash = hash(r′)
            if !(r′hash ∈ missed_hashes)
                push!(missed_hashes, r′hash)
            end
        end
    end
    return missed
end
find_missing!(missed, missed_hashes, r::Number, vs, vs′hash; kwargs...) = missed

function find_missing(
    eqs::Vector,
    vhash::Vector{UInt},
    vs′hash::Vector{UInt};
    get_adjoints = true,
)
    missed = []
    missed_hashes = UInt[]
    for i = 1:length(eqs)
        find_missing!(
            missed,
            missed_hashes,
            eqs[i].rhs,
            vhash,
            vs′hash;
            get_adjoints = get_adjoints,
        )
    end
    return missed
end

"""
    complete(de::MeanfieldEquations)

From a set of differential equation of averages, find all averages that are missing
and derive the corresponding equations of motion. Uses [`find_missing`](@ref)
and [`meanfield`](@ref) to do so.

Optional arguments
==================

*`order=de.order`: The order at which the [`cumulant_expansion`](@ref) is performed
    on the newly derived equations. If `nothing`, the order is inferred from the
    existing equations.
*`filter_func=nothing`: Custom function that specifies whether some averages should
    be ignored when completing a system. This works by calling `filter!(filter_func, missed)`
    where `missed` is the vector resulting from [`find_missing`](@ref). Occurrences
    of averages for which `filter_func` returns `false` are substituted to 0.
*`extra_indices=Vector`: Used for indexed equations. Can be used to specify additional
    indices, that are needed for calculation.
*`kwargs...`: Further keyword arguments are passed on to [`meanfield`](@ref) and
    simplification.

see also: [`find_missing`](@ref), [`meanfield`](@ref)
"""
function MTK.complete(de::AbstractMeanfieldEquations; kwargs...)
    de_ = deepcopy(de)
    complete!(de_; kwargs...)
    return de_
end

"""
    complete!(de::MeanfieldEquations)

In-place version of [`complete`](@ref)
"""
function MTK.complete!(
    de::AbstractMeanfieldEquations;
    order = de.order,
    multithread = false,
    filter_func = nothing,
    mix_choice = maximum,
    simplify = true,
    kwargs...,
)
    vs = de.states
    order_lhs = maximum(get_order.(vs))
    order_rhs = 0
    for i = 1:length(de.equations)
        k = get_order(de.equations[i].rhs)
        k > order_rhs && (order_rhs = k)
    end
    if order === nothing
        order_ = max(order_lhs, order_rhs)
    else
        order_ = order
    end
    maximum(order_) >= order_lhs || error(
        "Cannot form cumulant expansion of derivative; you may want to use a higher order!",
    )

    if order_ != de.order
        for i = 1:length(de.equations)
            lhs = de.equations[i].lhs
            rhs = cumulant_expansion(
                de.equations[i].rhs,
                order_;
                mix_choice = mix_choice,
                simplify = simplify,
            )
            de.equations[i] = Symbolics.Equation(lhs, rhs)
        end
    end

    vhash = map(hash, vs)
    vs′ = map(_conj, vs)
    vs′hash = map(hash, vs′)
    filter!(!in(vhash), vs′hash)
    missed = find_missing(de.equations, vhash, vs′hash; get_adjoints = false)
    isnothing(filter_func) || filter!(filter_func, missed) # User-defined filter

    while !isempty(missed)
        ops_ = [SymbolicUtils.arguments(m)[1] for m in missed]
        me = meanfield(
            ops_,
            de.hamiltonian,
            de.jumps;
            Jdagger = de.jumps_dagger,
            rates = de.rates,
            simplify = simplify,
            multithread = multithread,
            order = order_,
            mix_choice = mix_choice,
            iv = de.iv,
            kwargs...,
        )

        _append!(de, me)

        vhash_ = hash.(me.states)
        vs′hash_ = hash.(_conj.(me.states))
        append!(vhash, vhash_)
        for i = 1:length(vhash_)
            vs′hash_[i] ∈ vhash_ || push!(vs′hash, vs′hash_[i])
        end

        missed = find_missing(me.equations, vhash, vs′hash; get_adjoints = false)
        isnothing(filter_func) || filter!(filter_func, missed) # User-defined filter
    end

    if !isnothing(filter_func)
        # Find missing values that are filtered by the custom filter function,
        # but still occur on the RHS; set those to 0
        missed = find_missing(de.equations, vhash, vs′hash; get_adjoints = false)
        filter!(!filter_func, missed)
        missed_adj = map(_adjoint, missed)
        subs = Dict(vcat(missed, missed_adj) .=> 0)
        for i = 1:length(de.equations)
            de.equations[i] = substitute(de.equations[i], subs)
            de.states[i] = de.equations[i].lhs
        end
    end
    return de
end

"""
    initial_values(eqs::MeanfieldEquations, state; level_map=nothing)

For a set of symbolic equations `eqs` compute the initial state average values
corresponding to the numeric quantum state `state` of the system. The quantum
state can either be of type `QuantumOpticsBase.StateVector` or `QuantumOpticsBase.Operator`.

See also: [`to_numeric`](@ref), [`numeric_average`](@ref)
"""
function initial_values(de::MeanfieldEquations, state; kwargs...)
    vs = de.states
    vals = eltype(state)[]
    for v ∈ vs
        push!(vals, numeric_average(v, state; kwargs...))
    end
    return vals
end

function initial_values(de::AbstractMeanfieldEquations, state; kwargs...)
    vs = de.states
    vals = eltype(state)[]
    for v ∈ vs
        push!(vals, numeric_average(v, state; kwargs...))
    end
    return vals
end


# Overload getindex to obtain solutions with averages
for T ∈ [:AbstractTimeseriesSolution, :AbstractNoTimeSolution]
    @eval function Base.getindex(sol::SciMLBase.$(T), avg::Average)
        ode_func = sol.prob.f
        t = MTK.get_iv(ode_func.sys)
        vars = SciMLBase.getsyms(sol)
        var = make_var(avg, t)

        if any(isequal(var), vars)
            # success, we found the symbol
            return getindex(sol, var)
        end

        # couldn't find the symbol, so let's assume we have a conjugate here
        var = make_var(_conj(avg), t)
        !any(isequal(var), vars) &&
            throw(ArgumentError("The average $avg isn't part of the system!"))
        return map(conj, getindex(sol, var))
    end
    @eval Base.getindex(sol::SciMLBase.$(T), op::QNumber) = getindex(sol, average(op))
end

# Internal functions

"""
    get_solution(sol, op::QTerm)
    get_solution(sol, op::QNumber)

Returns the result for the average of the operator expression `op` in the solution
`sol` of an ODE- or SteadyStateProblem, similar to `sol[op]`. It can also be used
for linear combinations of operators, which is not possible with `sol[op]`.
"""
get_solution(sol, op::QNumber) = sol[op]
function get_solution(sol, x)
    if length(sol[:, 1]) == 1 #SteadyStateProblem
        return x
    else
        return x*ones(length(sol))
    end
end
function get_solution(sol, op::QTerm)
    f = SymbolicUtils.operation(op)
    args = SymbolicUtils.arguments(op)
    if f===(+)
        sol_args = [get_solution(sol, args[i]) for i = 1:length(args)]
        return f(sol_args...)
    elseif f === (*)
        c = args[1]
        if length(args) > 2
            op_ = f(args[2:end]...)
        else
            op_ = args[2]
        end
        return c*sol[op_]
    end
end
function get_solution(sol, op::SymbolicUtils.BasicSymbolic{CNumber})
    f = SymbolicUtils.operation(op)
    args = SymbolicUtils.arguments(op)
    sol_args = [get_solution(sol, args[i]) for i = 1:length(args)]
    (f).(sol_args...)
end
function get_solution(sol, op::SymbolicUtils.BasicSymbolic{SQA.AvgSym})
    sol[op]
end
function get_solution(sol, op, dict::Dict)
    x = get_solution(sol, op)
    [substitute(x_, dict) for x_ in x]
end
function get_scale_solution(sol, op::Average, eqs; kwargs...)
    if op in eqs.states
        return sol[op]
    else
        ind_ = findfirst(x -> isscaleequal(op, x; kwargs...), eqs.states)
        if !=(ind_, nothing)
            return sol[eqs.states[ind_]]
        end
    end
    return sol[op]
end

# Takes the averages of all functions in the array rhs and simplifies them if the flag is set
function take_function_averages(rhs, simplify)
    if simplify
        rhs = map(SymbolicUtils.simplify, rhs)
    end
    rhs_avg = map(average, rhs)
    return rhs_avg, map(undo_average, rhs_avg)
end
function SQA._conj(v::Average)
    arg = v.arguments[1]
    adj_arg = adjoint(arg)
    if has_cluster(arg)
        aons, N = SQA.get_cluster_stuff(hilbert(arg))
        names = get_names(arg)
        return substitute_redundants(_average(adj_arg), aons, names)
    else
        return _average(adj_arg)
    end
end

"""
    modify_equations(eqs::MeanfieldEquations, f(lhs,rhs)::Function)

Modify the symbolic quations `eqs` with the function `f`. The function 
`f(lhs,rhs)` needs to return the desired new rhs of the equation. 
The first argument of `f` represents the lhs of each equation and the second 
argument is the rhs. For the lhs of the equation the operator is used in the 
function, not the average. All equations in `eqs.equations` are affected by `f`. 

For example, to add terms from an additional Hamiltonian `Hadd` with a second order
cumulant expansion, we have:

```
function f(lhs, rhs) 
    term = cumulant_expansion(average(commutator(1im*Hadd, lhs)), 2)
    return rhs + term
end
```

see also: [`modify_equations!`](@ref)
"""
function modify_equations(eqs::AbstractMeanfieldEquations, f::Function) 
    eqs_ = deepcopy(eqs)
    modify_equations!(eqs_, f)
    return eqs_
end

"""
    modify_equations!(eqs::MeanfieldEquations)

In-place version of [`modify_equations`](@ref)
"""
function modify_equations!(eqs::AbstractMeanfieldEquations, f::Function)
    eqs_vec = eqs.equations
    for i = 1:length(eqs_vec)
        lhs = eqs_vec[i].lhs
        rhs = eqs_vec[i].rhs
        lhs_op = undo_average(lhs)
        eqs_vec[i] = lhs ~ f(lhs_op, rhs)
    end
end
# TODO: modify_stochastic_equations() (proportional to dW); tests
