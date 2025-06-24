function _new_operator(op::IndexedOperator, h, aon = nothing; kwargs...)
    if !=(aon, nothing)
        if op.ind.hilb != h
            return IndexedOperator(
                _new_operator(op.op, h, aon; kwargs...),
                Index(h, op.ind.name, op.ind.range, aon),
            )
        end
        return IndexedOperator(_new_operator(op.op, h, aon; kwargs...), op.ind)
    end
    if op.ind.hilb != h
        return IndexedOperator(
            _new_operator(op.op, h; kwargs...),
            Index(h, op.ind.name, op.ind.range, op.ind.aon),
        )
    end
    return IndexedOperator(_new_operator(op.op, h; kwargs...), op.ind)
end
function _new_operator(nOp::NumberedOperator, h, aon = nothing; kwargs...)
    if !=(aon, nothing)
        return NumberedOperator(_new_operator(nOp.op, h, aon; kwargs...), nOp.numb)
    end
    return NumberedOperator(_new_operator(nOp.op, h; kwargs...), nOp.numb)
end
function _new_operator(sum::SingleSum, h, aon = nothing; kwargs...)
    newsum_index = sum.sum_index
    if sum.sum_index.hilb != h
        newsum_index = Index(h, sum.sum_index.name, sum.sum_index.range, sum.sum_index.aon)
    end
    newSumNonEquals = Index[]
    for ind in sum.non_equal_indices
        if ind.hilb != h
            push!(newSumNonEquals, Index(h, ind.name, ind.range, ind.aon))
        end
    end
    return SingleSum(
        _new_operator(sum.term, h, aon; kwargs...),
        newsum_index,
        newSumNonEquals,
    )
end
function _new_operator(sum::DoubleSum, h, aon = nothing; kwargs...)
    newsum_index = sum.sum_index
    if sum.sum_index.hilb != h
        newsum_index = Index(h, sum.sum_index.name, sum.sum_index.range, sum.sum_index.aon)
    end
    newSumNonEquals = Index[]
    for ind in sum.NEI
        if ind.hilb != h
            push!(newSumNonEquals, Index(h, ind.name, ind.range, ind.aon))
        end
    end
    inner = _new_operator(sum.innerSum, h, aon; kwargs...)
    return DoubleSum(inner, newsum_index, newSumNonEquals)
end
function _new_operator(sum::IndexedAverageSum, h, aon = nothing; kwargs...)
    newsum_index = sum.sum_index
    if sum.sum_index.hilb != h
        newsum_index = Index(h, sum.sum_index.name, sum.sum_index.range, sum.sum_index.aon)
    end
    newSumNonEquals = Index[]
    for ind in sum.non_equal_indices
        if ind.hilb != h
            push!(newSumNonEquals, Index(h, ind.name, ind.range, ind.aon))
        end
    end
    return IndexedAverageSum(
        _new_operator(sum.term, h, aon; kwargs...),
        newsum_index,
        newSumNonEquals,
    )
end
_new_operator(sym::BasicSymbolic{IndexedAverageSum}, h, aon = nothing; kwargs...) =
    _new_operator(TermInterface.metadata(sym)[IndexedAverageSum], h, aon; kwargs...)

function _new_operator(sum::IndexedAverageDoubleSum, h, aon = nothing; kwargs...)
    newsum_index = sum.sum_index
    if sum.sum_index.hilb != h
        newsum_index = Index(h, sum.sum_index.name, sum.sum_index.range, sum.sum_index.aon)
    end
    newSumNonEquals = Index[]
    for ind in sum.non_equal_indices
        if ind.hilb != h
            push!(newSumNonEquals, Index(h, ind.name, ind.range, ind.aon))
        end
    end
    inner = _new_operator(sum.innerSum, h, aon; kwargs...)
    return IndexedAverageDoubleSum(inner, newsum_index, newSumNonEquals)
end
_new_operator(sym::BasicSymbolic{IndexedAverageDoubleSum}, h, aon = nothing; kwargs...) =
    _new_operator(TermInterface.metadata(sym)[IndexedAverageDoubleSum], h, aon; kwargs...)
_new_operator(sym::BasicSymbolic{IndexedVariable}, h, aon = nothing; kwargs...) =
    _new_operator(TermInterface.metadata(sym)[IndexedVariable], h, aon; kwargs...)
_new_operator(sym::BasicSymbolic{DoubleIndexedVariable}, h, aon = nothing; kwargs...) =
    _new_operator(TermInterface.metadata(sym)[DoubleIndexedVariable], h, aon; kwargs...)
_new_operator(sym::BasicSymbolic{SpecialIndexedAverage}, h, aon = nothing; kwargs...) =
    _new_operator(TermInterface.metadata(sym)[SpecialIndexedAverage], h, aon; kwargs...)
function _new_operator(sym::SpecialIndexedAverage, h, aon = nothing; kwargs...)
    indexMap =
        [(_new_index(f, h, aon), _new_index(l, h, aon)) for (f, l) in sym.indexMapping]
    return SpecialIndexedAverage(_new_operator(sym.term, h, aon; kwargs...), indexMap)
end
function _new_operator(sym::DoubleIndexedVariable, h, aon = nothing; kwargs...)
    if sym.ind1.hilb != h
        newInd1 = Index(h, sym.ind1.name, sym.ind1.range, sym.ind1.aon)
    end
    if sym.ind2.hilb != h
        newInd2 = Index(h, sym.ind2.name, sym.ind2.range, sym.ind2.aon)
    end
    return DoubleIndexedVariable(sym.name, newInd1, newInd2)
end
_new_operator(x::Vector, h, aon = nothing; kwargs...) = _new_operator.(x, h, aon; kwargs...)
function _new_operator(sym::IndexedVariable, h, aon = nothing; kwargs...)
    if sym.ind.hilb != h
        newInd = Index(h, sym.ind.name, sym.ind.range, sym.ind.aon)
    end
    return IndexedVariable(sym.name, newInd)
end

function _new_index(ind::Index, h, aon = nothing)
    if isnothing(aon)
        return Index(h, ind.name, ind.range, ind.aon)
    end
    return Index(h, ind.name, ind.range, aon)
end
function _new_indices(Inds::Vector, h)
    Inds_ = deepcopy(Inds)
    for i = 1:length(Inds)
        Inds_[i] = Index(h, Inds[i].name, Inds[i].range, Inds[i].aon)
    end
    return Inds_
end

"""
    IndexedCorrelationFunction(op1,op2,de0;steady_state=false,add_subscript=0,mix_choice=maximum)

The first-order two-time correlation function of two operators.

The first-order two-time correlation function of `op1` and `op2` evolving under
the system `de0`. The keyword `steady_state` determines whether the original
system `de0` was evolved up to steady state. The arguments `add_subscript`
defines the subscript added to the name of `op2` representing the constant time.

Note that the correlation function is stored in the first index of the underlying
system of equations.

This is the indexed-version of the [`CorrelationFunction`](@ref) and allows [`IndexedOperator`](@ref)
entities as argument-values. This function will automatically be called by [`CorrelationFunction`](@ref),
when the original system `de0` contains any types of [`Index`](@ref) entities.

See also: [`CorrelationFunction`](@ref)
"""
function IndexedCorrelationFunction(
    op1,
    op2,
    de0::AbstractMeanfieldEquations;
    steady_state = false,
    add_subscript = 0,
    filter_func = nothing,
    mix_choice = maximum,
    iv = nothing,
    order = nothing,
    extra_indices::Vector = [:i, :j, :k, :l, :m, :n, :p, :q, :r, :s, :t],
    simplify = true,
    kwargs...,
)
    if iv === nothing
        iv_ = SymbolicUtils.Sym{Real}(:τ)
        iv_ = SymbolicUtils.setmetadata(iv_, Symbolics.VariableSource, (:parameters, :τ))
        iv_ = SymbolicUtils.setmetadata(iv_, MTK.MTKVariableTypeCtx, MTK.PARAMETER)
    else
        iv_ = iv
    end
    h1 = hilbert(op1)
    h2 = _new_hilbert(hilbert(op2), acts_on(op2))
    h = h1⊗h2

    extras_ = _new_indices(get_all_indices(de0.states), h) #the extra indices used in the equations beforehand

    H0 = de0.hamiltonian
    J0 = de0.jumps
    Jd0 = de0.jumps_dagger

    op1_ = _new_operator(op1, h)
    op2_ = _new_operator(op2, h, length(h.spaces); add_subscript = add_subscript)
    op2_0 = _new_operator(op2, h)
    H = _new_operator(H0, h)
    J = [_new_operator(j, h) for j in J0]
    Jd = [_new_operator(j, h) for j in Jd0]
    rates_ = [_new_operator(r, h) for r in de0.rates]
    lhs_new = [_new_operator(l, h) for l in de0.states]

    order_ = if order === nothing
        if de0.order === nothing
            order_lhs = maximum(get_order(l) for l in de0.states)
            order_corr = get_order(op1_*op2_)
            max(order_lhs, order_corr)
        elseif de0.order isa Vector
            ind_ = acts_on(op1)
            vcat(de0.order, [de0.order[ind_]])
        else
            de0.order
        end
    else
        order
    end

    op_ = op1_*op2_

    varmap = make_varmap(lhs_new, de0.iv)

    de0_ = begin
        eqs = Symbolics.Equation[]
        eqs_op = Symbolics.Equation[]
        ops = map(undo_average, lhs_new)
        for i = 1:length(de0.equations)
            rhs = _new_operator(de0.equations[i].rhs, h)
            rhs_op = _new_operator(de0.operator_equations[i].rhs, h)
            push!(eqs, Symbolics.Equation(lhs_new[i], rhs))
            push!(eqs_op, Symbolics.Equation(ops[i], rhs_op))
        end
        IndexedMeanfieldEquations(
            eqs,
            eqs_op,
            lhs_new,
            ops,
            H,
            J,
            Jd,
            rates_,
            de0.iv,
            varmap,
            order_,
        )
    end

    de = indexed_meanfield(
        [op_],
        H,
        J;
        Jdagger = Jd,
        rates = rates_,
        iv = iv_,
        order = order_,
    )
    indexed_complete_corr!(
        de,
        length(h.spaces),
        lhs_new,
        order_,
        steady_state,
        de0_;
        filter_func = filter_func,
        mix_choice = mix_choice,
        simplify = simplify,
        extra_indices = extras_,
        kwargs...,
    )

    return CorrelationFunction(op1_, op2_, op2_0, de0_, de, steady_state)
end
function scale(corr::CorrelationFunction; kwargs...)
    de = scaleME(corr.de; kwargs...)
    de0 = scale(corr.de0; kwargs...)
    de_ = substituteIntoCorrelation(de, de0; scaling = true, kwargs...)
    op1_sc = scale_term(corr.op1)
    op2_sc = scale_term(corr.op2)
    op2_0_sc = scale_term(corr.op2_0)
    return CorrelationFunction(op1_sc, op2_sc, op2_0_sc, de0, de_, corr.steady_state)
end
function evaluate(
    corr::CorrelationFunction;
    op1_aon = 1,
    op2_aon = 2,
    limits = nothing,
    kwargs...,
)
    me = corr.de
    # eqs = [Symbolics.Equation(insert_index(insert_index(eq.lhs,corr.op1.ind,op1_aon),corr.op2.ind,op2_aon),
    #             insert_index(insert_index(eq.rhs,corr.op1.ind,op1_aon),corr.op2.ind,op2_aon))
    #                 for eq in (me.equations)]
    if !isempty(get_indices(corr.op2))
        eqs = [
            Symbolics.Equation(
                insert_index(eq.lhs, corr.op2.ind, op2_aon),
                insert_index(eq.rhs, corr.op2.ind, op2_aon),
            ) for eq in (me.equations)
        ]
    else
        eqs = me.equations
    end

    op1_eval =
        isempty(get_indices(corr.op1)) ? corr.op1 :
        insert_index(corr.op1, corr.op1.ind, op1_aon)
    op2_eval =
        isempty(get_indices(corr.op2)) ? corr.op2 :
        insert_index(corr.op2, corr.op2.ind, op2_aon)
    op2_0_eval =
        isempty(get_indices(corr.op2_0)) ? corr.op2_0 :
        insert_index(corr.op2_0, corr.op2_0.ind, op2_aon)
    if !=(limits, nothing) && limits isa Pair
        limits_ = Dict{SymbolicUtils.BasicSymbolic,Int64}(first(limits) => last(limits));
        limits = limits_
    end
    if limits === nothing
        limits = Dict{SymbolicUtils.BasicSymbolic,Int64}();
    end
    de0 = evaluate(corr.de0; limits = limits, kwargs...)
    states = getLHS.(eqs)
    operats = undo_average.(states)
    varmap = make_varmap(states, me.iv)
    de = IndexedMeanfieldEquations(
        eqs,
        me.operator_equations,
        states,
        operats,
        me.hamiltonian,
        me.jumps,
        me.jumps_dagger,
        me.rates,
        me.iv,
        varmap,
        me.order,
    )
    de_ = evalME(de; limits = limits, kwargs...)
    de_ = substituteIntoCorrelation(de_, de0; scaling = false, kwargs...)
    return CorrelationFunction(op1_eval, op2_eval, op2_0_eval, de0, de_, corr.steady_state)
end
evaluate(corr::CorrelationFunction, op1_ao, op2_ao; kwargs...) =
    evaluate(corr; op1_aon = op1_ao, op2_aon = op2_ao, kwargs...)

#this function is almost similar to the subst_reds function -> maybe merge together?
function substituteIntoCorrelation(me, de0; scaling::Bool = false, kwargs...)
    to_sub = find_missing(me)
    de_states = [me.states; de0.states]
    to_sub = inorder!.(to_sub)
    filter!(x->x ∉ de_states, to_sub)
    to_sub = inorder!.(to_sub)
    to_insert = Vector{Any}(nothing, length(to_sub))
    if scaling
        states = de_states
        counter = 1
        while counter <= length(to_sub)
            elem = to_sub[counter]
            ind_ = findfirst(x -> isscaleequal(elem, x; kwargs...), states)
            if !=(ind_, nothing)
                to_insert[counter] = states[ind_]
                counter = counter + 1
                continue
            end
            ind_ = findfirst(x -> isscaleequal((_inconj(elem)), x; kwargs...), states)
            if !=(ind_, nothing)
                to_insert[counter] = conj(states[ind_])
                counter = counter + 1
                continue
            end
            deleteat!(to_insert, counter)
            deleteat!(to_sub, counter)
        end
    else
        counter = 1
        while counter <= length(to_sub)
            elem = to_sub[counter]
            if _inconj(elem) in de_states
                to_insert[counter] = conj(_inconj(to_sub[counter]))
                counter = counter + 1
            else
                to_insert[counter] = 0
                counter = counter + 1
            end
        end
    end
    subs = Dict(to_sub .=> to_insert)
    eqs = [
        Symbolics.Equation(
            inorder!(substitute(eq.lhs, subs)),
            inorder!(substitute(eq.rhs, subs)),
        ) for eq in me.equations
    ]
    return IndexedMeanfieldEquations(
        eqs,
        me.operator_equations,
        me.states,
        me.operators,
        me.hamiltonian,
        me.jumps,
        me.jumps_dagger,
        me.rates,
        me.iv,
        me.varmap,
        me.order,
    )
end

#the function below is quite similar to the indexed_complete function
function indexed_complete_corr!(
    de,
    aon0,
    lhs_new,
    order,
    steady_state,
    de0;
    mix_choice = maximum,
    simplify::Bool = true,
    filter_func = nothing,
    extra_indices::Vector = [],
    kwargs...,
)

    vs = de.states
    H = de.hamiltonian
    J = de.jumps
    Jd = de.jumps_dagger
    rates = de.rates

    extras = copy(extra_indices)
    maxNumb = maximum(length.(get_indices.(de0.operators)))

    sort!(extra_indices)

    if isempty(extra_indices)
        error("can not complete equations with empty extra_indices!")
    end

    for ind in extra_indices
        if typeof(ind) != typeof(extra_indices[1])
            error(
                "Cannot use extra_indices of different types. Use either only Symbols or Index-Objects!",
            )
        end
    end

    if de.order isa Vector
        order_max = maximum(de.order)
    else
        order_max = de.order
    end

    if order_max > maxNumb && order_max - maxNumb > length(extra_indices)
        error(
            "Too few extra_indices provided! Please make sure that for higher orders of cumulant expansion,
          you also use the extra_indices argument to provide additional indices for calculation. The Number of
          extra_indices provided should be at least $(order_max - maxNumb).
      ",
        )
    end

    vhash = map(hash, vs)
    vs′ = map(_inconj, vs)
    vs′hash = map(hash, vs′)
    filter!(!in(vhash), vs′hash)
    missed = find_missing(de.equations, vhash, vs′hash; get_adjoints = false)
    missed = find_missing_sums(missed, de; extra_indices = extra_indices)
    missed = find_missing_Dsums(missed, de; extra_indices = extra_indices)

    missed = inorder!.(missed)
    isnothing(filter_func) || filter!(filter_func, missed) # User-defined filter

    filter!(
        x -> filterComplete_corr(x, de.states, de0.states, false, steady_state; kwargs...),
        missed,
    )

    missed = inorder!.(missed)

    missed = unique(missed) #no duplicates

    vhash_new = map(hash, lhs_new)
    vhash_new′ = map(hash, _adjoint.(lhs_new))
    filter!(!in(vhash_new), vhash_new′)

    function _filter_aon(x) # Filter values that act only on Hilbert space representing system at time t0
        aon = acts_on(x)
        if aon0 in aon
            length(aon)==1 && return false
            return true
        end
        if steady_state # Include terms without t0-dependence only if the system is not in steady state
            h = hash(x)
            return !(h∈vhash_new || h∈vhash_new′)
        else
            return true
        end
    end
    filter!(_filter_aon, missed)
    isnothing(filter_func) || filter!(filter_func, missed) # User-defined filter

    missed = unique(missed)
    while !isempty(missed)
        ops_ = [SymbolicUtils.arguments(m)[1] for m in missed]
        me = indexed_meanfield(
            ops_,
            H,
            J;
            Jdagger = Jd,
            rates = rates,
            simplify = simplify,
            order = order,
            iv = de.iv,
            kwargs...,
        )

        _append!(de, me)

        vhash_ = hash.(me.states)
        vs′hash_ = hash.(_inconj.(me.states))
        append!(vhash, vhash_)
        for i = 1:length(vhash_)
            vs′hash_[i] ∈ vhash_ || push!(vs′hash, vs′hash_[i])
        end

        missed = find_missing(me.equations, vhash, vs′hash; get_adjoints = false)
        missed = find_missing_sums(missed, de; extra_indices = extra_indices)
        missed = find_missing_Dsums(missed, de; extra_indices = extra_indices)

        missed = inorder!.(missed)
        isnothing(filter_func) || filter!(filter_func, missed) # User-defined filter

        filter!(
            x -> filterComplete_corr(
                x,
                de.states,
                de0.states,
                false,
                steady_state;
                kwargs...,
            ),
            missed,
        )

        missed = inorder!.(missed)

        for i = 1:length(missed)
            minds = get_indices(missed[i])
            newMinds = copy(minds)
            for ind1 in minds
                extras_=filterExtras(ind1, extras)
                for k = 1:length(extras_)
                    if findall(x->isequal(x, ind1), extras_)[1] > k && extras_[k] ∉ newMinds #this might go somewhat easier, maybe delete ind2 out of extras after each replacement somehow
                        missed[i] = change_index(missed[i], ind1, extras_[k])
                        newMinds = get_indices(missed[i])
                        break
                    elseif findall(x->isequal(x, ind1), extras_)[1] <= k
                        break
                    end
                end
            end
        end
        filter!(_filter_aon, missed)
        isnothing(filter_func) || filter!(filter_func, missed) # User-defined Filter
        filter!(
            x -> filterComplete_corr(
                x,
                de.states,
                de0.states,
                false,
                steady_state;
                kwargs...,
            ),
            missed,
        )
        missed = unique(missed) #no duplicates
        missed = elimRed!(missed)
        missed = inorder!.(missed)
    end

    if !isnothing(filter_func)
        # Find missing values that are filtered by the custom filter function,
        # but still occur on the RHS; set those to 0
        missed = find_missing(de.equations, vhash, vs′hash; get_adjoints = false)
        missed = find_missing_sums(
            missed,
            de;
            extra_indices = extra_indices,
            checking = false,
            scaling = false,
        )
        missed = find_missing_Dsums(
            missed,
            de;
            extra_indices = extra_indices,
            checking = false,
            scaling = false,
        )
        missed = find_missing_sums(missed, de; checking = false, scaling = false)
        missed = find_missing_Dsums(missed, de; checking = false, scaling = false)

        missed_ = inorder!.(missed)
        missed = vcat(missed, missed_)
        filter!(!filter_func, missed)
        missed_adj = map(_inconj, missed)
        subs = Dict(vcat(missed, missed_adj) .=> 0)
        for i = 1:length(de.equations)
            de.equations[i] = substitute(de.equations[i], subs; fold = true)
            de.states[i] = de.equations[i].lhs
        end
    end

    return de
end

function filterComplete_corr(x, states1, states2, scaling, steady_state; kwargs...)
    if steady_state
        return (
            isNotIn(x, states1, scaling; kwargs...) &&
            isNotIn(_inconj(x), states1, scaling; kwargs...) &&
            isNotIn(x, states2, scaling; kwargs...) &&
            isNotIn(_inconj(x), states2, scaling; kwargs...)
        )
    else
        return (
            isNotIn(x, states1, scaling; kwargs...) &&
            isNotIn(_inconj(x), states1, scaling; kwargs...)
        )
    end
end

CorrelationFunction(op1, op2, de0::IndexedMeanfieldEquations; kwargs...) =
    IndexedCorrelationFunction(op1, op2, de0; kwargs...)
