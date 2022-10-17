function _new_operator(op::IndexedOperator,h,aon=acts_on(op)) 
    if op.ind.hilb != h
        return IndexedOperator(Transition(h,op.op.name,op.op.i,op.op.j,aon;op.op.metadata),Index(h,op.ind.name,op.ind.rangeN,op.ind.specHilb))
    end
    return IndexedOperator(Transition(h,op.op.name,op.op.i,op.op.j,aon;op.op.metadata),op.ind)
end
_new_operator(nOp::NumberedOperator,h,aon=acts_on(nOp)) = NumberedOperator(Transition(h,nOp.op.name,nOp.op.i,nOp.op.j,aon;nOp.op.metadata),nOp.numb)
function _new_operator(sum::IndexedSingleSum,h,aon) 
    newSumIndex = sum.sumIndex
    if sum.sumIndex.hilb != h
        newSumIndex = Index(h,sum.sumIndex.name,sum.sumIndex.rangeN,ind.specHilb)
    end
    newSumNonEquals = Index[]
    for ind in sum.nonEqualIndices
        if ind.hilb != h
            push!(newSumNonEquals,Index(h,ind.name,ind.rangeN,ind.specHilb))
        end
    end
    return IndexedSingleSum(_new_operator(sum.term,h),newSumIndex,newSumNonEquals)
end
function _new_operator(sum::IndexedSingleSum,h) 
    newSumIndex = sum.sumIndex
    if sum.sumIndex.hilb != h
        newSumIndex = Index(h,sum.sumIndex.name,sum.sumIndex.rangeN,sum.sumIndex.specHilb)
    end
    newSumNonEquals = Index[]
    for ind in sum.nonEqualIndices
        if ind.hilb != h
            push!(newSumNonEquals,Index(h,ind.name,ind.rangeN,ind.specHilb))
        end
    end
    return IndexedSingleSum(_new_operator(sum.term,h),newSumIndex,newSumNonEquals)
    
end
function _new_indices(Inds::Vector,h)
    Inds_ = copy(Inds)
    for i=1:length(Inds)
        Inds_[i] = Index(h,Inds[i].name,Inds[i].rangeN,Inds[i].specHilb)
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
entities as argument-values.

See also: [`CorrelationFunction`](@ref)
"""
function IndexedCorrelationFunction(op1,op2,de0::AbstractMeanfieldEquations;
    steady_state=false, add_subscript=0,
    filter_func=nothing, mix_choice=maximum,
    iv=SymbolicUtils.Sym{Real}(:τ),
    order=nothing,
    extra_indices::Vector=[:i,:j,:k,:l,:m,:n,:p,:q,:r,:s,:t],
    scaling::Bool=false,
    simplify=true, kwargs...)
    h1 = hilbert(op1)
    h2 = _new_hilbert(hilbert(op2), acts_on(op2))
    h = h1⊗h2

    allInds = getAllIndices(de0)
    filter!(x -> x ∉ getIndName.(allInds),extra_indices)

    extra_inds = copy(extra_indices)
    extras_ = _new_indices(getAllIndices(de0.states),h) #the extra indices used in the equations beforehand

    H0 = de0.hamiltonian
    J0 = de0.jumps
    Jd0 = de0.jumps_dagger

    op1_ = _new_operator(op1, h)
    op2_ = _new_operator(op2, h, length(h.spaces); add_subscript=add_subscript)
    op2_0 = _new_operator(op2, h)
    H = _new_operator(H0, h)
    J = [_new_operator(j, h) for j in J0]
    Jd = [_new_operator(j, h) for j in Jd0]
    lhs_new = [_new_operator(l, h) for l in de0.states]

    order_ = if order===nothing
        if de0.order===nothing
            order_lhs = maximum(get_order(l) for l in de0.states)
            order_corr = get_order(op1_*op2_)
            max(order_lhs, order_corr)
        else
            de0.order
        end
    else
        order
    end

    if length(extras_) < order_ && !isemtpy(extras_)
        if length(extras_) + length(extra_indis) < order_
            error("For higher order, more extra_indices are required!")
        end
        for elem in extra_inds
            if elem isa Symbol
                push!(extras_,Index(extras_[1].hilb,elem,extras_[1].rangeN,extras[1].specHilb))
            elseif elem isa Index
                push!(extras_,elem)
            end
        end
    end
    unique!(extras_)
    if isempty(extras_)
        extras_ = extra_inds
    end

    op_ = op1_*op2_
    @assert get_order(op_) <= order_

    varmap = make_varmap(lhs_new, de0.iv)

    de0_ = begin
        eqs = Symbolics.Equation[]
        eqs_op = Symbolics.Equation[]
        ops = map(undo_average, lhs_new)
        for i=1:length(de0.equations)
            rhs = _new_operator(de0.equations[i].rhs, h)
            rhs_op = _new_operator(de0.operator_equations[i].rhs, h)
            push!(eqs, Symbolics.Equation(lhs_new[i], rhs))
            push!(eqs_op, Symbolics.Equation(ops[i], rhs_op))
        end
        IndexedMeanfieldEquations(eqs,eqs_op,lhs_new,ops,H,J,Jd,de0.rates,de0.iv,varmap,order_)
    end

    de = indexed_meanfield([op_],H,J;Jdagger=Jd,rates=de0.rates,iv=iv,order=order_)
    indexed_complete_corr!(de, length(h.spaces), lhs_new, order_, steady_state, de0_;
            filter_func=filter_func,
            mix_choice=mix_choice,
            simplify=simplify,
            extra_indices=extras_,
            scaling=scaling,
            kwargs...) 
    if scaling
        de = scaleME(de)
        de0_ = scaleME(de0_)
        de = subst_reds(de;scaling=true)
        de0_ = subst_reds(de0_;scaling=true)
    end
    de = substituteIntoCorrelation(de,de0_;scaling=scaling)
    
    return CorrelationFunction(op1_, op2_, op2_0, de0_, de, steady_state)
end

function substituteIntoCorrelation(me,de0;scaling::Bool=false)
    neweqs = []
    for eq in me.equations
        push!(neweqs,Symbolics.Equation(eq.lhs,subst_reds(eq.rhs,de0.states;scaling=scaling)))
    end
    return IndexedMeanfieldEquations(neweqs,me.operator_equations,me.states,me.operators,me.hamiltonian,me.jumps,me.jumps_dagger,me.rates,me.iv,me.varmap,me.order)
end

#the function below is quite similar to the indexed_complete function
function indexed_complete_corr!(de,aon0,lhs_new,order,steady_state,de0;
        mix_choice=maximum,
        simplify::Bool=true,
        filter_func=nothing,
        extra_indices::Vector=[],
        scaling::Bool=false,
        kwargs...)
    vs = de.states
    H = de.hamiltonian
    J = de.jumps
    Jd = de.jumps_dagger
    rates = de.rates

    extras = copy(extra_indices)
    maxNumb = maximum(length.(getIndices.(de0.operators)))

    sort!(extra_indices)

    if isempty(extra_indices)
        error("can not complete equations with empty extra_indices!")
    end

    for ind in extra_indices
        if typeof(ind) != typeof(extra_indices[1])
            error("Cannot use extra_indices of different types. Use either only Symbols or Index-Objects!")
        end
    end

    if containsMultiple(getAllIndices(de.states))
        if extra_indices[1] isa Symbol
            error("It is not possible to complete equations, containing indices, that act on different hilbertspaces using Symbols as
            extra_indices. For this case use specific Indices.")
        end
        #maybe write also a check that checks for the indices being correct/enough
    end

    if de.order > maxNumb && de.order - maxNumb > length(extra_indices)
        error("Too few extra_indices provided! Please make sure that for higher orders of cumulant expansion, 
            you also use the extra_indices argument to provide additional indices for calculation. The Number of
            extra_indices provided should be at least $(de.order - maxNumb).
        ")
    end

    vhash = map(hash, vs)
    vs′ = map(_conj, vs)
    vs′hash = map(hash, vs′)
    filter!(!in(vhash), vs′hash)
    missed = find_missing(de.equations, vhash, vs′hash; get_adjoints=false)
    
    missed = find_missing_sums(missed,de;extra_indices=extra_indices,scaling=scaling)
    missed = findMissingSpecialTerms(missed,de;scaling=scaling)
    
    missed = sortByIndex.(missed)
    isnothing(filter_func) || filter!(filter_func, missed) # User-defined filter

    filter!(x -> filterComplete_corr(x,de.states,de0.states,scaling), missed)

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
        me = indexed_meanfield(ops_,H,J;
            Jdagger=Jd,
            rates=rates,
            simplify=simplify,
            order=order,
            iv=de.iv,
            kwargs...)

        _append!(de, me)

        vhash_ = hash.(me.states)
        vs′hash_ = hash.(_conj.(me.states))
        append!(vhash, vhash_)
        for i=1:length(vhash_)
            vs′hash_[i] ∈ vhash_ || push!(vs′hash, vs′hash_[i])
        end

        missed = find_missing(me.equations, vhash, vs′hash; get_adjoints=false)
        
        missed = find_missing_sums(missed,de;extra_indices=extra_indices,scaling=scaling)
        missed = findMissingSpecialTerms(missed,de;scaling=scaling)
        
        missed = sortByIndex.(missed)
        isnothing(filter_func) || filter!(filter_func, missed) # User-defined filter
    
        filter!(x -> filterComplete_corr(x,de.states,de0.states,scaling), missed)
       
        for i = 1:length(missed)
            minds = getIndices(missed[i])
            newMinds = copy(minds)
            for ind1 in minds
                extras_=filterExtras(ind1,extras)
                for k = 1:length(extras_)
                    if findall(x->isequal(x,ind1),extras_)[1] > k && extras_[k] ∉ newMinds #this might go somewhat easier, maybe delete ind2 out of extras after each replacement somehow
                        missed[i] = change_index(missed[i],ind1,extras_[k])
                        newMinds = getIndices(missed[i])
                        break
                    elseif findall(x->isequal(x,ind1),extras_)[1] <= k
                        break
                    end
                end
            end
        end
        filter!(_filter_aon, missed)
        isnothing(filter_func) || filter!(filter_func, missed) # User-defined Filter
        filter!(x -> filterComplete_corr(x,de.states,de0.states,scaling), missed)
        missed = unique(missed) #no duplicates
        missed = elimRed(missed;scaling=scaling)
    end

    if !isnothing(filter_func)
        # Find missing values that are filtered by the custom filter function,
        # but still occur on the RHS; set those to 0
        missed = find_missing(de.equations, vhash, vs′hash; get_adjoints=false)
        if order != 1
            missed = find_missing_sums(missed,de;extra_indices=extra_indices,checking=false,scaling=false)
            missed = findMissingSpecialTerms(missed,de;scaling=false)
        end
        missed_ = sortByIndex.(missed)
        missed = vcat(missed,missed_)
        filter!(!filter_func, missed)
        missed_adj = map(_adjoint, missed)
        subs = Dict(vcat(missed, missed_adj) .=> 0)
        for i=1:length(de.equations)
            de.equations[i] = substitute(de.equations[i], subs; fold=true)
            de.states[i] = de.equations[i].lhs
        end
    end

    return de
end

filterComplete_corr(x,states1,states2,scaling) = (isNotIn(getOps(x;scaling=scaling),getOps.(states1;scaling=scaling),scaling) && isNotIn(getOps(sortByIndex(_conj(x));scaling=scaling),getOps.(states1;scaling=scaling),scaling) 
    && isNotIn(getOps(x;scaling=scaling),getOps.(states2;scaling=scaling),scaling)&& isNotIn(getOps(sortByIndex(_conj(x));scaling=scaling),getOps.(states2;scaling=scaling),scaling))

CorrelationFunction(op1,op2,de0::IndexedMeanfieldEquations; kwargs...) = IndexedCorrelationFunction(op1,op2,de0;kwargs...)