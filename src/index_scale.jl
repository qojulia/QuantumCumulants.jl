#file for functions regarding scaling
"""
    scaleME(me::IndexedMeanfieldEquations)

Function, that evaluates a given [`MeanfieldEquations`](@ref) entity and returns again equations,
where indices have been inserted and sums evaluated, regarding the same relations, as done when calculating
with oparators using a [`ClusterSpace`](@ref).

# Arguments
*`me::MeanfieldEquations`: A [`MeanfieldEquations`](@ref) entity, which shall be evaluated.

"""
function scaleME(me::IndexedMeanfieldEquations; kwargs...)
    newEqs = []
    for eq in me.equations
        tempEq = scaleEq(eq; kwargs...)
        if tempEq.lhs in getLHS.(newEqs)
            continue
        elseif isNotIn(tempEq.lhs, getLHS.(newEqs), true; kwargs...) &&
               isNotIn(_inconj(tempEq.lhs), getLHS.(newEqs), true; kwargs...)
            push!(newEqs, tempEq)
        end
    end
    vs = getLHS.(newEqs)
    varmap = make_varmap(vs, me.iv)
    ops = undo_average.(vs)
    return IndexedMeanfieldEquations(
        newEqs,
        me.operator_equations,
        vs,
        ops,
        me.hamiltonian,
        me.jumps,
        me.jumps_dagger,
        me.rates,
        me.iv,
        varmap,
        me.order,
    )
end
scale_term(sym::BasicSymbolic{IndexedAverageSum}; kwargs...) =
    scale_term(SymbolicUtils.metadata(sym)[IndexedAverageSum]; kwargs...)
function scale_term(sum::IndexedAverageSum; h = nothing, kwargs...)
    if !=(h, nothing)
        if !(h isa Vector)
            h = [h]
        end
        if !(sum.sum_index.aon in h)
            return IndexedAverageSum(
                scale_term(sum.term; h = h, kwargs...),
                sum.sum_index,
                sum.non_equal_indices,
            )
        end
    end
    NEI = sum.non_equal_indices
    # NEI = filter(x->x in get_indices(sum.term),sum.non_equal_indices)
    prefact = sum.sum_index.range - length(NEI)
    term_ = scale_term(sum.term; h = h, kwargs...)
    return prefact*term_
end
scale_term(sym::BasicSymbolic{IndexedAverageDoubleSum}; kwargs...) =
    scale_term(SymbolicUtils.metadata(sym)[IndexedAverageDoubleSum]; kwargs...)
function scale_term(sum::IndexedAverageDoubleSum; h = nothing, kwargs...)
    if !=(h, nothing)
        if !(h isa Vector)
            h = [h]
        end
        if !(sum.sum_index.aon in h)
            return IndexedAverageDoubleSum(
                scale_term(sum.innerSum; h = h, kwargs...),
                sum.sum_index,
                sum.non_equal_indices,
            )
        end
    end
    NEI = sum.non_equal_indices
    # NEI = filter(x->x in get_indices(sum.innerSum),sum.non_equal_indices)
    prefact = sum.sum_index.range - length(NEI)
    term_ = scale_term(sum.innerSum; h = h, kwargs...)
    return prefact*term_
end
function scale_term(pow::SymbolicUtils.Pow; kwargs...)
    args = arguments(pow)
    return scale_term(args[1]; kwargs...)^(args[2])
end
function scale_term(add::BasicSymbolic{<:CNumber}; h = nothing, kwargs...)
    if iscall(add)
        op = operation(add)
        if op === +
            args = arguments(add)
            adds = Vector{Any}(nothing, length(args))
            for i = 1:length(args)
                adds[i] = scale_term(args[i]; h = h, kwargs...)
            end
            return sum(adds)
        elseif op === *
            mul = add
            mults = []
            for arg in arguments(mul)
                if arg isa BasicSymbolic{DoubleIndexedVariable}
                    # Double numbered variables need extra computation, since they depend on the number of indices within an average
                    # for example Γij * <σi> * <σj> -> Γ11, since both averages only have 1 index
                    # however Γij <σi * σj> -> Γ12 = Γ21 for scaled terms
                    meta = SymbolicUtils.metadata(arg)[DoubleIndexedVariable]
                    inds = []
                    for arg in arguments(mul)
                        push!(inds, get_indices(arg))
                    end
                    if !=(h, nothing)
                        if !(h isa Vector)
                            h = [h]
                        end
                        filter!(x -> x.aon in h, inds)
                    end
                    order = maximum(length.(inds))
                    if order <= 1
                        push!(mults, DoubleNumberedVariable(meta.name, 1, 1))
                    else
                        push!(mults, DoubleNumberedVariable(meta.name, 1, 2))
                    end
                else
                    push!(mults, scale_term(arg; h = h, kwargs...))
                end
            end
            if length(mults) == 1
                return mults[1]
            elseif isempty(mults)
                return 0
            else
                # filter!(x->!=(x,nothing),mults)
                return *(mults...)
            end
        elseif op === ^
            args = arguments(add)
            return scale_term(args[1])^(args[2])
        end
    end
    return add
end

function scale_term(x::Average; h = nothing, kwargs...)
    indices = get_indices(x)
    term = x
    if !=(h, nothing)
        if !(h isa Vector)
            h = [h]
        end
        filter!(x -> x.aon in h, indices)
    end
    dic = get_ind_dic(indices)
    for a in keys(dic)
        for i = 1:length(dic[a])
            term = insert_index(term, dic[a][i], i)
        end
    end
    return term
end
scale_term(x::IndexedOperator; h = nothing, kwargs...) = NumberedOperator(x.op, 1)
function scale_term(x::QMul; h = nothing, kwargs...)
    indices = get_indices(x)
    term = x
    if !=(h, nothing)
        if !(h isa Vector)
            h = [h]
        end
        filter!(x -> x.aon in h, indices)
    end
    dic = get_ind_dic(indices)
    for a in keys(dic)
        for i = 1:length(dic[a])
            term = insert_index(term, dic[a][i], i)
        end
    end
    return term
end
function get_ind_dic(inds)
    dic = Dict()
    aons = unique(SQA.get_aon.(inds))
    for a in aons
        push!(dic, a => filter(x->isequal(a, x.aon), inds))
    end
    return dic
end
scale_term(x::BasicSymbolic{SpecialIndexedAverage}; kwargs...) =
    scale_term(SymbolicUtils.metadata(x)[SpecialIndexedAverage].term; kwargs...) #this is fine, since the intrinsic conditions on the indices go away automatically from the scaling
scaleEq(eq::Symbolics.Equation; kwargs...) =
    Symbolics.Equation(scale_term(eq.lhs; kwargs...), scale_term(eq.rhs; kwargs...))
function scale_term(sym::BasicSymbolic{IndexedVariable}; h = nothing, kwargs...)
    meta = SymbolicUtils.metadata(sym)[IndexedVariable]
    if !=(h, nothing)
        if meta.ind.aon in h
            return SingleNumberedVariable(meta.name, 1)
        else
            return sym
        end
    end
    return SingleNumberedVariable(meta.name, 1)
end
function scale_term(sym::BasicSymbolic{DoubleIndexedVariable}; h = nothing, kwargs...)
    meta = SymbolicUtils.metadata(sym)[DoubleIndexedVariable]
    if !=(h, nothing)
        if meta.ind1.aon in h
            return DoubleNumberedVariable(meta.name, 1, meta.ind2)
        elseif meta.ind2.aon in h
            return DoubleNumberedVariable(meta.name, meta.ind1, 1)
        else
            return sym
        end
    end
    if meta.ind1 != meta.ind2
        return DoubleNumberedVariable(meta.name, 1, 2)
    elseif meta.ind1 == meta.ind2
        return DoubleNumberedVariable(meta.name, 1, 1)
    end
    return sym
end
scale_term(x; kwargs...) = x

SymbolicUtils.substitute(avrg::BasicSymbolic{SpecialIndexedAverage}, subs; fold = false) =
    SymbolicUtils.substitute(
        SymbolicUtils.metadata(avrg)[SpecialIndexedAverage].term,
        subs;
        fold = fold,
    )#SpecialIndexedAverage(SymbolicUtils.substitute(term.metadata.term,subs;fold=fold),term.metadata.indexMapping)
function SymbolicUtils.substitute(sum::BasicSymbolic{IndexedAverageSum}, subs; fold = false)
    meta = SymbolicUtils.metadata(sum)[IndexedAverageSum]
    subTerm = SymbolicUtils.substitute(meta.term, subs; fold = fold)
    if SymbolicUtils._iszero(subTerm)
        return 0
    elseif subTerm isa SQA.symbolics_terms
        return IndexedAverageSum(subTerm, meta.sum_index, meta.non_equal_indices)
    else
        return (meta.sum_index.range - length(meta.non_equal_indices)) * subTerm
    end
end
SymbolicUtils.substitute(Dsum::BasicSymbolic{IndexedAverageDoubleSum}, subs; fold = false) =
    SymbolicUtils.substitute(
        SymbolicUtils.metadata(Dsum)[IndexedAverageDoubleSum],
        subs;
        fold = fold,
    )
function SymbolicUtils.substitute(Dsum::IndexedAverageDoubleSum, subs; fold = false)
    inner = SymbolicUtils.substitute(Dsum.innerSum, subs; fold = fold)
    if SymbolicUtils._iszero(inner)
        return 0
    elseif inner isa BasicSymbolic{IndexedAverageSum}
        return IndexedAverageDoubleSum(inner, Dsum.sum_index, Dsum.non_equal_indices)
    else
        return (Dsum.sum_index.range - length(Dsum.non_equal_indices)) * inner
    end
end




#function to split sums into different clusters, keeping
#assume: all the indices that are not equal to the summation index are in the same Sum
#for anything else than the assumption, there needs an extra argument, to specify where or how the extra indices are handled
"""
    split_sums(term::SymbolicUtils.Symbolic,amount::Union{<:SymbolicUtils.Sym,<:Int64})
    split_sums(me::MeanfieldEquations,amount)

Function, that splits sums inside a given term. The sums are split into a number of equal sums, specified in the `amount` argument,
where in only one of the sums the dependencies for the indices (non equal indices) is considered.

# Arguments
*`me::MeanfieldEquations`: A [`MeanfieldEquations`](@ref) entity, which shall be evaluated, can also be any symbolic expression.
*`amount::Union{<:SymbolicUtils.Sym,<:Int64}`: A Number or Symbolic determining, in how many terms a sum is split

"""
function split_sums(
    term::SymbolicUtils.BasicSymbolic,
    ind::Index,
    amount::Union{<:SymbolicUtils.BasicSymbolic,<:Int64},
)
    if term isa Average
        return term
    end
    if SymbolicUtils.iscall(term)
        args = []
        op = SymbolicUtils.operation(term)
        for i = 1:length(arguments(term))
            push!(args, split_sums(arguments(term)[i], ind, amount))
        end
        if isempty(args)
            return 0
        elseif length(args) == 1
            return args[1]
        end
        return op(args...)
    end
    if term isa BasicSymbolic{IndexedAverageSum}
        meta = SymbolicUtils.metadata(term)[IndexedAverageSum]
        term_ = meta.term
        sumInd = meta.sum_index
        if isequal(ind, sumInd)
            ind2 =
                Index(sumInd.hilb, sumInd.name, (meta.sum_index.range/amount), sumInd.aon)
            term_ = change_index(term_, ind, ind2)
            if isempty(meta.non_equal_indices)
                return (amount)*IndexedAverageSum(term_, ind2, Index[])
            end
            extrasum = IndexedAverageSum(term_, ind2, meta.non_equal_indices)
            return extrasum + (amount-1)*IndexedAverageSum(term_, ind2, Index[])
        end
    end
    if term isa BasicSymbolic{IndexedAverageDoubleSum}
        Dsum = SymbolicUtils.metadata(term)[IndexedAverageDoubleSum]
        if isequal(ind, Dsum.sum_index)
            #create multiple doubleSums with the same innerSum
            ind2 = Index(
                Dsum.sum_index.hilb,
                Dsum.sum_index.name,
                (Dsum.sum_index.range/amount),
                Dsum.sum_index.aon,
            )
            if isempty(Dsum.non_equal_indices)
                return amount*IndexedAverageDoubleSum(Dsum.innerSum, ind2, Index[])
            end
            extraSum = IndexedAverageDoubleSum(Dsum.innerSum, ind2, Dsum.non_equal_indices)
            return extraSum + (amount-1)*Σ(Dsum.innerSum, ind2, Index[])
        elseif isequal(
            ind,
            SymbolicUtils.metadata(Dsum.innerSum)[IndexedAverageSum].sum_index,
        )
            #create doublesum of split innersums
            innerSums = split_sums(Dsum.innerSum, ind, amount)
            return IndexedAverageDoubleSum(
                innerSums,
                Dsum.sum_index,
                Dsum.non_equal_indices,
            )
        end
    end
    return term
end
split_sums(x::Symbolics.Equation, ind::Index, amount) =
    Symbolics.Equation(x.lhs, split_sums(x.rhs, ind, amount))
function split_sums(me::AbstractMeanfieldEquations, ind::Index, amount)
    newEqs = Vector{Union{Missing,Symbolics.Equation}}(missing, length(me.equations))
    for i = 1:length(me.equations)
        newEqs[i] = split_sums(me.equations[i], ind, amount)
    end
    newEqs = filter(x -> !isequal(x, missing), newEqs)
    vs = getLHS.(newEqs)
    varmap = make_varmap(vs, me.iv)
    return IndexedMeanfieldEquations(
        newEqs,
        me.operator_equations,
        vs,
        me.operators,
        me.hamiltonian,
        me.jumps,
        me.jumps_dagger,
        me.rates,
        me.iv,
        varmap,
        me.order,
    )
end
split_sums(x, ind, amount) = x

function scale(args...; kwargs...)
    return subst_reds_scale(scale_term(args...; kwargs...); kwargs...)
end
