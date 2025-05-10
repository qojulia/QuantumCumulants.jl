function _new_operator(op::IndexedOperator, h, aon=nothing; kwargs...)
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
function _new_operator(nOp::NumberedOperator, h, aon=nothing; kwargs...)
    if !=(aon, nothing)
        return NumberedOperator(_new_operator(nOp.op, h, aon; kwargs...), nOp.numb)
    end
    return NumberedOperator(_new_operator(nOp.op, h; kwargs...), nOp.numb)
end
function _new_operator(sum::SingleSum, h, aon=nothing; kwargs...)
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
        _new_operator(sum.term, h, aon; kwargs...), newsum_index, newSumNonEquals
    )
end
function _new_operator(sum::DoubleSum, h, aon=nothing; kwargs...)
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
function _new_operator(sum::IndexedAverageSum, h, aon=nothing; kwargs...)
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
        _new_operator(sum.term, h, aon; kwargs...), newsum_index, newSumNonEquals
    )
end
function _new_operator(sym::BasicSymbolic{IndexedAverageSum}, h, aon=nothing; kwargs...)
    _new_operator(SymbolicUtils.metadata(sym)[IndexedAverageSum], h, aon; kwargs...)
end

function _new_operator(sum::IndexedAverageDoubleSum, h, aon=nothing; kwargs...)
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
function _new_operator(
    sym::BasicSymbolic{IndexedAverageDoubleSum}, h, aon=nothing; kwargs...
)
    _new_operator(SymbolicUtils.metadata(sym)[IndexedAvergeDoubleSum], h, aon; kwargs...)
end
function _new_operator(sym::BasicSymbolic{IndexedVariable}, h, aon=nothing; kwargs...)
    _new_operator(SymbolicUtils.metadata(sym)[IndexedVariable], h, aon; kwargs...)
end
function _new_operator(sym::BasicSymbolic{DoubleIndexedVariable}, h, aon=nothing; kwargs...)
    _new_operator(SymbolicUtils.metadata(sym)[DoubleIndexedVariable], h, aon; kwargs...)
end
function _new_operator(sym::BasicSymbolic{SpecialIndexedAverage}, h, aon=nothing; kwargs...)
    _new_operator(SymbolicUtils.metadata(sym)[SpecialIndexedAverage], h, aon; kwargs...)
end
function _new_operator(sym::SpecialIndexedAverage, h, aon=nothing; kwargs...)
    indexMap = [
        (_new_index(f, h, aon), _new_index(l, h, aon)) for (f, l) in sym.indexMapping
    ]
    return SpecialIndexedAverage(_new_operator(sym.term, h, aon; kwargs...), indexMap)
end
function _new_operator(sym::DoubleIndexedVariable, h, aon=nothing; kwargs...)
    if sym.ind1.hilb != h
        newInd1 = Index(h, sym.ind1.name, sym.ind1.range, sym.ind1.aon)
    end
    if sym.ind2.hilb != h
        newInd2 = Index(h, sym.ind2.name, sym.ind2.range, sym.ind2.aon)
    end
    return DoubleIndexedVariable(sym.name, newInd1, newInd2)
end
_new_operator(x::Vector, h, aon=nothing; kwargs...) = _new_operator.(x, h, aon; kwargs...)
function _new_operator(sym::IndexedVariable, h, aon=nothing; kwargs...)
    if sym.ind.hilb != h
        newInd = Index(h, sym.ind.name, sym.ind.range, sym.ind.aon)
    end
    return IndexedVariable(sym.name, newInd)
end

function _new_index(ind::Index, h, aon=nothing)
    if isnothing(aon)
        return Index(h, ind.name, ind.range, ind.aon)
    end
    return Index(h, ind.name, ind.range, aon)
end
function _new_indices(Inds::Vector, h)
    Inds_ = deepcopy(Inds)
    for i in 1:length(Inds)
        Inds_[i] = Index(h, Inds[i].name, Inds[i].range, Inds[i].aon)
    end
    return Inds_
end
