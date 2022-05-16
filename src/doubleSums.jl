#Base file for defining DoubleIndexedSums

include("indexing.jl")

struct IndexedDoubleSum <:QTerm
    innerSum::IndexedSingleSum
    sumIndex::Index
    NEI::Vector{Index}
    function IndexedDoubleSum(innerSum,sumIndex,NEI)
        if typeof(innerSum) == QAdd
            sums = []
            for arg in innerSum.arguments
                push!(sums, IndexedDoubleSum(arg,sumIndex,NEI))
            end
            return +(sums...)
        end
        if typeof(innerSum) <: SymbolicUtils.Add
            args = arguments(innerSum)
            sums = []
            for arg in args
                push!(sums, IndexedDoubleSum(arg,sumIndex,NEI))
            end
            return +(sums...)
        end
        if typeof(innerSum) == IndexedSingleSum
            if innerSum.sumIndex == sumIndex
                error("summation index is the same as the index of the inner sum")
            else
                extraterm = 0
                NEI_ = copy(NEI)
                for index in innerSum.nonEqualIndices
                    if index != sumIndex && index ∉ NEI 
                        extraterm = IndexedSingleSum(changeIndex(innerSum.term,sumIndex,index),innerSum.sumIndex,innerSum.nonEqualIndices)
                        push!(NEI_,index)
                    end 
                end
                indicesToOrder = sort([innerSum.sumIndex,sumIndex],by=getIndName)
                newargs = orderByIndex(innerSum.term.args_nc,indicesToOrder)
                qmul = merge_commutators(innerSum.term.arg_c,newargs)
                innerSum_ = IndexedSingleSum(qmul,innerSum.sumIndex,innerSum.nonEqualIndices)
                if typeof(innerSum_) == IndexedSingleSum
                    if extraterm == 0
                        return new(innerSum_,sumIndex,NEI_)
                    end
                    return new(innerSum_,sumIndex,NEI_) + extraterm
                else
                    return IndexedDoubleSum(innerSum_,sumIndex,NEI_)
                end
            end
        else
            return IndexedSingleSum(innerSum,sumIndex,NEI)
        end
    end
end
function IndexedDoubleSum(term::QAdd, sumIndex::Index,NEI::Vector{Index})
    sums = []
    for arg in term.arguments
        push!(sums,IndexedDoubleSum(arg,sumIndex,Index[]))
    end
    return +(sums...)
end

IndexedDoubleSum(term::QSym,ind::Index,NEI::Vector{Index}) = IndexedSingleSum(term,ind,NEI)
IndexedDoubleSum(term::SNuN,ind::Index,NEI::Vector{Index}) = IndexedSingleSum(term,ind,NEI)

#In this constructor the NEI is considered so, that all indices given in ind are unequal to any of the NEI
function IndexedDoubleSum(term::QMul,ind::Vector{Index},NEI::Vector{Index})
    if length(ind) != 2
        error("Can only create Double-Sum with 2 indices!")
    end
    return IndexedDoubleSum(IndexedSingleSum(term,ind[1],NEI),ind[2],NEI)
end
function IndexedDoubleSum(term::QMul,outerInd::Index,innerInd::Index,ne::Bool)
    if ne
        innerSum = IndexedSingleSum(term,innerInd,[outerInd])
        return IndexedDoubleSum(innerSum,outerInd,[])
    else
        innerSum = IndexedSingleSum(term,innerInd,[])
        return IndexedDoubleSum(innerSum,outerInd,[])
    end
end

hilbert(elem::IndexedDoubleSum) = hilbert(elem.sumIndex)
#adds
+(sum1::IndexedDoubleSum, sum2::IndexedDoubleSum) = (sum1.sumIndex == sum2.sumIndex) && checkInnerSums(sum1,sum2) ? 0 : QAdd([sum1,sum2])
+(elem::SNuN,sum::IndexedDoubleSum) = QAdd([elem,sum])
+(op::QNumber,sum::IndexedDoubleSum) = QAdd([op,sum])
function +(op::QAdd,sum::IndexedDoubleSum)
    args = copy(op.arguments)
    push!(args,sum)
    return QAdd(args)
end
#multiplications
function *(elem::SNuN, sum::IndexedDoubleSum)
    if elem == 0
        return 0
    end
    return IndexedDoubleSum(elem*sum.innerSum,sum.sumIndex,sum.NEI)
end
function *(sum::IndexedDoubleSum,elem::SNuN)
    if elem == 0
        return 0
    end
    return IndexedDoubleSum(sum.innerSum*elem,sum.sumIndex,sum.NEI)
end
function *(elem,sum::IndexedDoubleSum)
    NEI = copy(sum.NEI)
    if typeof(elem) == IndexedOperator || typeof(elem) == IndexedVariable
        if elem.ind != sum.sumIndex && elem.ind ∉ NEI
            push!(NEI,elem.ind)
        end
    end
    return IndexedDoubleSum(elem*sum.innerSum,sum.sumIndex,NEI)
end
*(elem::QNumber,sum::IndexedDoubleSum) = IndexedDoubleSum(elem*sum.innerSum,sum.sumIndex,sum.NEI)
*(sum::IndexedDoubleSum,elem::QNumber) = IndexedDoubleSum(sum.innerSum*elem,sum.sumIndex,sum.NEI)
    
SymbolicUtils.istree(a::IndexedDoubleSum) = false
checkInnerSums(sum1::IndexedDoubleSum, sum2::IndexedDoubleSum) = ((sum1.innerSum + sum2.innerSum) == 0)
reorder(dsum::IndexedDoubleSum,indexMapping::Vector{Tuple{Index,Index}}) = IndexedDoubleSum(reorder(dsum.innerSum,indexMapping),dsum.sumIndex,dsum.NEI)
#Base functions
function Base.show(io::IO,elem::IndexedDoubleSum)
    write(io,"Σ", "($(elem.sumIndex.name)=1:$(elem.sumIndex.rangeN))")
    if !(isempty(elem.NEI))
        write(io,"($(elem.sumIndex.name)≠")
        for i = 1:length(elem.NEI)
            write(io, "$(elem.NEI[i].name))")
        end
    end
    show(io,elem.innerSum)
end