#Base file for defining DoubleIndexedSums

"""

    IndexedDoubleSum <: QTerm

Defines a symbolic summation over another [`IndexedSingleSum`](@ref), using one [`Index`](@ref) entity. This corresponds to a double-summation over a multiplication of terms.

Fields:
======

* innerSum: A [`IndexedSingleSum`](@ref) entity.
* sumIndex: The index, for which the (outer) summation will go over.
* NEI: (optional) A vector of indices, for which the (outer) summation-index can not be equal with.

"""
struct IndexedDoubleSum <:QTerm
    innerSum::IndexedSingleSum
    sumIndex::Index
    NEI::Vector{Index}
    function IndexedDoubleSum(innerSum,sumIndex,NEI)
        if innerSum isa QAdd
            sums = []
            for arg in innerSum.arguments
                push!(sums, IndexedDoubleSum(arg,sumIndex,NEI))
            end
            isempty(sums) && return 0
            length(sums) == 1 && return sums[1]
            return +(sums...)
        end
        if innerSum isa IndexedSingleSum
            if innerSum.sumIndex == sumIndex
                error("summation index is the same as the index of the inner sum")
            else
                extraterm = 0
                NEI_ = copy(NEI)
                for index in getIndices(innerSum.term)
                    if sumIndex in innerSum.nonEqualIndices && isequal(index,innerSum.sumIndex)
                        (innerSum.sumIndex ∉ NEI_) && push!(NEI_,index)
                        continue
                    end
                    if index != sumIndex && index ∉ NEI && isequal(index.specHilb,sumIndex.specHilb)
                        extraterm = IndexedSingleSum(change_index(innerSum.term,sumIndex,index),innerSum.sumIndex,innerSum.nonEqualIndices)
                        push!(NEI_,index)
                    end 
                end
                if innerSum.term isa QMul
                    # put terms of the outer index in front
                    indicesToOrder = sort([innerSum.sumIndex,sumIndex],by=getIndName)
                    newargs = order_by_index(innerSum.term.args_nc,indicesToOrder)
                    qmul = 0
                    if length(newargs) == 1
                        qmul = *(innerSum.term.arg_c,newargs[1])
                    else
                        qmul = *(innerSum.term.arg_c,newargs...)
                    end
                    sort!(NEI_)
                    innerSum_ = IndexedSingleSum(qmul,innerSum.sumIndex,innerSum.nonEqualIndices)
                    if typeof(innerSum_) == IndexedSingleSum
                        if extraterm == 0
                            return new(innerSum_,sumIndex,NEI_)
                        end
                        return new(innerSum_,sumIndex,NEI_) + extraterm
                    else
                        return IndexedDoubleSum(innerSum_,sumIndex,NEI_)
                    end
                else
                    sort!(NEI)
                    return new(innerSum,sumIndex,NEI)
                end
            end
        else
            sort!(NEI)
            return IndexedSingleSum(innerSum,sumIndex,NEI)
        end
    end
end
IndexedDoubleSum(x,ind::Index) = IndexedDoubleSum(x,ind,Index[])

#In this constructor the NEI is considered so, that all indices given in ind are unequal to any of the NEI
function IndexedDoubleSum(term::QMul,ind::Vector{Index},NEI::Vector{Index})
    if length(ind) != 2
        error("Can only create Double-Sum with 2 indices!")
    end
    return IndexedDoubleSum(IndexedSingleSum(term,ind[1],NEI),ind[2],NEI)
end
function IndexedDoubleSum(term::QMul,outerInd::Index,innerInd::Index;non_equal::Bool=false)
    if non_equal
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
    if isequal(elem,0)
        return 0
    end
    return IndexedDoubleSum(elem*sum.innerSum,sum.sumIndex,sum.NEI)
end
function *(sum::IndexedDoubleSum,elem::SNuN)
    if isequal(elem,0)
        return 0
    end
    return IndexedDoubleSum(sum.innerSum*elem,sum.sumIndex,sum.NEI)
end
function *(qmul::QMul,sum::IndexedDoubleSum)
    newsum = sum
    for i = length(qmul.args_nc):-1:1 #elemt wise multiplication
        newsum = qmul.args_nc[i]*newsum
    end
    return qmul.arg_c*newsum
end
function *(sum::IndexedDoubleSum,qmul::QMul)
    newsum = sum
    for arg in qmul.args_nc #element wise multiplication
        newsum = newsum*arg
    end
    return newsum*qmul.arg_c
end
function *(elem,sum::IndexedDoubleSum)
    NEI = copy(sum.NEI)
    if typeof(elem) == IndexedOperator || typeof(elem) == IndexedVariable
        if elem.ind != sum.sumIndex && elem.ind ∉ NEI
            if ((sum.sumIndex.specHilb != sum.innerSum.sumIndex.specHilb) && isequal(elem.ind.specHilb,sum.sumIndex.specHilb))
                push!(NEI,elem.ind) 
                addterm = IndexedSingleSum(elem*change_index(sum.innerSum.term,sum.sumIndex,elem.ind),sum.innerSum.sumIndex,sum.innerSum.nonEqualIndices)
                return IndexedDoubleSum(elem*sum.innerSum,sum.sumIndex,NEI) + addterm
            end
        end
    end
    return IndexedDoubleSum(elem*sum.innerSum,sum.sumIndex,NEI)
end
function *(sum::IndexedDoubleSum,elem)
    NEI = copy(sum.NEI)
    if typeof(elem) == IndexedOperator || typeof(elem) == IndexedVariable
        if elem.ind != sum.sumIndex && elem.ind ∉ NEI
            if ((sum.sumIndex.specHilb != sum.innerSum.sumIndex.specHilb) && isequal(elem.ind.specHilb,sum.sumIndex.specHilb))
                push!(NEI,elem.ind)
                addterm = IndexedSingleSum(change_index(sum.innerSum.term,sum.sumIndex,elem.ind)*elem,sum.innerSum.sumIndex,sum.innerSum.nonEqualIndices)
                return IndexedDoubleSum(sum.innerSum*elem,sum.sumIndex,NEI) + addterm
            end
        end
    end
    return IndexedDoubleSum(sum.innerSum*elem,sum.sumIndex,NEI)
end
    
SymbolicUtils.istree(a::IndexedDoubleSum) = false
SymbolicUtils.arguments(a::IndexedDoubleSum) = SymbolicUtils.arguments(a.innerSum)
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
Base.isequal(a::IndexedDoubleSum,b::IndexedDoubleSum) = isequal(a.innerSum,b.innerSum) && isequal(a.sumIndex,b.sumIndex) && isequal(a.NEI,b.NEI)
_to_expression(x::IndexedDoubleSum) = :( IndexedDoubleSum($(_to_expression(x.innerSum)),$(x.sumIndex.name),$(x.sumIndex.rangeN),$(writeNEIs(x.NEI))))

function *(sum1::IndexedSingleSum,sum2::IndexedSingleSum; ind=nothing)
    if sum1.sumIndex != sum2.sumIndex
        term = sum1.term*sum2.term
        return IndexedDoubleSum(IndexedSingleSum(term,sum1.sumIndex,sum1.nonEqualIndices),sum2.sumIndex,sum2.nonEqualIndices)
    else
        if !(ind isa Index)
            error("Specification of an extra Index is needed!")
        end
        term2 = change_index(sum2.term,sum2.sumIndex,ind)
        return IndexedDoubleSum(IndexedSingleSum(sum1.term*term2,sum1.sumIndex,sum1.nonEqualIndices),ind,sum1.nonEqualIndices)

    end
end

