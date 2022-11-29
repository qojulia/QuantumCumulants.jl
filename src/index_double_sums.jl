#Base file for defining DoubleIndexedSums

"""

    DoubleSum <: QTerm

Defines a symbolic summation over another [`SingleSum`](@ref), using one [`Index`](@ref) entity. This corresponds to a double-summation over a multiplication of terms.

Fields:
======

* innerSum: A [`SingleSum`](@ref) entity.
* sum_index: The index, for which the (outer) summation will go over.
* NEI: (optional) A vector of indices, for which the (outer) summation-index can not be equal with.

"""
struct DoubleSum{M} <:QTerm
    innerSum::SingleSum
    sum_index::Index
    NEI::Vector{Index}
    metadata::M
    function DoubleSum(innerSum::SingleSum,sum_index::Index,NEI::Vector,metadata)
        try
            return new{typeof(metadata)}(innerSum,sum_index,NEI,metadata)
        catch e
            println("Could not create DoubleSum with input: term= $(innerSum) ; sum_index=c$(sum_index) ; NEI= $(NEI) ; metadata= $(metadata)")
            rethrow(e)
        end

    end
end
function DoubleSum(innerSum::SingleSum,sum_index::Index,NEI;metadata=NO_METADATA)
        if innerSum.sum_index == sum_index
            error("summation index is the same as the index of the inner sum")
        else
            extraterm = 0
            NEI_ = copy(NEI)
            for index in get_indices(innerSum.term)
                if sum_index in innerSum.non_equal_indices && isequal(index,innerSum.sum_index)
                    (innerSum.sum_index ∉ NEI_) && push!(NEI_,index)
                    continue
                end
                if index != sum_index && index ∉ NEI && isequal(index.aon,sum_index.aon)
                    extraterm = SingleSum(change_index(innerSum.term,sum_index,index),innerSum.sum_index,innerSum.non_equal_indices)
                    push!(NEI_,index)
                end 
            end
            if innerSum.term isa QMul
                # put terms of the outer index in front
                indicesToOrder = sort([innerSum.sum_index,sum_index],by=getIndName)
                newargs = order_by_index(innerSum.term.args_nc,indicesToOrder)
                qmul = 0
                if length(newargs) == 1
                    qmul = *(innerSum.term.arg_c,newargs[1])
                else
                    qmul = *(innerSum.term.arg_c,newargs...)
                end
                sort!(NEI_)
                innerSum_ = SingleSum(qmul,innerSum.sum_index,innerSum.non_equal_indices)
                if innerSum_ isa SingleSum
                    if extraterm == 0
                        return DoubleSum(innerSum_,sum_index,NEI_,metadata)
                    end
                    return DoubleSum(innerSum_,sum_index,NEI_,metadata) + extraterm
                else
                    return DoubleSum(innerSum_,sum_index,NEI_;metadata=metadata)
                end
            else
                sort!(NEI)
                return DoubleSum(innerSum,sum_index,NEI,metadata)
            end
        end
end
function DoubleSum(innerSum::IndexedAdd,sum_index::Index,NEI;metadata=NO_METADATA)
    sums = [DoubleSum(arg,sum_index,NEI;metadata=metadata) for arg in arguments(innerSum)]
    isempty(sums) && return 0
    length(sums) == 1 && return sums[1]
    return +(sums...)
end
DoubleSum(x,ind::Index,NEI;metadata=NO_METADATA) = SingleSum(x,ind,NEI)
DoubleSum(x,ind::Index;metadata=NO_METADATA) = DoubleSum(x,ind,Index[])

#In this constructor the NEI is considered so, that all indices given in ind are unequal to any of the NEI
function DoubleSum(term::QMul,ind::Vector{Index},NEI::Vector{Index};metadata=NO_METADATA)
    if length(ind) != 2
        error("Can only create Double-Sum with 2 indices!")
    end
    return DoubleSum(SingleSum(term,ind[1],NEI),ind[2],NEI;metadata=metadata)
end
function DoubleSum(term::QMul,outerInd::Index,innerInd::Index;non_equal::Bool=false,metadata=NO_METADATA)
    if non_equal
        innerSum = SingleSum(term,innerInd,[outerInd])
        return DoubleSum(innerSum,outerInd,[];metadata=metadata)
    else
        innerSum = SingleSum(term,innerInd,[])
        return DoubleSum(innerSum,outerInd,[];metadata=metadata)
    end
end


hilbert(elem::DoubleSum) = hilbert(elem.sum_index)
#multiplications
*(elem::SNuN, sum::DoubleSum) = DoubleSum(elem*sum.innerSum,sum.sum_index,sum.NEI)
*(sum::DoubleSum,elem::SNuN) = DoubleSum(sum.innerSum*elem,sum.sum_index,sum.NEI)
*(sum::DoubleSum,qmul::QMul) = qmul.arg_c*(*(sum,qmul.args_nc...))
function *(qmul::QMul,sum::DoubleSum)

    sum_ = sum
    for i = length(qmul.args_nc):-1:1
        sum_ = qmul.args_nc[i] * sum_
    end
    return qmul.arg_c*sum_ 
end

function *(elem::IndexedObSym,sum::DoubleSum)

    NEI = copy(sum.NEI)
    if elem.ind != sum.sum_index && elem.ind ∉ NEI
        if ((sum.sum_index.aon != sum.innerSum.sum_index.aon) && isequal(elem.ind.aon,sum.sum_index.aon))
            push!(NEI,elem.ind) 
            addterm = SingleSum(elem*change_index(sum.innerSum.term,sum.sum_index,elem.ind),sum.innerSum.sum_index,sum.innerSum.non_equal_indices)
            return DoubleSum(elem*sum.innerSum,sum.sum_index,NEI) + addterm
        end
    end
    return DoubleSum(elem*sum.innerSum,sum.sum_index,NEI)
end
function *(sum::DoubleSum,elem::IndexedObSym)
    NEI = copy(sum.NEI)
    if elem.ind != sum.sum_index && elem.ind ∉ NEI
        if ((sum.sum_index.aon != sum.innerSum.sum_index.aon) && isequal(elem.ind.aon,sum.sum_index.aon))
            push!(NEI,elem.ind)
            addterm = SingleSum(change_index(sum.innerSum.term,sum.sum_index,elem.ind)*elem,sum.innerSum.sum_index,sum.innerSum.non_equal_indices)
            return DoubleSum(sum.innerSum*elem,sum.sum_index,NEI) + addterm
        end
    end
    return DoubleSum(sum.innerSum*elem,sum.sum_index,NEI)
end
*(sum::DoubleSum,x) = DoubleSum(sum.innerSum*x,sum.sum_index,sum.NEI)
*(x,sum::DoubleSum) = DoubleSum(x*sum.innerSum,sum.sum_index,sum.NEI) 


SymbolicUtils.istree(a::DoubleSum) = false
SymbolicUtils.arguments(a::DoubleSum) = SymbolicUtils.arguments(a.innerSum)
checkInnerSums(sum1::DoubleSum, sum2::DoubleSum) = ((sum1.innerSum + sum2.innerSum) == 0)
reorder(dsum::DoubleSum,indexMapping::Vector{Tuple{Index,Index}}) = DoubleSum(reorder(dsum.innerSum,indexMapping),dsum.sum_index,dsum.NEI)

#Base functions
function Base.show(io::IO,elem::DoubleSum)
    write(io,"Σ", "($(elem.sum_index.name)=1:$(elem.sum_index.range))")
    if !(isempty(elem.NEI))
        write(io,"($(elem.sum_index.name)≠")
        for i = 1:length(elem.NEI)
            write(io, "$(elem.NEI[i].name))")
        end
    end
    show(io,elem.innerSum)
end
Base.isequal(a::DoubleSum,b::DoubleSum) = isequal(a.innerSum,b.innerSum) && isequal(a.sum_index,b.sum_index) && isequal(a.NEI,b.NEI)
_to_expression(x::DoubleSum) = :( DoubleSum($(_to_expression(x.innerSum)),$(x.sum_index.name),$(x.sum_index.range),$(writeNEIs(x.NEI))))

function *(sum1::SingleSum,sum2::SingleSum; ind=nothing)
    if sum1.sum_index != sum2.sum_index
        term = sum1.term*sum2.term
        return DoubleSum(SingleSum(term,sum1.sum_index,sum1.non_equal_indices),sum2.sum_index,sum2.non_equal_indices)
    else
        if !(ind isa Index)
            error("Specification of an extra Index is needed!")
        end
        term2 = change_index(sum2.term,sum2.sum_index,ind)
        return DoubleSum(SingleSum(sum1.term*term2,sum1.sum_index,sum1.non_equal_indices),ind,sum1.non_equal_indices)

    end
end

