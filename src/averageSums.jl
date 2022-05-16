#Main file for manipulating indexed averages and sums over averages.
include("doubleSums.jl")

#some of these imports and usings can probably be removed (I just copied all of them)
import SymbolicUtils
import SymbolicUtils: substitute

import Symbolics
import TermInterface

import SciMLBase

import ModelingToolkit
const MTK = ModelingToolkit

using Combinatorics: partitions, combinations

using LaTeXStrings
using Latexify

const NO_METADATA = SymbolicUtils.NO_METADATA

struct AvgSum <: CNumber end

const SNuN = Union{<:SymbolicUtils.Symbolic{<:Number}, <:Number}
const Average_Sum = SymbolicUtils.Term{<:AvgSum}
const indornum = Union{<:Index,<:Int64}

const sum_average = begin # Symbolic function for averages
    T = SymbolicUtils.FnType{Tuple{CNumber}, AvgSum}
    SymbolicUtils.Sym{T}(:avg)
end

const symbolics_terms = Union{<:SymbolicUtils.Term,<:SymbolicUtils.Mul}

abstract type numberedVariable <: CNumber end

struct IndexedAverageSum <: CNumber
    term::symbolics_terms
    sumIndex::Index
    nonEqualIndices::Vector{indornum}
    function IndexedAverageSum(term,sumIndex,nonEqualIndices)
        if typeof(term) <: SymbolicUtils.Add
            newterm = 0.0
            for elem in arguments(term)
                newterm = newterm + IndexedAverageSum(elem,sumIndex,nonEqualIndices)
            end
            return newterm
        end
        if !(typeof(term) <: symbolics_terms)
            return (sumIndex.rangeN - length(nonEqualIndices)) * term
        end
        prefact = 1.0    #assume term is of type SymbolicUtils.Mul
        metadata = new(term,sumIndex,nonEqualIndices)
        neis_sym = ""
        if !(isempty(nonEqualIndices))
            neis_sym = string("(",neis_sym)
            neis_sym = string(neis_sym, "$(sumIndex.name)≠")
            neis_sym = string(neis_sym, writeNEIs(nonEqualIndices))
            neis_sym = string(neis_sym,")")
        end
        if typeof(arguments(term)[1]) <: Number # put numbers outside of sum (for easier evaluation)
            prefact = arguments(term)[1]
            deleteat!(arguments(term),1)
        end
        return prefact*SymbolicUtils.Sym{Parameter, IndexedAverageSum}(Symbol("∑($(sumIndex.name)=1:$(sumIndex.rangeN))$(neis_sym)$(term)"), metadata) #Symbol("∑($(sumIndex.name)=1:$(sumIndex.rangeN))$(neis_sym)$(term)")
    end
end

struct IndexedAverageDoubleSum <: CNumber
    innerSum::Sym{Parameter, IndexedAverageSum}
    sumIndex::Index
    nonEqualIndices::Vector{indornum}
    function IndexedAverageDoubleSum(term,sumIndex,nonEqualIndices)
        if typeof(term) <: SymbolicUtils.Add
            newterm = 0.0
            for elem in arguments(term)
                newterm = newterm + IndexedAverageDoubleSum(elem,sumIndex,nonEqualIndices)
            end
            return newterm
        end
        if typeof(term) == Sym{Parameter,IndexedAverageSum}
            metadata = new(term,sumIndex,nonEqualIndices)
            neis_sym = ""
            if !(isempty(nonEqualIndices))
                neis_sym = string("(",neis_sym)
                neis_sym = string(neis_sym, "$(sumIndex.name)≠")
                neis_sym = string(neis_sym, writeNEIs(nonEqualIndices))
                neis_sym = string(neis_sym,")")
            end
            return SymbolicUtils.Sym{Parameter, IndexedAverageDoubleSum}(Symbol("∑($(sumIndex.name):=1:$(sumIndex.rangeN))$(neis_sym)$(String(term.name))"), metadata)
        end
        if typeof(term) <: SymbolicUtils.Mul
            args = arguments(term)
            param = 1.0
            if typeof(args[1]) <: Number #put numbers out in front
                param = args[1]
                deleteat!(args,1)
            end
            if length(args) == 1 && typeof(args[1]) == Sym{Parameter, IndexedAverageSum}
                return param*IndexedAverageDoubleSum(args[1],sumIndex,nonEqualIndices)
            else
                for arg in args
                    if typeof(arg) == Sym{Parameter, IndexedAverageSum}
                        error("cannot convert multiplication of Sums into double sums")
                        return 0
                    end
                end
            end
        end
        return IndexedAverageSum(term,sumIndex,nonEqualIndices)
    end
end

#For representing in average terms
struct NumberedOperator <:QNumber
    op::Transition
    numb::Int64
    function NumberedOperator(op,numb)
        if numb <= 0
            error("can not create numbered-operator with negative or 0 number")
            return 0
        end
        if typeof(op) <: SNuN
            return op
        end
        return new(op,numb)
    end
end
#Variables
struct SingleNumberedVariable <: numberedVariable
    name::Symbol
    numb::Int64
    function SingleNumberedVariable(name,numb)
        metadata=source_metadata(:Parameter, name)
        s = SymbolicUtils.Sym{Parameter, typeof(metadata)}(Symbol("$(name)$(numb)"), metadata)
        return SymbolicUtils.setmetadata(s, MTK.MTKParameterCtx, true)
    end
end
struct DoubleNumberedVariable <: numberedVariable
    name::Symbol
    numb1::indornum
    numb2::indornum
    function DoubleNumberedVariable(name,numb1,numb2)
        if typeof(numb1) == typeof(numb2) && typeof(numb1) == Int64
            metadata = source_metadata(:Parameter, name)
            s = SymbolicUtils.Sym{Parameter, typeof(metadata)}(Symbol("$(name)$(numb1)$(numb2)"), metadata)
            return SymbolicUtils.setmetadata(s, MTK.MTKParameterCtx, true)
        else
            return SymbolicUtils.Sym{Parameter, numberedVariable}(Symbol("$(name)$(numb1)$(numb2)"), new(name,numb1,numb2))
        end
    end
end
#Special averages
struct NumberedIndexedAverage <: CNumber
    term::symbolics_terms
    indexMapping::Vector{Tuple{Index,Int64}} #not to be confused with later definition of index mapping here if l => 2 means l ≠ 2
    function NumberedIndexedAverage(term,indexMapping)
        if length(indexMapping) == 0
            return term
        elseif (typeof(term) <: QMul && term.arg_c === 0) || SymbolicUtils._iszero(term)
            return 0
        else
            Neqs = writeNeqs(indexMapping)
            return SymbolicUtils.Sym{Parameter, NumberedIndexedAverage}(Symbol("$(Neqs)$(term)"), new(term,indexMapping))
        end
    end
end
struct SpecialIndexedAverage <: CNumber #An average-Term with special condition, for example l ≠ k; needed for correct calculus of Double indexed Sums
    term::symbolics_terms
    indexMapping::Vector{Tuple{indornum,indornum}}
    function SpecialIndexedAverage(term,indexMapping)
        if isempty(indexMapping)
            return term
        end
        if typeof(term) <: SymbolicUtils.Add
            adds = []
            for arg in arguments(term)
                push!(adds,SpecialIndexedAverage(arg,indexMapping))
            end
            return +(adds...)
        end
        metadata = new(term,indexMapping)
        neis = writeIndexNEIs(indexMapping)
        return SymbolicUtils.Sym{Parameter, SpecialIndexedAverage}(Symbol("$(neis)$(term)"), metadata)
    end
end
average(indOp::IndexedOperator) = _average(indOp)
average(x::SpecialIndexedTerm) = SpecialIndexedAverage(average(x.term),x.indexMapping)

function average(indSum::IndexedSingleSum; kwargs...)
    return IndexedAverageSum(average(indSum.term),indSum.sumIndex,indSum.nonEqualIndices)
end
function undo_average(a::IndexedAverageSum)
    return IndexedSingleSum(undo_average(a.term),a.sumIndex,a.nonEqualIndices)
end
function undo_average(a::Sym{Parameter,IndexedAverageSum})
    return undo_average(a.metadata)
end
function average(indDSum::IndexedDoubleSum)
    return IndexedAverageDoubleSum(average(indDSum.innerSum),indDSum.sumIndex,indDSum.NEI)
end
function undo_average(a::Sym{Parameter,IndexedAverageDoubleSum})
    return undo_average(a.metadata)
end
function undo_average(a::IndexedAverageDoubleSum)
    return IndexedDoubleSum(undo_average(a.innerSum),a.sumIndex,a.nonEqualIndices)
end

#define calculus for numbered operators -> break it down into QNuber multiplication
*(numOp::NumberedOperator, qmul::QMul) = QMul(qmul.arg_c,vcat(numOp,qmul.args_nc))
*(qmul::QMul, numOp::NumberedOperator) = QMul(qmul.arg_c,vcat(qmul.args_nc,numOp))
*(numOp1::NumberedOperator,numOp2::NumberedOperator) = numOp1.numb == numOp2.numb ? NumberedOperator(numOp1.op*numOp2.op,numOp1.numb) : QMul(1,[numOp1,numOp2])
*(elem::SNuN, numOp::NumberedOperator) = QMul(elem,[numOp])
*(numOp::NumberedOperator,elem::SNuN) = QMul(elem,[numOp])
*(a::Create,b::NumberedOperator) = merge_commutators(1,[a,b])
*(b::NumberedOperator,a::Create) = merge_commutators(1,[b,a])
*(a::Destroy,b::NumberedOperator) = merge_commutators(1,[a,b])
*(b::NumberedOperator,a::Destroy) = merge_commutators(1,[b,a])
*(a::IndexedOperator,b::NumberedOperator) = merge_commutators(1,[a,b])
*(b::NumberedOperator,a::IndexedOperator) = merge_commutators(1,[b,a])

get_order(a::Sym{Parameter,IndexedAverageSum}) = get_order(a.metadata.term)
cumulant_expansion(a::IndexedAverageSum,order::Int) = IndexedAverageSum(simplifyMultiplifcation(cumulant_expansion(a.term,order)),a.sumIndex,a.nonEqualIndices) #not used (?)
SymbolicUtils.istree(a::IndexedAverageSum) = false
SymbolicUtils.promote_symtype(::typeof(sum_average), ::Type{IndexedAverageSum}) = AvgSym
SymbolicUtils.promote_symtype(::typeof(+), ::Type{IndexedAverageSum}) = AvgSym
SymbolicUtils._iszero(a::IndexedAverageSum) = SymbolicUtils._iszero(a.term)
SymbolicUtils._isone(a::IndexedAverageSum) = SymbolicUtils._isone(a.term)
IndexedAverageSum(x::Number) = x
SymbolicUtils.istree(a::SymbolicUtils.Sym{Parameter,SpecialIndexedAverage}) = false
SymbolicUtils.istree(a::IndexedAverageDoubleSum) = false
SymbolicUtils.istree(a::SymbolicUtils.Sym{Parameter,IndexedAverageDoubleSum}) = false
get_order(a::Sym{Parameter,IndexedAverageDoubleSum}) = get_order(a.metadata.innerSum)
get_order(a::SymbolicUtils.Sym{Parameter,SpecialIndexedAverage}) = get_order(a.metadata.term)

average(x::NumberedOperator) = _average(x)
hilbert(x::NumberedOperator) = hilbert(x.op)
Base.adjoint(x::NumberedOperator) = NumberedOperator(Base.adjoint(x.op),x.numb)
has_cluster(x::NumberedOperator) = has_cluster(x.op)
acts_on(x::NumberedOperator) = acts_on(x.op) + x.numb
get_order(x::NumberedOperator) = get_order(x.op)

#Functions for easier symbol creation in Constructor
function writeNEIs(neis::Vector{indornum})
    syms = ""
    for i = 1:length(neis)
        syms = typeof(neis[i]) == Index ? join([syms,neis[i].name]) : join([syms,neis[i]])
        if i != length(neis)
            syms = join([syms,","])
        end
    end
    return syms
end
function writeNEIs(neis::Vector{Index})
    syms = ""
    for i = 1:length(neis)
        syms = join([syms,neis[i].name])
        if i != length(neis)
            syms = join([syms,","])
        end
    end
    return syms
end
function writeIndexNEIs(neis::Vector{Tuple{indornum,indornum}})
    syms = ""
    syms = join([syms,"("])
    for i = 1:length(neis)
        if typeof(first(neis[i])) == Index
            syms = join([syms,first(neis[i]).name])
        else
            syms = join([syms,first(neis[i])])
        end
        syms = join([syms,"≠"])
        if typeof(last(neis[i])) == Index
            syms = join([syms,last(neis[i]).name])
        else
            syms = join([syms,last(neis[i])])
        end
        if i != length(neis)
            syms = join([syms,","])
        end
    end
    syms = join([syms,")"])
    return syms
end
writeIndexNEIs(neis::Vector{Tuple{Index,Index}}) = writeIndexNEIs(convert(Vector{Tuple{indornum,indornum}},neis))
function writeNeqs(vec::Vector{Tuple{Index,Int64}})
    syms = ""
    for i = 1:length(vec)
        syms = join([syms, "("])
        syms = join([syms,first(vec[i]).name])
        syms = join([syms,"≠"])
        syms = join([syms,last(vec[i])])
        syms = join([syms,")"])
    end
    return syms
end

#Base functions
function Base.hash(a::IndexedAverageSum, h::UInt)
    return hash(IndexedAverageSum, hash(a.term, hash(a.sumIndex, hash(a.nonEqualIndices,h))))
end 
Base.isless(a::IndexedAverageSum,b::IndexedAverageSum) = a.sumIndex < b.sumIndex
function Base.isequal(a::IndexedAverageSum, b::IndexedAverageSum)
    isequal(a.sumIndex,b.sumIndex) || return false
    isequal(a.term, b.term) || return false
    isequal(a.nonEqualIndices,b.nonEqualIndices) || return false
    return true
end
Base.isequal(a::IndexedAverageSum,b) = false
Base.isequal(::IndexedAverageSum, ::SymbolicUtils.Symbolic) = false
function Base.isequal(nVal1::Sym{Parameter,numberedVariable},nVal2::Sym{Parameter,numberedVariable}) 
    if typeof(nVal1) == typeof(nVal2) && typeof(nVal1) == SingleNumberedVariable
        return (nVal1.name == nVal2.name) && (nVal1.numb == nVal2.numb)
    elseif typeof(nVal1) == typeof(nVal2) && typeof(nVal1) == DoubleNumberedVariable
        return (nVal1.name == nVal2.name) && (nVal1.numb1 == nVal2.numb1) && (nVal1.numb2 == nVal2.numb2)
    end
    return false
end
Base.:(==)(nVal1::Sym{Parameter,numberedVariable},nVal2::Sym{Parameter,numberedVariable}) = (nVal1.name == nVal2.name) && (nVal1.numb == nVal2.numb)
function cumulant_expansion(x::SymbolicUtils.Sym{Parameter,IndexedAverageSum},order::Integer;simplify=true,kwargs...)
    sum = x.metadata
    return IndexedAverageSum(simplifyMultiplication(cumulant_expansion(sum.term,order;simplify,kwargs...)),sum.sumIndex,sum.nonEqualIndices)
end
function cumulant_expansion(x::SymbolicUtils.Sym{Parameter,IndexedAverageDoubleSum},order::Integer;simplify=true,kwargs...)
    inner = cumulant_expansion(x.metadata.innerSum,order;simplify,kwargs...)
    return IndexedAverageDoubleSum(inner,x.metadata.sumIndex,x.metadata.nonEqualIndices)
end
cumulant_expansion(a::SymbolicUtils.Sym{Parameter,SpecialIndexedAverage},order::Int;simplify=true,kwargs...) = SpecialIndexedAverage(cumulant_expansion(a.metadata.term,order;simplify=true,kwargs...),a.metadata.indexMapping)
SymbolicUtils.arguments(op::Sym{Parameter,IndexedAverageSum}) = arguments(op.metadata)
SymbolicUtils.arguments(op::IndexedAverageSum) = op.term
SymbolicUtils.arguments(op::Sym{Parameter, IndexedAverageDoubleSum}) = arguments(op.metadata)
SymbolicUtils.arguments(op::IndexedAverageDoubleSum) = op.innerSum

#evaluate functions are not needed anymore, are now included into the insert functions, still kept here for references
function evaluateTerm(indDSum::SymbolicUtils.Sym{Parameter, IndexedAverageDoubleSum}, indMap::Dict{Index,Int64})
    innerSum = indDSum.metadata.innerSum
    return IndexedAverageSum(evaluateTerm(innerSum,indMap))
end
function evaluateTerm(term::SymbolicUtils.Term{AvgSym, Nothing}, indMap::Dict{Index,Int64})
    args = arguments(term)
    newargs = []
    args_ = 0
    if typeof(args[1]) <: QMul #Qmul
        for arg in args[1].args_nc
            if typeof(arg) == IndexedOperator && arg.ind ∈ keys(indMap) #index of operator in the indexmapping
                push!(newargs, NumberedOperator(arg.op,indMap[arg.ind]))
            else    #index of operator not in the indexmapping, or has no index
                push!(newargs,arg) 
            end
        end
    else #single op
        if typeof(args[1]) == IndexedOperator && args[1].ind ∈ keys(indMap)
            return average(NumberedOperator(args[1].op,indMap[args[1].ind]))
        else
            return average(args[1])
        end
    end
    return average(*(newargs...))
end
#Evaluate SingleSum
function evaluateTerm(indSum::SymbolicUtils.Sym{Parameter, IndexedAverageSum}, indMap::Dict{Index,Int64})
    adds = []
    term = indSum.metadata.term
    args = arguments(term)
    sumInd = indSum.metadata.sumIndex
    NEIs = indSum.metadata.nonEqualIndices
    NEIvals = []
    for ind in NEIs #ordnung anzahl indices in lhs
        if ind in keys(indMap)
            push!(NEIvals,indMap[ind])
        end
    end
    for i = 1:sumInd.rangeN #ordnung N <- !!
        if i in NEIvals
            continue
        end 
        NEValMapping = Tuple{Index,Int64}[]
        args_after = []
        for arg in args
            meta = 0
            if typeof(arg) == SymbolicUtils.Sym{Parameter, IndexedVariable} # symbols like gₖ 
                meta = arg.metadata
                if meta.ind == sumInd # k = j
                    push!(args_after,SingleNumberedVariable(meta.name,i))
                elseif meta.ind ∈ keys(indMap) # k => Number
                    push!(args_after,SingleNumberedVariable(meta.name,indMap[meta.ind]))
                else #neither
                    push!(args_after, arg)
                end
            elseif typeof(arg) == SymbolicUtils.Term{AvgSym, Nothing} #At this stage the terms making it into the if clausle are average symbols like ⟨a∗σ21j⟩
                indices_  = getIndices(arg) #all indices that are in the term
                for ind in indices_ #all indices that are not in the indMap and are not the summation index
                    if ind == sumInd || ind ∈ keys(indMap)
                        continue
                    else
                        push!(NEValMapping,(ind,i)) #save those indices for later (create terms like (k≠1)<σₖ*σ₁>)
                    end
                end
                newargs = []
                args_ = 0
                avrg_args = arguments(arg)[1]
                if typeof(avrg_args) <: QMul
                    args_ = avrg_args.args_nc
                else
                    args_ = [avrg_args]
                end
                for avrg_arg in args_ #loop over arguments inside the average
                    if typeof(avrg_arg) == IndexedOperator && avrg_arg.ind == sumInd #this is the interesting case, index of operator the same as the summation index
                        push!(newargs, NumberedOperator(avrg_arg.op,i))
                    elseif typeof(avrg_arg) == IndexedOperator && avrg_arg.ind ∈ keys(indMap) #index of operator in the indexmapping
                        push!(newargs, NumberedOperator(avrg_arg.op,indMap[avrg_arg.ind]))
                    else    #index of operator not in the indexmapping, or has no index
                        push!(newargs,avrg_arg) 
                    end
                end
                if length(newargs) == 1
                    push!(args_after, average(newargs[1]))
                else
                    push!(args_after, average(*(newargs...))) #stick everything back together
                end
            else
                push!(args_after, arg)
            end
        end
         push!(adds, NumberedIndexedAverage(*(args_after...),NEValMapping))#multiply everything, that was originaly in the sum
        
    end
    if isempty(adds)
        return 0
    end
    return +(adds...) #add everything back up
end
function evaluateTerm(term::SymbolicUtils.Sym{Parameter, IndexedVariable}, indMap::Dict{Index,Int64})
    indVar = term.metadata
    if indVar.ind ∈ keys(indMap)
        return SingleNumberedVariable(indVar.name, indMap[indVar.ind])
    end
    return term
end
function evaluateTerm(term::SymbolicUtils.Mul, indMap::Dict{Index,Int64})
    args = arguments(term)
    newterms = []
    for arg in args
        push!(newterms, evaluateTerm(arg, indMap))
    end
    return *(newterms...)
end
function evaluateTerm(term::SymbolicUtils.Pow, indMap::Dict{Index,Int64})
    avg = arguments(term)[1]
    pow = arguments(term)[2]
    return evaluateTerm(avg,indMap)^pow
end
function evaluateTerm(term::Vector,indMap::Dict{Index,Int64})
    result = []
    for op in term
        newops = []
        if op.ind ∉ keys(indMap)
            for i = 1:op.ind.rangeN
                push!(newops, NumberedOperator(op.op,i))
            end
        else
            push!(newops,NumberedOperator(op.op,indMap[op.ind]))
        end
        vcat(result,newops)
    end
    return result
end
function evaluateTerm(x,indMap)
    try
        args = arguments(x)
        f = operation(x)
        newterms = []
        for arg in args
            push!(newterms, evaluateTerm(arg, indMap))
        end
        return f(newterms...)
    catch e
        return x
    end
end
evaluateTerm(term::SNuN, indMap::Dict{Index,Int64}) = term
evaluateTerm(x) = evaluateTerm(x, Dict{Index,Int64}())
#function that evaluates several equations given by one indexed equation, using index into number substituation
#this function requires, that the range of an index is a number, not a symbol
function evaluateEquation(eq::Equation)
    lhs = eq.lhs
    rhs = eq.rhs
    rhs_ = []
    lhs_ = []
    if containsIndexedOps(lhs)
        ind_ = getIndices(lhs)
        maxRange = 1
        for ind in ind_
            maxRange = maxRange * ind.rangeN
        end
        for i = 1:maxRange
            tuples = []
            factor = 1
            for k = 1:length(ind_)
                if k != 1
                    factor = factor * ind_[k-1].rangeN
                end
                push!(tuples, (ind_[k] => (((i ÷ factor) % ind_[k].rangeN)+1) ))    #(((i ÷ factor) % ind_[k])+1) corresponds in this case for the assotiated index the corresponding index in the initial index-list should have for the i-th equation
            end
            Mapping = Dict(tuples) #create correstponding index-mapping

            if length(unique(values(Mapping))) != length(values(Mapping)) #true if values(Mapping) contains duplicated values
                continue
            end

            push!(lhs_, orderTermsByNumber(evaluateTerm(lhs, Mapping)))
            push!(rhs_, orderTermsByNumber(evaluateTerm(rhs, Mapping)))
        end
    else
        rhs_ = [orderTermsByNumber(evaluateTerm(rhs,Dict{Index,Int64}()))]
        lhs_ = [orderTermsByNumber(evaluateTerm(lhs,Dict{Index,Int64}()))]
    end
    eqs = [Symbolics.Equation(l,r) for (l,r)=zip(lhs_,rhs_)]
    return eqs
end
#function for evaluating indexed-meanfield equations
function evaluateMeanfieldEquations(me::MeanfieldEquations)
    newEqs = []
    newOpEqs = []
    eqs = copy(me.equations)
    for eq in eqs
        evals = evaluateEquation(eq) #returns vector of equations -> need to iterate again
        for eq_ in evals
            if (eq_.lhs ∉ getLHS.(newEqs)) && (_conj(eq_.lhs) ∉ getLHS.(newEqs)) && (orderTermsByNumber(_conj(eq_.lhs)) ∉ orderTermsByNumber.(getLHS.(newEqs))) && (orderTermsByNumber(eq.lhs) ∉ orderTermsByNumber.(getLHS.(newEqs)))
                push!(newEqs, eq_)
            end
        end
    end
    vs = []
    for eq in newEqs
        push!(vs,eq.lhs)
    end
    varmap = make_varmap(vs, me.iv)
    newME = MeanfieldEquations(newEqs,newOpEqs,vs,me.operators,me.hamiltonian,me.jumps,me.jumps_dagger,me.rates,me.iv,varmap,me.order)
    return newME
end

#this is the new method, insert values directly into the average before calculating anything, simplifies evaluation afterwards extremely
#function for inserting index, k -> 1,2,...,N
function insertIndex(sum::SymbolicUtils.Sym{Parameter,IndexedAverageSum}, ind::Index, value::Int64)
    if ind == sum.metadata.sumIndex
        error("cannot exchange summation index with number!")
    end
    if ind in sum.metadata.nonEqualIndices
        newNEI = filter(x-> x != ind,sum.metadata.nonEqualIndices)
        push!(newNEI,value)
        return IndexedAverageSum(insertIndex(sum.metadata.term,ind,value),sum.metadata.sumIndex,newNEI)
    else
        return IndexedAverageSum(insertIndex(sum.metadata.term,ind,value),sum.metadata.sumIndex,sum.metadata.nonEqualIndices)
    end
end
function insertIndex(sum::SymbolicUtils.Sym{Parameter,IndexedAverageDoubleSum}, ind::Index,value::Int64)
    inner = insertIndex(sum.metadata.innerSum,ind,value)
    return IndexedAverageDoubleSum(inner,sum.metadata.sumIndex,sum.metadata.nonEqualIndices)
end
function insertIndex(term::SymbolicUtils.Mul, ind::Index, value::Int64)
    args = []
    for arg in arguments(term)
        push!(args,insertIndex(arg,ind,value))
    end
    return *(args...)
end
function insertIndex(term::SymbolicUtils.Add,ind::Index,value::Int64)
    args = []
    for arg in arguments(term)
        push!(args,insertIndex(arg,ind,value))
    end
    return +(args...)
end
function insertIndex(term::SymbolicUtils.Pow,ind::Index,value::Int64)
    return insertIndex(arguments(term)[1],ind,value)^(arguments(term)[2])
end
function insertIndex(term::SymbolicUtils.Term{AvgSym,Nothing},ind::Index,value::Int64)
    newargs = []
    if typeof(arguments(term)[1]) <: QMul
        for arg in arguments(term)[1].args_nc
            push!(newargs,insertIndex(arg,ind,value))
        end
        return average(QMul(1,newargs))
    else
        return average(insertIndex(arguments(term)[1],ind,value))
    end
end
function insertIndex(term_::Sym{Parameter, DoubleIndexedVariable},ind::Index,value::Int64)
    term = term_.metadata
    if term.ind1 == ind && term.ind2 == ind
        return DoubleNumberedVariable(term.name,value,value)
    elseif term.ind1 == ind
        return DoubleNumberedVariable(term.name,value,term.ind2)
    elseif term.ind2 == ind
        return DoubleNumberedVariable(term.name,term.ind1,value)
    end
    return term_
end
function insertIndex(term::SymbolicUtils.Sym{Parameter,numberedVariable},ind::Index,value::Int64)
    if typeof(term.metadata) == SingleNumberedVariable
        return term
    end
    data = term.metadata
    if typeof(data.numb1) == Index && data.numb1 == ind
        return DoubleNumberedVariable(data.name,value,data.numb2)
    elseif typeof(data.numb2) == Index && data.numb2 == ind
        return DoubleNumberedVariable(data.name,data.numb1,value)
    end
    return term
end
function insertIndex(term::SymbolicUtils.Sym{Parameter,SpecialIndexedAverage},ind::Index,value::Int64)
    meta = term.metadata
    newterm = insertIndex(meta.term,ind,value)
    newMapping = Tuple{indornum,indornum}[]
    for tuple in meta.indexMapping
        if first(tuple) == ind
            if last(tuple) == value
                return 0
            end
            push!(newMapping,(value,last(tuple)))
        elseif last(tuple) == ind
            if first(tuple) == value
                return 0
            end
            push!(newMapping,(first(tuple),value))
        else
            push!(newMapping,tuple)
        end
    end
    filter!(x -> !(typeof(first(x)) == Int64 && typeof(last(x)) == Int64),newMapping)
    return SpecialIndexedAverage(newterm,newMapping)
end
insertIndex(eq::Symbolics.Equation,ind::Index,value::Int64) = Symbolics.Equation(insertIndex(eq.lhs,ind,value),insertIndex(eq.rhs,ind,value))
insertIndex(term::IndexedOperator,ind::Index,value::Int64) = term.ind == ind ? NumberedOperator(term.op,value) : term
insertIndex(term::SymbolicUtils.Sym{Parameter,IndexedVariable},ind::Index,value::Int64) = term.metadata.ind == ind ? SingleNumberedVariable(term.metadata.name,value) : term
insertIndex(x,ind::Index,value::Int64) = x
function insertIndices(eq::Symbolics.Equation,map::Dict{Index,Int64};mapping::Dict{Symbol,Int64}=Dict{Symbol,Int64}())
    eq_ = eq
    while !isempty(map)
        pair = first(map)
        eq_ = insertIndex(eq_,first(pair),last(pair))
        delete!(map,first(pair))
    end
    return evalEq(eq_;mapping) #return finished equation
end
function evalEquation(eq::Symbolics.Equation,arr,indices;mapping::Dict{Symbol,Int64}=Dict{Symbol,Int64}())
    if !(isempty(indices))
        eqs = Vector{Any}(nothing, length(arr))
        #Threads.@threads 
        for i=1:length(arr)
            dict = Dict(indices .=> arr[i])
            eq_ = orderTermsByNumber(insertIndices(eq,dict;mapping))
            eqs[i] = eq_
        end
        return filter(x -> x != nothing,eqs)
    else
        return [evalEq(eq;mapping)]
    end
end
function evalME(me::MeanfieldEquations;mapping::Dict{Symbol,Int64}=Dict{Symbol,Int64}())#this is still pretty slow
    indices = nothing
    for eq in me.equations
        if containsIndexedOps(eq.lhs) && length(getIndices(eq.lhs)) == me.order
            indices = getIndices(eq.lhs)
            break
        end
    end
    #the maximum nummber of equations should be something like: order*(numberOfAtoms)^(order)
    #in case for 2nd order and 30 atoms there were 1663 equations
    #for the allocation it is assumed, that the first index given has the highest range
    range = 0
    if indices[1].rangeN in keys(mapping)
        range = mapping[indices[1].rangeN]
    else
        range = indices[1].rangeN
    end
    newEqs = Vector{Union{Missing,Symbolics.Equation}}(missing,length(indices)*(range+5)^(length(indices))) #preallocation for newEqs
    ranges = []
    arrays = []
    for ind in indices
        if ind.rangeN in keys(mapping)
            push!(ranges,1:mapping[ind.rangeN])
        else
            push!(ranges,1:ind.rangeN)
        end
        push!(arrays,unique(sort.(collect.(filter(x -> length(x) == length(unique(x)),collect(Iterators.product(ranges...)))))))
    end
    # Threads.@threads
    for i=1:length(me.equations)
        ord = length(getIndices(me.equations[i].lhs))
        if ord == 0
            evals = evalEquation(me.equations[i],[],[];mapping)
        else
            evals = evalEquation(me.equations[i],arrays[ord],indices[1:ord];mapping)
        end
        for eq_ in evals #might be able to reduce this loop into one of the other loops
            if (eq_.lhs ∉ getLHS.(newEqs)) && (_conj(eq_.lhs) ∉ getLHS.(newEqs))
                push!(newEqs, eq_)
            end
        end
        #newEqs = vcat(newEqs,evals)
        #vs = vcat(vs, getLHS.(evals))
    end
    newEqs = filter(x -> !isequal(x,missing), newEqs)
    vs = getLHS.(newEqs)
    varmap = make_varmap(vs, me.iv)
    return MeanfieldEquations(newEqs,me.operator_equations,vs,me.operators,me.hamiltonian,me.jumps,me.jumps_dagger,me.rates,me.iv,varmap,me.order)
end
#TODO change the adds to a sum of adds -> instead of +(adds...) do something like sum(term_i for i = 1:N)
function evalTerm(sum::SymbolicUtils.Sym{Parameter,IndexedAverageSum};mapping::Dict{Symbol,Int64}=Dict{Symbol,Int64}())
    adds = []
    rangeEval = 0
    if sum.metadata.sumIndex.rangeN in keys(mapping)
        rangeEval = mapping[sum.metadata.sumIndex.rangeN]
    else
        rangeEval = sum.metadata.sumIndex.rangeN
    end
    for i = 1:rangeEval
        if i in sum.metadata.nonEqualIndices
            continue
        end
        push!(adds,orderTermsByNumber(insertIndex(sum.metadata.term,sum.metadata.sumIndex,i)))
    end
    if isempty(adds)
        return 0
    elseif length(adds) == 1
        return adds[1]
    end
    return +(adds...)
end
function evalTerm(sum::SymbolicUtils.Sym{Parameter,IndexedAverageDoubleSum};mapping::Dict{Symbol,Int64})
    return evalTerm(IndexedAverageDoubleSum(evalTerm(sum.metadata.innerSum;mapping),sum.metadata.sumIndex,sum.metadata.nonEqualIndices);mapping)
end
function evalTerm(term::SymbolicUtils.Mul;mapping::Dict{Symbol,Int64}=Dict{Symbol,Int64}()) 
    mults = []
    for arg in arguments(term)
        push!(mults,evalTerm(arg;mapping))
    end
    if isempty(mults)
        return 0
    elseif length(mults) == 1
        return mults[1]
    end
    return *(mults...)
end
function evalTerm(term::SymbolicUtils.Add;mapping::Dict{Symbol,Int64}=Dict{Symbol,Int64}()) 
    adds = []
    for arg in arguments(term)
        push!(adds,evalTerm(arg;mapping))
    end
    if isempty(adds)
        return 0
    elseif length(adds) == 1
        return adds[1]
    end
    return +(adds...)
end
evalTerm(x;mapping::Dict{Symbol,Int64}) = x
evalEq(eq::Symbolics.Equation;mapping::Dict{Symbol,Int64}=Dict{Symbol,Int64}()) = Symbolics.Equation(eq.lhs,evalTerm(eq.rhs;mapping))

function getLHS(eq::Symbolics.Equation)
    return eq.lhs
end
getLHS(x) = []

#functions to order terms inside the sumy by their index-number, used for checking if averages already exist in LHS of the equations
function orderTermsByNumber(qmul::QMul)
    args_nc = qmul.args_nc
    newargs = sort(args_nc,by=getNumber)
    return QMul(1,newargs)
end
function orderTermsByNumber(term1::Term{AvgSym, Nothing})
    if typeof(arguments(term1)[1]) <: QMul
        return average(orderTermsByNumber(arguments(term1)[1]))
    end
    return term1
end
function orderTermsByNumber(mul::SymbolicUtils.Mul)
    args = arguments(mul)
    args_ = []
    for arg in args
        push!(args_, orderTermsByNumber(arg))
    end
    return *(args_...)
end
function orderTermsByNumber(add::SymbolicUtils.Add)
    args = arguments(add)
    args_ = []
    for arg in args
        push!(args_, orderTermsByNumber(arg))
    end
    return +(args_...)
end
orderTermsByNumber(eq::Symbolics.Equation) = Symbolics.Equation(eq.lhs,orderTermsByNumber(eq.rhs))
orderTermsByNumber(x) = x
getNumber(op::NumberedOperator) = op.numb
getNumber(x) = 0

Base.:(==)(term1::Term{AvgSym, Nothing},term2::Term{AvgSym, Nothing}) = isequal(arguments(term1), arguments(term2))

#Insert values no longer needed to do explicitly, now integrated, when defining an ODEProblem
function insertValuesIntoMeanfieldEquation(me::MeanfieldEquations,valMapping::Dict{Sym{Parameter, numberedVariable},SNuN})
    newME = deepcopy(me)
    newEqs = Vector{Union{Missing,Symbolics.Equation}}(missing,length(me.equations))
    newOpEqs = []
    Threads.@threads for i=1:length(me.equations)
        newEqs[i] = insertValuesIntoEquation(me.equations[i],valMapping)
    end
    newME = MeanfieldEquations(newEqs,me.operator_equations,me.states,me.operators,me.hamiltonian,me.jumps,me.jumps_dagger,me.rates,me.iv,me.varmap,me.order)
    return newME
end
function insertValuesIntoEquation(eq::Symbolics.Equation,valMapping::Dict{Sym{Parameter, numberedVariable},SNuN})
    return Symbolics.Equation(eq.lhs,insertValues(eq.rhs,valMapping))
end
function insertValues(term::SymbolicUtils.Mul,valMapping::Dict{Sym{Parameter, numberedVariable},SNuN})
    mults = []
    for arg in arguments(term)
        push!(mults,insertValues(arg,valMapping))
    end
    if length(mults) == 1
        return mults[1]
    end
    return *(mults...)
end
function insertValues(term::SymbolicUtils.Add,valMapping::Dict{Sym{Parameter, numberedVariable},SNuN})
    adds = []
    for arg in arguments(term)
        push!(adds, insertValues(arg,valMapping))
    end
    if length(adds) == 1
        return adds[1]
    end
    return +(adds...)
end
function insertValues(val::SymbolicUtils.Sym{Parameter,numberedVariable},valMapping::Dict{Sym{Parameter, numberedVariable},SNuN})
    if length(valMapping) == 1
        if val.metadata.name == first(keys(valMapping)).metadata.name
            return first(values(valMapping))
        end
    end
    if val in keys(valMapping)
        return valMapping[val]
    end
    return val
end
insertValues(x,valMapping) = x

#Value map creation, for easier inserting into the ODEProblem
function createValueMap(sym::Sym{Parameter, IndexedVariable}, values::Vector)
    iVar = sym.metadata
    if iVar.ind.rangeN != length(values)
        error("different length of index-range and given values!")
    end
    dict = Dict{Sym{Parameter, Base.ImmutableDict{DataType, Any}},Float64}()
    for i = 1:iVar.ind.rangeN
        push!(dict,(SingleNumberedVariable(iVar.name,i) => values[i]))
    end
    return dict
end
function createValueMap(sym::Sym{Parameter, IndexedVariable}, value::Number)
    iVar = sym.metadata
    dict = Dict{Sym{Parameter, Base.ImmutableDict{DataType, Any}},Float64}()
    push!(dict,(SingleNumberedVariable(iVar.name,1) => value))
    return dict
end
function createValueMap(sym::Sym{Parameter,DoubleIndexedVariable},values::Matrix)
    dict = Dict{Sym{Parameter, Base.ImmutableDict{DataType, Any}},Float64}()
    var = sym.metadata
    for i = 1:var.ind1.rangeN
        for j = 1:var.ind2.rangeN
            push!(dict,(DoubleNumberedVariable(var.name,i,j) => values[i,j]))
        end
    end
    return dict
end

#functions for checking if indices occure in specific terms
function containsIndexedOps(term::SymbolicUtils.Term{AvgSym, Nothing})
    arg_ = arguments(term)
    found = false
    if typeof(arg_[1]) <: QMul
        for arg in arg_[1].args_nc
            if typeof(arg) == IndexedOperator
                found = true
            end
        end
    else
        return typeof(arg_[1]) == IndexedOperator 
    end
    return found
end
function getIndices(term::SymbolicUtils.Term{AvgSym, Nothing})
    arg_ = arguments(term)
    indices = []
    if typeof(arg_[1]) <: QMul
        for arg in arg_[1].args_nc
            if typeof(arg) == IndexedOperator && arg.ind ∉ indices
                push!(indices, arg.ind)
            end
        end
    else
        if typeof(arg_[1]) == IndexedOperator 
            return [arg_[1].ind]
        end
    end
    return indices
end
function getIndices(term::QMul)
    indices = []
    for arg in term.args_nc
        if typeof(arg) == IndexedOperator && arg.ind ∉ indices
            push!(indices,arg.ind)
        end
    end
    return indices
end
getIndices(a::QNumber) = typeof(a) == IndexedOperator ? [a.ind] : []


#simplify functions not "really" needed, they are nice to have, since equation creation of Symbolics sometimes does not simplify certain terms
#function to reduce multiplication of numbers with a sum into just a sum of multiplication
function simplifyMultiplication(term::SymbolicUtils.Mul)
    args = arguments(term)
    sum = nothing
    sumArg = 0
    found = false
    #check where sum is, if sum is there
    for i = 1:length(args)
        if typeof(args[i]) <: SymbolicUtils.Add
            found = true
            sumArg = i
            break
        end
    end 
    found || return term #no add-terms were found inside the multiplication
    adds = []
    lefts = []
    rights = []
    #split multiplication into left and right terms
    for i = 1:length(args)
        if i < sumArg
            push!(lefts,args[i])
        elseif i > sumArg
            push!(rights,args[i])
        end
    end
    for arg in arguments(args[sumArg]) #arguments in the sum
        push!(adds, simplifyMultiplication(*(lefts...) * arg *(rights...)))
    end
    return +(adds...)
end
function simplifyAdd(term::SymbolicUtils.Add)
    adds = []
    for arg in arguments(term)
        if typeof(arg) <: SymbolicUtils.Mul
            push!(adds,simplifyMultiplication(arg))
        else
            push!(adds,arg)
        end
    end
    return +(adds...)
end
function simplifyEquation(eq::Symbolics.Equation)
    if typeof(eq.rhs) <: SymbolicUtils.Add
        return Symbolics.Equation(eq.lhs,simplifyAdd(eq.rhs))
    elseif typeof(eq.rhs) <: SymbolicUtils.Mul
        return Symbolics.Equation(eq.lhs,simplifyMul(eq.rhs))
    else
        return eq
    end
end
function simplifyMeanfieldEquations(me::MeanfieldEquations)
    eq_after = []
    op_after = []
    for eq in me.equations
        push!(eq_after,simplifyEquation(eq))
    end
    for eq in eq_after
        push!(op_after,undo_average(eq))
    end
    return MeanfieldEquations(eq_after,op_after,me.states,me.operators,me.hamiltonian,me.jumps,me.jumps_dagger,me.rates,me.iv,me.varmap,me.order)
end 

#functions for simplifying the indexedComplete function
function getOps(sum::SymbolicUtils.Sym{Parameter,IndexedAverageSum})
    term = sum.metadata.term
    if typeof(term) == Term{AvgSym, Nothing}
        return Any[arguments(term)[1].args_nc]
    end
    if typeof(term) <: SymbolicUtils.Mul
        ops = Any[]
        for arg in arguments(term)
            if typeof(arg) == Term{AvgSym, Nothing}
                push!(ops,arguments(arg)[1].args_nc)
            end
        end
        return ops
    end
end
function getAvrgs(sum::SymbolicUtils.Sym{Parameter,IndexedAverageSum})
    term = sum.metadata.term
    if typeof(term) == Term{AvgSym, Nothing}
        return Any[term]
    end
    if typeof(term) <: SymbolicUtils.Mul
        ops = Any[]
        for arg in arguments(term)
            if typeof(arg) == Term{AvgSym, Nothing}
                push!(ops,arg)
            end
        end
        return ops
    end
end

#not sure if needed, I think not
function evaluateTerm(term::SymbolicUtils.Add, indMap::Dict{Index,Int64})
    addterms = []
    for arg in arguments(term)
        push!(addterms,evaluateTerm(arg, indMap))
    end
    return +(addterms...)
end

function Base.show(io::IO,indSum::IndexedAverageSum) 
    write(io, "Σ", "($(indSum.sumIndex.name)", "=1:$(indSum.sumIndex.rangeN))",)
    if !(isempty(indSum.nonEqualIndices))
        write(io,"($(indSum.sumIndex.name) ≠ ")
        for i = 1:length(indSum.nonEqualIndices)
            write(io, "$(indSum.nonEqualIndices[i].name)")
            if i == length(indSum.nonEqualIndices)
                write(io,")")
            else
                write(io,",")
            end
        end
    end
    Base.show(io,indSum.term)
end
function Base.show(io::IO,indSum::IndexedAverageDoubleSum)
    write(io, "Σ", "($(indSum.sumIndex.name)", "=1:$(indSum.sumIndex.rangeN))",)
    if !(isempty(indSum.nonEqualIndices))
        write(io,"($(indSum.sumIndex.name) ≠ ")
        for i = 1:length(indSum.nonEqualIndices)
            write(io, "$(indSum.nonEqualIndices[i].name)")
            if i == length(indSum.nonEqualIndices)
                write(io,")")
            else
                write(io,",")
            end
        end
    end
    Base.show(io,indSum.innerSum)
end
function Base.show(io::IO, numbOp::NumberedOperator)
    Base.show(io,numbOp.op)
    Base.show(io,numbOp.numb)
end