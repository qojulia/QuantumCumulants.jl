#Main file for manipulating indexed averages and sums over averages.
using ModelingToolkit

struct EvaledMeanfieldEquations <: AbstractMeanfieldEquations
    equations::Vector{Symbolics.Equation}
    operator_equations::Vector{Symbolics.Equation}
    states::Vector
    operators::Vector{QNumber}
    hamiltonian::QNumber
    jumps::Vector
    jumps_dagger
    rates::Vector
    iv::SymbolicUtils.Sym
    varmap::Vector{Pair}
    order::Union{Int,Vector{<:Int},Nothing}
end
function Base.show(io::IO,de::EvaledMeanfieldEquations)
    write(io,"Evaluated Meanfield equations with: ")
    write(io, "$(length(de.equations))")
    write(io, " number of equations")
end
@latexrecipe function f(de::EvaledMeanfieldEquations)
    return de
end
function plotME(me::EvaledMeanfieldEquations)
    return MeanfieldEquations(me.equations,me.operator_equations,me.states,me.operators,me.hamiltonian,me.jumps,me.jumps_dagger,me.rates,me.iv,me.varmap,me.order)
end


const symbolics_terms = Union{<:Average,<:SymbolicUtils.Mul,<:SymbolicUtils.Sym}
"""
    numberedVariable <: CNumber

abstract type for numbered Variables.
"""
abstract type numberedVariable <: CNumber end
"""

    IndexedAverageSum <: CNumber

Defines a symbolic summation over an average, or a multiplication of several averages, using one [`Index`](@ref) entity.

Fields:
======

* term: A multiplication of average terms.
* sumIndex: The index, for which the summation will go over.
* nonEqualIndices: (optional) A vector of indices, for which the summation-index can not be equal with.

"""
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
        if sumIndex ∉ getIndices(term)
            return (sumIndex.rangeN - length(nonEqualIndices)) * term
        end
        neis_sym = ""
        if !(isempty(nonEqualIndices))
            neis_sym = string("(",neis_sym)
            neis_sym = string(neis_sym, "$(sumIndex.name)≠")
            neis_sym = string(neis_sym, writeNEIs(nonEqualIndices))
            neis_sym = string(neis_sym,")")
        end
        prefact = 1.0
        if term isa SymbolicUtils.Mul && arguments(term)[1] isa Number # put numbers outside of sum (for easier evaluation)
            prefact = arguments(term)[1]
            args_nc = arguments(term)[2:end]
            if length(args_nc) == 1
                term = args_nc[1]
            else
                term = *(args_nc...)
            end
        end
        metadata = new(term,sumIndex,nonEqualIndices)
        return prefact*SymbolicUtils.Sym{Parameter, IndexedAverageSum}(Symbol("∑($(sumIndex.name)=1:$(sumIndex.rangeN))$(neis_sym)$(term)"), metadata) #Symbol("∑($(sumIndex.name)=1:$(sumIndex.rangeN))$(neis_sym)$(term)")
    end
end

"""

    IndexedAverageDoubleSum <: CNumber

Defines a symbolic summation over an [`IndexedAverageSum`](@ref), using a [`Index`](@ref) entity. This schematically represent a double-sum over a multiplication of averages.

Fields:
======

* innerSum: An [`IndexedAverageSum`](@ref) entity.
* sumIndex: The index, for which the (outer) summation will go over.
* nonEqualIndices: (optional) A vector of indices, for which the (outer) summation-index can not be equal with.

"""
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
"""

    NumberedOperator <: QNumber

Defines an operator, associated with a Number. Commutator-relations are calculated using these numbers, as a sort of an index.

Fields:
======

* op: An Operator, either a [`Transition`](@ref), a [`Destroy`](@ref) or a [`Create`](@ref) can be defined.
* numb: An Integer Number.

"""
struct NumberedOperator <:QNumber
    op::indexable
    numb::Int64
    function NumberedOperator(op,numb)
        if numb <= 0
            error("can not create numbered-operator with negative or 0 number")
            return 0
        end
        if typeof(op) <: SNuN
            return op
        end
        if SymbolicUtils.istree(op)
            f = SymbolicUtils.operation(op)
            args = []
            for arg in SymbolicUtils.arguments(op)
                push!(args,NumberedOperator(arg,numb))
            end
            if isempty(args)
                return 0
            elseif length(args) == 1
                return args[1]
            end
            return f(args...)
        end
        return new(op,numb)
    end
end
#Variables
"""

    SingleNumberedVariable <: numberedVariable

Defines an variable, associated with a Number. Used by [`value_map`](@ref)

Fields:
======

* name: The name of the variable.
* numb: An Integer Number.

"""
struct SingleNumberedVariable <: numberedVariable
    name::Symbol
    numb::Int64
    function SingleNumberedVariable(name,numb)
        metadata=source_metadata(:Parameter, name)
        s = SymbolicUtils.Sym{Parameter, typeof(metadata)}(Symbol("$(name)_$(numb)"), metadata)
        return SymbolicUtils.setmetadata(s, MTK.MTKParameterCtx, true)
    end
end
"""

    DoubleNumberedVariable <: numberedVariable

Defines an variable, associated with two Numbers. Used by [`value_map`](@ref)

Fields:
======

* name: The name of the variable.
* numb1: An Integer Number.
* numb2: Another Integer Number.

"""
struct DoubleNumberedVariable <: numberedVariable
    name::Symbol
    numb1::indornum
    numb2::indornum
    function DoubleNumberedVariable(name,numb1,numb2)
        if typeof(numb1) == typeof(numb2) && typeof(numb1) == Int64
            metadata = source_metadata(:Parameter, name)
            s = SymbolicUtils.Sym{Parameter, typeof(metadata)}(Symbol("$(name)_{$(numb1)$(numb2)}"), metadata)
            return SymbolicUtils.setmetadata(s, MTK.MTKParameterCtx, true)
        else
            return SymbolicUtils.Sym{Parameter, numberedVariable}(Symbol("$(name)_{$(numb1)$(numb2)}"), new(name,numb1,numb2))
        end
    end
end
#Special averages
struct NumberedIndexedAverage <: CNumber
    term::symbolics_terms
    indexMapping::Vector{Tuple{Index,Int64}} #not to be confused with later definition of index mapping here if (l,2) means l ≠ 2
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
            return sum(adds)
        end
        if typeof(term) == SymbolicUtils.Term{AvgSym,Nothing}
            if SymbolicUtils._iszero(arguments(term)[1])
                return 0
            end
            metadata = new(term,indexMapping)
            neis = writeIndexNEIs(indexMapping)
            return SymbolicUtils.Sym{Parameter, SpecialIndexedAverage}(Symbol("$(neis)$(term)"), metadata)
        end
        if typeof(term) <: SymbolicUtils.Mul
            prefac = 1.0
            term_ = []
            isempty(arguments(term)) && return 0
            for arg in arguments(term)
                if typeof(arg) <: Number
                    prefac = prefac * arg
                else
                    push!(term_,arg)
                end
            end
            isempty(term_) && return prefac

            term = length(term_) == 1 ? term_[1] : *(term_...)
            if typeof(term) <: symbolics_terms
                if typeof(term) == SymbolicUtils.Term{AvgSym,Nothing}
                    if SymbolicUtils._iszero(arguments(term)[1])
                        return 0
                    end
                end
                metadata = new(term,indexMapping)
                neis = writeIndexNEIs(indexMapping)
                return prefac * SymbolicUtils.Sym{Parameter, SpecialIndexedAverage}(Symbol("$(neis)$(term)"), metadata)
            else
                return prefac * SpecialIndexedAverage(term,indexMapping)
            end
        end
        return term
    end
end

const AvgSums = Union{SymbolicUtils.Sym{Parameter,IndexedAverageSum},SymbolicUtils.Sym{Parameter,IndexedAverageDoubleSum},SymbolicUtils.Sym{Parameter,SpecialIndexedAverage},IndexedAverageSum,IndexedAverageDoubleSum,SpecialIndexedTerm}

function average(indOp::IndexedOperator) 
    if SymbolicUtils._iszero(indOp)
        return 0
    end
    return _average(indOp)
end
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
undo_average(a::SymbolicUtils.Sym{Parameter,SpecialIndexedAverage}) = reorder(undo_average(a.metadata.term),a.metadata.indexMapping)

#define calculus for numbered operators -> break it down into QNuber multiplication
*(numOp::NumberedOperator, qmul::QMul) = merge_commutators(qmul.arg_c,vcat(numOp,qmul.args_nc))
*(qmul::QMul, numOp::NumberedOperator) = merge_commutators(qmul.arg_c,vcat(qmul.args_nc,numOp))
function *(numOp1::NumberedOperator,numOp2::NumberedOperator) 
    if numOp1.op isa Create || numOp1.op isa Destroy || numOp2.op isa Create || numOp2.op isa Destroy
        return merge_commutators(1,[numOp1,numOp2])
    end
    return numOp1.numb == numOp2.numb ? NumberedOperator(numOp1.op*numOp2.op,numOp1.numb) : QMul(1,[numOp1,numOp2])
end
*(elem::SNuN, numOp::NumberedOperator) = merge_commutators(elem,[numOp])
*(numOp::NumberedOperator,elem::SNuN) = merge_commutators(elem,[numOp])
*(a::Create,b::NumberedOperator) = merge_commutators(1,[a,b])
*(b::NumberedOperator,a::Create) = merge_commutators(1,[b,a])
*(a::Destroy,b::NumberedOperator) = merge_commutators(1,[a,b])
*(b::NumberedOperator,a::Destroy) = merge_commutators(1,[b,a])
*(a::IndexedOperator,b::NumberedOperator) = merge_commutators(1,[a,b])
*(b::NumberedOperator,a::IndexedOperator) = merge_commutators(1,[b,a])

#Symbolics functions
get_order(a::Sym{Parameter,IndexedAverageSum}) = get_order(a.metadata.term)
cumulant_expansion(a::IndexedAverageSum,order::Int) = IndexedAverageSum(simplifyMultiplifcation(cumulant_expansion(a.term,order)),a.sumIndex,a.nonEqualIndices) #not used (?)
SymbolicUtils.istree(a::IndexedAverageSum) = false
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
acts_on(x::NumberedOperator) = acts_on(x.op)
get_order(x::NumberedOperator) = get_order(x.op)

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
Base.isequal(a::SymbolicUtils.Sym{Parameter,IndexedAverageSum},b::SymbolicUtils.Sym{Parameter,IndexedAverageSum}) = isequal(a.metadata,b.metadata)
function Base.isequal(a::IndexedAverageSum, b::IndexedAverageSum)
    isequal(a.sumIndex,b.sumIndex) || return false
    isequal(a.term, b.term) || return false
    isequal(a.nonEqualIndices,b.nonEqualIndices) || return false
    return true
end
Base.isequal(a::SymbolicUtils.Sym{Parameter,IndexedAverageSum},x) = false
Base.isequal(a::IndexedAverageSum,b) = false
Base.isequal(::Sym{Parameter, IndexedAverageSum}, ::Sym) = false
Base.isequal(::IndexedAverageSum, ::SymbolicUtils.Symbolic) = false
Base.isequal(a::SymbolicUtils.Sym{Parameter,SpecialIndexedAverage},b::SymbolicUtils.Sym{Parameter,SpecialIndexedAverage}) = isequal(a.metadata.term,b.metadata.term) && isequal(a.metadata.indexMapping,b.metadata.indexMapping)
function Base.isequal(nVal1::Sym{Parameter,numberedVariable},nVal2::Sym{Parameter,numberedVariable}) 
    if typeof(nVal1) == typeof(nVal2) && typeof(nVal1.metadata) == SingleNumberedVariable
        return (nVal1.metadata.name == nVal2.metadata.name) && (nVal1.metadata.numb == nVal2.metadata.numb)
    elseif typeof(nVal1) == typeof(nVal2) && typeof(nVal1.metadata) == DoubleNumberedVariable
        return (nVal1.metadata.name == nVal2.metadata.name) && (nVal1.metadata.numb1 == nVal2.metadata.numb1) && (nVal1.metadata.numb2 == nVal2.metadata.numb2)
    end
    return false
end
Base.isequal(::SymbolicUtils.Sym{Parameter,IndexedAverageSum},::SymbolicUtils.Sym{Parameter,IndexedAverageDoubleSum}) = false
Base.isequal(::Sym{Parameter, IndexedAverageSum}, ::SymbolicUtils.Symbolic) = false
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
SymbolicUtils.arguments(op::IndexedAverageSum) = arguments(op.term)
SymbolicUtils.arguments(op::Sym{Parameter, IndexedAverageDoubleSum}) = arguments(op.metadata)
SymbolicUtils.arguments(op::IndexedAverageDoubleSum) = op.innerSum
SymbolicUtils.arguments(op::Sym{Parameter,SpecialIndexedAverage}) = arguments(op.metadata)
SymbolicUtils.arguments(op::SpecialIndexedAverage) = arguments(op.term)

#this is the new method, insert values directly into the average before calculating anything, simplifies evaluation afterwards extremely
#function for inserting index, k -> 1,2,...,N
"""
    inserIndex(term,ind,value)

Function, that inserts an integer value for a index in a specified term. 
This function creates Numbered- Variables/Operators/Sums upon calls.

Examples
========

    insert_index(σⱼ²¹,j,1) = σ₁²¹ 

"""
function insert_index(sum::SymbolicUtils.Sym{Parameter,IndexedAverageSum}, ind::Index, value::Int64)
    if ind == sum.metadata.sumIndex
        error("cannot exchange summation index with number!")
    end
    if ind in sum.metadata.nonEqualIndices
        newNEI = filter(x-> x != ind,sum.metadata.nonEqualIndices)
        push!(newNEI,value)
        return IndexedAverageSum(insert_index(sum.metadata.term,ind,value),sum.metadata.sumIndex,newNEI)
    else
        return IndexedAverageSum(insert_index(sum.metadata.term,ind,value),sum.metadata.sumIndex,sum.metadata.nonEqualIndices)
    end
end
function insert_index(sum::SymbolicUtils.Sym{Parameter,IndexedAverageDoubleSum}, ind::Index,value::Int64)
    inner = insert_index(sum.metadata.innerSum,ind,value)
    return IndexedAverageDoubleSum(inner,sum.metadata.sumIndex,sum.metadata.nonEqualIndices)
end
function insert_index(term::SymbolicUtils.Mul, ind::Index, value::Int64)
    args = []
    for arg in arguments(term)
        push!(args,insert_index(arg,ind,value))
    end
    return *(args...)
end
function insert_index(term::SymbolicUtils.Add,ind::Index,value::Int64)
    args = []
    for arg in arguments(term)
        push!(args,insert_index(arg,ind,value))
    end
    return sum(args)
end
function insert_index(term::SymbolicUtils.Pow,ind::Index,value::Int64)
    return insert_index(arguments(term)[1],ind,value)^(arguments(term)[2])
end
function insert_index(term::SymbolicUtils.Term{AvgSym,Nothing},ind::Index,value::Int64)
    newargs = []
    if typeof(arguments(term)[1]) <: QMul
        for arg in arguments(term)[1].args_nc
            push!(newargs,insert_index(arg,ind,value))
        end
        if isempty(newargs)
            return 0
        end
        if length(newargs) == 1
            qmul = newargs[1]
        else
            qmul = *(newargs...)
        end
        return average(qmul)
    else
        return average(insert_index(arguments(term)[1],ind,value))
    end
end
function insert_index(term_::Sym{Parameter, DoubleIndexedVariable},ind::Index,value::Int64)
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
function insert_index(term::SymbolicUtils.Sym{Parameter,numberedVariable},ind::Index,value::Int64)
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
function insert_index(term::SymbolicUtils.Sym{Parameter,SpecialIndexedAverage},ind::Index,value::Int64)
    meta = term.metadata
    newterm = insert_index(meta.term,ind,value)
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
function insert_index(qmul::QMul,ind::Index,to::Int64)
    args_after = []
    for arg in qmul.args_nc
        push!(args_after,insert_index(arg,ind,to))
    end
    sort!(args_after, by=getNumber)
    return *(qmul.arg_c,args_after...)
end
insert_index(eq::Symbolics.Equation,ind::Index,value::Int64) = Symbolics.Equation(insert_index(eq.lhs,ind,value),insert_index(eq.rhs,ind,value))
insert_index(term::IndexedOperator,ind::Index,value::Int64) = term.ind == ind ? NumberedOperator(term.op,value) : term
insert_index(term::SymbolicUtils.Sym{Parameter,IndexedVariable},ind::Index,value::Int64) = term.metadata.ind == ind ? SingleNumberedVariable(term.metadata.name,value) : term
insert_index(x,ind::Index,value::Int64) = x
"""
    insertIndices(eq::Symbolics.Equation,map::Dict{Index,Int64};mapping::Dict{SymbolicUtils.Sym,Int64}=Dict{Symbol,Int64}())

Function, that inserts an integer value for a index in a specified Equation. This function creates Numbered- Variables/Operators/Sums upon calls.
Mainly used by [`evalEquation`](@ref).

# Arguments
*`eq::Symbolics.Equation`: The equation for which the indices will be inserted.
*`map::Dict{Index,Int64}`: A dictionary, which contains specifications for the insertions
    the entry (i => 5) would result in all `i` indices being replaced with the number 5.

# Optional argumentes
*`mapping::Dict{SymbolicUtils.Sym,Int64}=Dict{Symbol,Int64}()`: A seperate dictionary, to
    specify any symbolic limits used when [`Index`](@ref) entities were defined. This needs
    to be specified, when the equation contains summations, for which the upper bound is given
    by a Symbolic.

"""
function insertIndices(eq::Symbolics.Equation,map::Dict{Index,Int64};mapping::Dict{SymbolicUtils.Sym,Int64}=Dict{Symbol,Int64}())
    eq_ = eq
    while !isempty(map)
        pair = first(map)
        eq_ = insert_index(eq_,first(pair),last(pair))
        delete!(map,first(pair))
    end
    return evalEq(eq_;mapping) #return finished equation
end
function evalEquation(eq::Symbolics.Equation,arr,indices;mapping::Dict{SymbolicUtils.Sym,Int64}=Dict{Symbol,Int64}())
    if !(isempty(indices))
        eqs = Vector{Any}(nothing, length(arr))
        #Threads.@threads 
        for i=1:length(arr)
            dict = Dict(indices .=> arr[i])
            eq_ = orderTermsByNumber(insertIndices(eq,dict;mapping))
            eqs[i] = eq_
        end
        return filter(x -> !=(x,nothing),eqs)
    else
        return [evalEq(eq;mapping)]
    end
end
"""
    evalME(me::MeanfieldEquations;mapping::Dict{SymbolicUtils.Sym,Int64}=Dict{SymbolicUtils.Sym,Int64}())

Function, that evaluates a given [`MeanfieldEquations`](@ref) entity and returns again equations,
where indices have been inserted and sums evaluated.

# Arguments
*`me::MeanfieldEquations`: A [`MeanfieldEquations`](@ref) entity, which shall be evaluated.

# Optional argumentes
*`mapping::Dict{SymbolicUtils.Sym,Int64}=Dict{Symbol,Int64}()`: A seperate dictionary, to
    specify any symbolic limits used when [`Index`](@ref) entities were defined. This needs
    to be specified, when the equations contain summations, for which the upper bound is given
    by a Symbolic.

"""
function evalME(me::AbstractMeanfieldEquations;mapping::Dict{SymbolicUtils.Sym,Int64}=Dict{SymbolicUtils.Sym,Int64}())#this is still pretty slow
    indices = nothing
    for eq in me.equations
        if containsIndexedOps(eq.lhs) && length(getIndices(eq.lhs)) == me.order
            indices = getIndices(eq.lhs)
            break
        end
    end
    if indices === nothing
        indices = []
        indVecs = getIndices.(me.states)
        for indvec in indVecs
            for ind in indvec
                if ind ∉ indices
                    push!(indices,ind)
                end
            end
        end
    end

    range = 0
    for ind in indices
        if typeof(ind.rangeN) != Int64 && ind.rangeN ∉ keys(mapping)
            error("Please provide numbers for the upper-limits used: $(ind.rangeN); you can do this by calling: evaluate(me;mapping=[Dictionary with corresponding numbers for limits])")
        end
    end
    
    sort!(indices)

    ranges = []
    arrays = []
    for ind in indices
        if ind.rangeN in keys(mapping)
            push!(ranges,1:mapping[ind.rangeN])
        else
            push!(ranges,1:ind.rangeN)
        end
    end
    maxRange = maximum(ranges)[end]
    newEqs = Vector{Union{Missing,Symbolics.Equation}}(missing,length(indices)*(maxRange*5)^length(indices))
    
    counter = 1
    for i = 1:length(me.equations)
        inds = getIndices(me.equations[i].lhs)
        if isempty(inds)
            evals = evalEquation(me.equations[i],[],[];mapping)
            for eq_ in evals
                if (eq_.lhs ∉ getLHS.(newEqs)) && (_conj(eq_.lhs) ∉ getLHS.(newEqs))
                    newEqs[counter] = eq_
                    counter = counter+1
                end
            end
            continue
        end

        nLvls=[]
        others=[]
        for ind in inds
            if ind.specHilb isa NLevelSpace
                push!(nLvls,ind)
            else
                push!(others,ind)
            end
        end
        ranges_nLvl = []
        for ind in nLvls
            if ind.rangeN in keys(mapping)
                push!(ranges_nLvl,1:mapping[ind.rangeN])
            else
                push!(ranges_nLvl,1:ind.rangeN)
            end
        end
        ranges_other = []
        for ind in others
            if ind.rangeN in keys(mapping)
                push!(ranges_other,1:mapping[ind.rangeN])
            else
                push!(ranges_other,1:ind.rangeN)
            end
        end

        arr1 = unique(sort.(collect.(filter(x -> length(x) == length(unique(x)),collect(Iterators.product(ranges_nLvl...))))))
        arr2 = unique(sort.(collect.(collect(Iterators.product(ranges_other...)))))

        arrf = appendEach(arr1,arr2)

        evals = evalEquation(me.equations[i],arrf,[nLvls;others];mapping)
        for eq_ in evals
            if (eq_.lhs ∉ getLHS.(newEqs)) && (_conj(eq_.lhs) ∉ getLHS.(newEqs))
                newEqs[counter] = eq_
                counter = counter+1
            end
        end
                    
    end
    
    newEqs = filter(x -> !isequal(x,missing), newEqs)
    vs = getLHS.(newEqs)
    varmap = make_varmap(vs, me.iv)
    return EvaledMeanfieldEquations(newEqs,me.operator_equations,vs,me.operators,me.hamiltonian,me.jumps,me.jumps_dagger,me.rates,me.iv,varmap,me.order)
 end
function appendEach(vec1,vec2)
    isempty(vec1) && return vec2
    isempty(vec2) && return vec1
    vecf = Vector{Vector{Int64}}(undef,length(vec1)*length(vec2))
    counter = 1
    for v1 in vec1
        for v2 in vec2
            vecf[counter] = [v1;v2]
            counter=counter+1
        end
    end
    return vecf
end
function eval_term(sum_::SymbolicUtils.Sym{Parameter,IndexedAverageSum};mapping::Dict{SymbolicUtils.Sym,Int64}=Dict{SymbolicUtils.Sym,Int64}())
    rangeEval = 0
    if sum_.metadata.sumIndex.rangeN in keys(mapping)
        rangeEval = mapping[sum_.metadata.sumIndex.rangeN]
    else
        if typeof(sum_.metadata.sumIndex.rangeN) <: SymbolicUtils.Mul
            args = arguments(sum_.metadata.sumIndex.rangeN)
            for i=1:length(args)
                if args[i] in keys(mapping)
                    args[i] = mapping[args[i]]
                end
            end
            rangeEval = *(args...)
        else
            rangeEval = sum_.metadata.sumIndex.rangeN
        end
    end
    adds = Vector{Any}(nothing,rangeEval)
    for i = 1:rangeEval
        if i in sum_.metadata.nonEqualIndices
            continue
        end
        adds[i]=orderTermsByNumber(insert_index(sum_.metadata.term,sum_.metadata.sumIndex,i))
    end
    filter!(x->x!=nothing,adds)
    if isempty(adds)
        return 0
    elseif length(adds) == 1
        return adds[1]
    end
    return sum(adds)
end
function eval_term(sum::SymbolicUtils.Sym{Parameter,IndexedAverageDoubleSum};mapping::Dict{SymbolicUtils.Sym,Int64})
    return eval_term(IndexedAverageDoubleSum(eval_term(sum.metadata.innerSum;mapping),sum.metadata.sumIndex,sum.metadata.nonEqualIndices);mapping)
end
function eval_term(term::SymbolicUtils.Mul;mapping::Dict{SymbolicUtils.Sym,Int64}=Dict{SymbolicUtils.Sym,Int64}()) 
    mults = []
    for arg in arguments(term)
        push!(mults,eval_term(arg;mapping))
    end
    if isempty(mults)
        return 0
    elseif length(mults) == 1
        return mults[1]
    end
    return *(mults...)
end
function eval_term(term::SymbolicUtils.Add;mapping::Dict{SymbolicUtils.Sym,Int64}=Dict{SymbolicUtils.Sym,Int64}()) 
    adds = []
    for arg in arguments(term)
        push!(adds,eval_term(arg;mapping))
    end
    if isempty(adds)
        return 0
    elseif length(adds) == 1
        return adds[1]
    end
    return sum(adds)
end
eval_term(x;mapping::Dict{SymbolicUtils.Sym,Int64}) = x
evalEq(eq::Symbolics.Equation;mapping::Dict{SymbolicUtils.Sym,Int64}=Dict{SymbolicUtils.Sym,Int64}()) = Symbolics.Equation(eq.lhs,eval_term(eq.rhs;mapping))

getLHS(eq::Symbolics.Equation) = eq.lhs
getLHS(x) = []

#functions to order terms inside the sumy by their index-number, used for checking if averages already exist in LHS of the equations
function orderTermsByNumber(qmul::QMul)
    args_nc = qmul.args_nc
    newargs = sort(args_nc,by=getNumber)
    if isempty(newargs)
        return qmul.arg_c
    end
    return *(qmul.arg_c,newargs...)
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
    return sum(args_)
end
orderTermsByNumber(eq::Symbolics.Equation) = Symbolics.Equation(eq.lhs,orderTermsByNumber(eq.rhs))
orderTermsByNumber(x) = x
Base.:(==)(term1::Term{AvgSym, Nothing},term2::Term{AvgSym, Nothing}) = isequal(arguments(term1), arguments(term2))

#Value map creation, for easier inserting into the ODEProblem
"""
    create_value_map(sym::Sym{Parameter, IndexedVariable}, values::Vector;mapping::Dict{SymbolicUtils.Sym,Int64}=Dict{SymbolicUtils.Sym,Int64}())
    create_value_map(sym::Sym{Parameter, IndexedVariable}, value::Number)
    create_value_map(sym::Sym{Parameter,DoubleIndexedVariable},values::Matrix;mapping::Dict{SymbolicUtils.Sym,Int64}=Dict{SymbolicUtils.Sym,Int64}())

Function, that creates a Dictionary, for which a indexedVariable is associated with a series of (number) values. The dictionary contains Symbols of either [`SingleNumberedVariable`](@ref)
or [`DoubleNumberedVariables`](@ref) as keys and the values as values. For a Single-indexed variable, one can
create such a mapping by giving a Vector of values, and for double-indexed variables by giving a Matrix. One can also create such a mapping, by using
only a single value, then all possible numbered-Variables are set to the same values.

# Arguments
*`sym`: Either a [`IndexedVariable`](@ref) or a [`DoubleIndexedVariable`](@ref)
*`values`: For a [`IndexedVariable`](@ref) either a vector or a single number, and for [`DoubleIndexedVariable`](@ref) a matrix.

# Optional argumentes
*`mapping::Dict{SymbolicUtils.Sym,Int64}=Dict{Symbol,Int64}()`: A seperate dictionary, to
    specify any symbolic limits used when [`Index`](@ref) entities were defined. This needs
    to be specified, when the equations contain summations, for which the upper bound is given
    by a Symbolic.

"""
function create_value_map(sym::Sym{Parameter, IndexedVariable}, values::Vector;mapping::Dict{SymbolicUtils.Sym,Int64}=Dict{SymbolicUtils.Sym,Int64}())
    iVar = sym.metadata
    if iVar.ind.rangeN isa SymbolicUtils.Sym
        if iVar.ind.rangeN in keys(mapping)
            range1 = mapping[iVar.ind.rangeN]
        else
            error("Can not evaluate without a mapping")
        end
    else
        range1 = iVar.ind.rangeN
    end
    if range1 != length(values)
        error("different length of index-range and given values!")
    end
    dict = Dict{Sym{Parameter, Base.ImmutableDict{DataType, Any}},ComplexF64}()
    for i = 1:range1
        push!(dict,(SingleNumberedVariable(iVar.name,i) => values[i]))
    end
    return dict
end
function create_value_map(sym::Sym{Parameter, IndexedVariable}, value::Number)
    iVar = sym.metadata
    dict = Dict{Sym{Parameter, Base.ImmutableDict{DataType, Any}},ComplexF64}()
    if iVar.ind.rangeN isa SymbolicUtils.Sym
        if iVar.ind.rangeN in keys(mapping)
            range1 = mapping[iVar.ind.rangeN]
        else
            error("Can not evaluate without a mapping")
        end
    else
        range1 = iVar.ind.rangeN
    end
    for i = 1:range1
        push!(dict,(SingleNumberedVariable(iVar.name,i) => value))
    end
    return dict
end
function create_value_map(sym::Sym{Parameter,DoubleIndexedVariable},values::Matrix;mapping::Dict{SymbolicUtils.Sym,Int64}=Dict{SymbolicUtils.Sym,Int64}())
    dict = Dict{Sym{Parameter, Base.ImmutableDict{DataType, Any}},ComplexF64}()
    var = sym.metadata
    if var.ind1.rangeN isa SymbolicUtils.Sym
        if var.ind1.rangeN in keys(mapping)
            range1 = mapping[var.ind1.rangeN]
        else
            error("Can not evaluate without a mapping")
        end
    else
        range1 = var.ind1.rangeN
    end
    if var.ind2.rangeN isa SymbolicUtils.Sym
        if var.ind2.rangeN in keys(mapping)
            range2 = mapping[var.ind2.rangeN]
        else
            error("Can not evaluate without a mapping")
        end
    else
        range2 = var.ind2.rangeN
    end
    for i = 1:range1
        for j = 1:range2
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
                break
            end
        end
    else
        return typeof(arg_[1]) == IndexedOperator 
    end
    return found
end
containsIndex(term::SymbolicUtils.Term{AvgSym,Nothing},ind::Index) = ind ∈ getIndices(term)

#simplify functions not "really" needed, they are nice to have, since equation creation of Symbolics sometimes does not simplify certain terms
#function to reduce multiplication of numbers with a sum into just a sum of multiplication
function simplifyMultiplication(term::SymbolicUtils.Mul)
    args = arguments(term)
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
    return sum(adds)
end
simplifyMultiplication(x) = x
function simplifyAdd(term::SymbolicUtils.Add)
    adds = []
    for arg in arguments(term)
        if typeof(arg) <: SymbolicUtils.Mul
            push!(adds,simplifyMultiplication(arg))
        else
            push!(adds,arg)
        end
    end
    return sum(adds)
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
function simplifyMeanfieldEquations(me::AbstractMeanfieldEquations)
    eq_after = []
    op_after = []
    for eq in me.equations
        push!(eq_after,simplifyEquation(eq))
    end
    for eq in eq_after
        push!(op_after,undo_average(eq))
    end
    if me isa MeanfieldEquations
        return MeanfieldEquations(eq_after,op_after,me.states,me.operators,me.hamiltonian,me.jumps,me.jumps_dagger,me.rates,me.iv,me.varmap,me.order)
    elseif me isa IndexedMeanfieldEquations
        return IndexedMeanfieldEquations(eq_after,op_after,me.states,me.operators,me.hamiltonian,me.jumps,me.jumps_dagger,me.rates,me.iv,me.varmap,me.order)
    end
end 

#functions for simplifying the indexed_complete function
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
function getAvrgs(sum::SymbolicUtils.Sym{Parameter,SpecialIndexedAverage})
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
getNumber(x::NumberedOperator) = [acts_on(x) + x.numb]
getNumber(x::QMul) = acts_on(x)
getNumber(x) = [acts_on(x)] # this is so that, any other operator still behaves the same as before

function _to_expression(x::NumberedOperator) 
    x.op isa Transition && return :( NumberedOperator($(x.op.name),$(x.numb),$(x.op.i),$(x.op.j)) )
    x.op isa Destroy && return :(NumberedDestroy($(x.op.name),$(x.numb)))
    x.op isa Create && return :(dagger(NumberedDestroy($(x.op.name),$(x.numb))))
end
_to_expression(x::SymbolicUtils.Sym{Parameter,IndexedAverageSum}) = :( IndexedAverageSum($(_to_expression(x.metadata.term)),$(x.metadata.sumIndex.name),$(x.metadata.sumIndex.rangeN),$(writeNEIs(x.metadata.nonEqualIndices))) )
_to_expression(x::SymbolicUtils.Sym{Parameter,SpecialIndexedAverage}) = :($(x.metadata.term))
_to_expression(x::SymbolicUtils.Sym{Parameter,IndexedAverageDoubleSum}) = :( IndexedAverageDoubleSum($(_to_expression(x.metadata.innerSum)),$(x.metadata.sumIndex.name),$(x.metadata.sumIndex.rangeN),$(writeNEIs(x.metadata.nonEqualIndices))) )

@latexrecipe function f(s_::SymbolicUtils.Sym{Parameter,IndexedAverageSum})
    s = s_.metadata
    neis = writeNEIs(s.nonEqualIndices)

    ex = latexify(s.term)
    sumString = nothing
    if neis != ""
        sumString = L"$\underset{%$(s.sumIndex.name) ≠%$(neis) }{\overset{%$(s.sumIndex.rangeN)}{\sum}}$ %$(ex)"
    else
        sumString = L"$\underset{%$(s.sumIndex.name)}{\overset{%$(s.sumIndex.rangeN)}{\sum}}$ %$(ex)"
    end
    return sumString
end
Base.isequal(x::Missing,y::SymbolicUtils.Symbolic) = false

SymbolicUtils.simplify(sym::SymbolicUtils.Sym{Parameter,SpecialIndexedAverage}) = SpecialIndexedAverage(SymbolicUtils.simplify(sym.metadata.term),sym.metadata.indexMapping)

#Numeric Conversion of NumberedOperators
to_numeric(op::NumberedOperator,b::QuantumOpticsBase.Basis; kwargs...) = to_numeric(op.op,b;kwargs...)# this does not make sense, since a numbered operator always needs to act on a composite basis (except the basis has trivial order)
function to_numeric(op::NumberedOperator,b::QuantumOpticsBase.CompositeBasis; kwargs...) 
    #keep in mind: this function does not have a check for hilbert-spaces, meaning it can produce wrong output
    aon = getNumber(op)[1] - 1
    op_ = _to_numeric(op.op,b.bases[aon];kwargs...)
    return QuantumOpticsBase.embed(b,aon,op_)
end