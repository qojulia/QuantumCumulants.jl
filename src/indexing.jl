#Main class for indexing, here indices, sums and indexed operators and variables are defined with the corresponding calculus.
#Many helping functions used in the different classes are also defined here.

#Keep in mind: I used alot of pattern-matching and object oriented methods, since this is how I like (learned) to code, 
#however this might lead to way more lines of code, than actually needed

#=
    ToDo:

    - correct display of averages of sums -> rewrite latexify
    - Implement some sort of function, that compares indexed-evaluated meanfield equations with the ones got from non-indexed procedure (maybe compare via ODESystem?)
=#

import SymbolicUtils
import SymbolicUtils: substitute

import Symbolics
import TermInterface

import SciMLBase

import ModelingToolkit
const MTK = ModelingToolkit

using Combinatorics: partitions, combinations

const NO_METADATA = SymbolicUtils.NO_METADATA

#this section can be removed once integrated in the QuantumCumulants.jl file
source_metadata(source, name) = 
    Base.ImmutableDict{DataType, Any}(Symbolics.VariableSource, (source, name))

include("hilbertspace.jl")
include("qnumber.jl")
include("cnumber.jl")
include("fock.jl")
include("nlevel.jl")
include("equations.jl")
include("meanfield.jl")
include("average.jl")
include("utils.jl")
include("diffeq.jl")
include("correlation.jl")
include("cluster.jl")
include("scale.jl")
include("latexify_recipes.jl")
include("printing.jl")

include("QuantumCumulants.jl")

const Summable = Union{<:QNumber,<:CNumber}

struct AvgSym <: CNumber end

struct Index #main tool
    hilb::HilbertSpace
    name::Symbol
    rangeN::Union{<:Symbol,<:Number}
end
struct IndexedVariable <: CNumber #just a symbel, that can be manipulated via the metadata field
    name::Symbol
    ind::Index
    function IndexedVariable(name,ind)
        metadata = new(name,ind)
        return SymbolicUtils.Sym{Parameter, IndexedVariable}(Symbol("$(name)$(ind.name)"), metadata)
    end
end
struct DoubleIndexedVariable <: CNumber #just a symbol, that can be manipulated via the metadata field
    name::Symbol
    ind1::Index
    ind2::Index
    canHaveSame::Bool
    function DoubleIndexedVariable(name,ind1,ind2,canHaveSame)
        if !(canHaveSame) && (ind1 == ind2)
            return 0
        end
        metadata = new(name,ind1,ind2,canHaveSame)
        return SymbolicUtils.Sym{Parameter, DoubleIndexedVariable}(Symbol("$(name)$(ind1.name)$(ind2.name)"), metadata)
    end
end
struct IndexedOperator <: QNumber #An operator with an index, for now only transition operators are possible to be declared like this
    op::Transition
    ind::Index
end
struct IndexedSingleSum <:QTerm #Sum with an index, the term inside the sum must be a multiplication, either a QMul or a Symbolic one
    term::Summable
    sumIndex::Index
    nonEqualIndices::Vector{Index}  #indces, not equal to the summation index
    function IndexedSingleSum(term,sumIndex,nonEqualIndices) #rather expensive constructor to make sure Sums are created as they should
        if (typeof(term) == QMul && (SymbolicUtils._iszero(term.arg_c) || term.arg_c == 0)) || term === 0
            return 0
        else
            if (typeof(term) == IndexedVariable || typeof(term) == IndexedOperator)
                if term.ind == sumIndex
                    return new(term,sumIndex,nonEqualIndices)
                else
                    return (sumIndex.rangeN - length(nonEqualIndices)) * term
                end
            end
            if typeof(term) == QAdd
                sums = []
                for arg in term.arguments
                    push!(sums, IndexedSingleSum(arg,sumIndex,nonEqualIndices))
                end
                return +(sums...)
            end
            if typeof(term) <: SymbolicUtils.Add
                sums = []
                for arg in arguments(term)
                    push!(sums, IndexedSingleSum(arg,sumIndex,nonEqualIndices))
                end
                return +(sums...)
            end
            NEI = Index[]
            NEI_ = copy(nonEqualIndices)
            for arg in term.args_nc
                if typeof(arg) == IndexedOperator || typeof(arg) == IndexedVariable
                    if arg.ind == sumIndex || arg.ind in NEI_
                        continue
                    else
                    push!(NEI,arg.ind)
                    push!(NEI_,arg.ind)
                    end
                end
            end
            if length(NEI) == 0 #NEI are newly found indices of all operators that do not have the summation index, or are not already in the non equals list
                #in this if-condition all operators always commute with the summation index (since there are no other indices left)
                args = copy(term.args_nc)
                args_ = orderByIndex(args,[sumIndex]) #here all operators in the sum comute with operators indexed with the summation index -> push them in front
                term_ = merge_commutators(term.arg_c,args_) #merge operators again, since their order in the sum has changed
                if term_ == 0 || SymbolicUtils._iszero(term_)
                    return 0
                end 
                return new(term_,sumIndex,NEI_)
            end
            addTerms = []
            for i = 1:length(NEI) #NEI are the newly found Indices of all the ops that do not have the summation index
                ind = NEI[i]
                #when adding a new index to the list of non equals, all the following insertions for the summation index have the
                #condition, that they now can no longer be equal to any of the already inserted indices
                #for example: in first iteration i -> j => i ≠ j
                #second iteration i ≠ j, i -> k => i ≠ k (in sum); j ≠ k (for the extra term)
                if length(addTerms) > 0
                    indexMapping = Tuple{Index,Index}[]
                    for j = 1:i
                        if i != j
                            push!(indexMapping,(NEI[j],NEI[i]))
                        end
                    end
                    push!(addTerms, reorder(changeIndex(term,sumIndex,ind),indexMapping))
                else
                    push!(addTerms,changeIndex(term,sumIndex,ind))
                end
            end
            args = copy(term.args_nc)
            args_ = orderByIndex(args,[sumIndex]) #here all operators in the sum comute with operators indexed with the summation index -> push them in front
            term_ = merge_commutators(term.arg_c,args_) #merge operators again, since their order in the sum has changed
            sort!(NEI_,by=getIndName)
            return +(IndexedSingleSum(term_,sumIndex,NEI_),addTerms...)
        end
    end
end
struct SpecialIndexedTerm <: QNumber    #A term, not in a sum, that has a condition on the indices, for example σⱼ*σₖ with condition j≠k
    term::QMul
    indexMapping::Vector{Tuple{Index,Index}}    #The conditions on indices are given via this tuple-vector, each tuple representing one condition (not to be confused with the numbered ones in averages)
    function SpecialIndexedTerm(term,indexMapping)
        if length(indexMapping) == 0
            return term
        elseif typeof(term) == IndexedOperator
            return term
        elseif typeof(term) == Destroy
            return term
        elseif typeof(term) <: SymbolicUtils.Add
            args = []
            for arg in arguments(term)
                push!(args,SpecialIndexedTerm(arg,indexMapping))
            end
            return +(args...)
        elseif typeof(term) == QAdd
            args = []
            for arg in arguments(term)
                push!(args,SpecialIndexedTerm(arg,indexMapping))
            end
            return +(args...)
        elseif SymbolicUtils._iszero(term) || term.arg_c === 0
            return 0
        else
            return new(term,indexMapping)
        end
    end
end

#Constructors:

#Operators
function IndexedOperator(op::QMul,ind::Index)
    arg_c = op.arg_c
    ops_nc = []
    for op_ in op.args_nc
        op_ind = IndexedOperator(op_,ind)
        push!(ops_nc,op_ind)
    end
    return QMul(arg_c,ops_nc)
end
function IndexedOperator(qadd::QAdd,ind::Index)
    terms = []
    for elem in qadd.arguments
        push!(terms,IndexedOperator(elem,ind))
    end
    return QAdd(terms)
end
function IndexedOperator(qadd::QAdd,ind::Index)
    terms = []
    for elem in qadd.arguments
        push!(terms,IndexedOperator(elem,ind))
    end
    return QAdd(terms)
end
IndexedOperator(op::SNuN,ind::Index) = op     #This is just declared, so one can ignore type-checking on numbers

#Sums
IndexedSingleSum(ops::Vector{Any},ind::Index,NEI::Vector{Index}) = IndexedSingleSum(QMul(1,ops),ind,NEI)
IndexedSingleSum(ops::QMul,ind::Index) = IndexedSingleSum(ops,ind,Index[])
IndexedSingleSum(ops::QAdd,ind::Index) = IndexedSingleSum(ops,ind,Index[])
IndexedSingleSum(op::QNumber,ind::Index) = IndexedSingleSum(op,ind,Index[])
IndexedSingleSum(ops::Number,ind::Index,NEI::Vector{Index}) = (ind.rangeN - length(NEI))*ops
IndexedSingleSum(term::QSym,ind::Index,NEI::Vector{Index}) = (ind.rangeN - length(NEI))*term
function IndexedSingleSum(term::SpecialIndexedTerm,ind::Index,NEI::Vector{Index})
    if length(term.indexMapping) == 0
        return IndexedSingleSum(term.term,ind,NEI)
    else
        NEI_ = copy(NEI)
        for tuple in term.indexMapping
            if first(tuple) == ind && last(tuple) ∉ NEI_
                push!(NEI_, last(tuple))
            elseif last(tuple) == ind && first(tuple) ∉ NEI_
                push!(NEI_, first(tuple))
            end
        end
        return IndexedSingleSum(term.term,ind,NEI_)
    end
end

#hilberts
hilbert(ind::Index) = ind.hilb
hilbert(op::IndexedOperator) = op.ind.hilb
hilbert(var::IndexedVariable) = var.ind.hilb
hilbert(indSum::IndexedSingleSum) = indSum.sumIndex.hilb
#hilbert(qmul::QMul) = hilbert(qmul.args_nc[1])
hilbert(x::SpecialIndexedTerm) = hilbert(x.term)

#Basic functions for indexed Operators
import Base: *, +, -

function +(sum1::IndexedSingleSum,sum2::IndexedSingleSum)
    if (sum1.sumIndex == sum2.sumIndex) && check_sign(sum1,sum2) && check_term(sum1,sum2) #check if summation of sums results in 0
        return 0
    end
    return QAdd([sum1,sum2])
end
function +(sum1::IndexedSingleSum,sum2::IndexedOperator)
    QAdd([sum1,sum2])
end
function +(sum1::IndexedSingleSum,sum2::IndexedVariable)
    QAdd([sum1,sum2])
end
function +(op1::IndexedOperator,op2::IndexedOperator)
    check_hilbert(op1,op2)
    return QAdd([op1,op2])
end
function +(op1::IndexedVariable,op2::IndexedOperator)
    check_hilbert(op1,op2)
    return QAdd([op1,op2])
end
#Number
function +(a::QNumber,op::IndexedOperator)
    check_hilbert(a,op)
    return QAdd([a,op])
end
function +(a::QNumber,op::IndexedVariable)
    check_hilbert(a,op)
    return QAdd([a,op])
end
function +(a::QNumber,op::IndexedSingleSum)
    return QAdd([a,op])
end
#QAdd
function +(qadd::QAdd, indO::IndexedOperator)
    args = copy(qadd.arguments)
    push!(args,indO)
    return QAdd(args)
end
function +(qadd::QAdd, indV::IndexedVariable)
    args = copy(qadd.arguments)
    push!(args,indV)
    return QAdd(args)
end
function +(qadd::QAdd, sum::IndexedSingleSum)
    args = copy(qadd.arguments)
    push!(args,sum)
    return QAdd(args)
end
function +(sum::IndexedSingleSum,qadd::QAdd)
    args = copy(qadd.arguments)
    push!(args,sum)
    return QAdd(args)
end
+(a::IndexedSingleSum, b::SNuN) = QAdd([a,b])
+(a::SNuN,b::IndexedSingleSum) = +(b,a)
#QMul
function +(qmul::QMul,sum::IndexedSingleSum)
    args = [qmul,sum]
    return QAdd(args)
end
function +(qmul::QMul,indO::IndexedOperator)
    args = [qmul,indO]
    return QAdd(args)
end
function +(qmul::QMul,indV::IndexedVariable)
    args = [qmul,indV]
    return QAdd(args)
end
# Special terms
function +(elem::SpecialIndexedTerm,x::QAdd)
    args = copy(x.arguments)
    push!(args,elem)
    return QAdd(args)
end
+(x::QAdd,elem::SpecialIndexedTerm) = +(elem,x)
function +(elem::SpecialIndexedTerm,x::SNuN)
    return QAdd([elem,x])
end
+(x::SNuN,elem::SpecialIndexedTerm) = +(elem,x)
function +(elem::SpecialIndexedTerm,x::QNumber)
    return QAdd([elem,x])
end
+(x::QNumber,elem::SpecialIndexedTerm) = +(elem,x)
+(a::SpecialIndexedTerm,b::SpecialIndexedTerm) = QAdd([a,b])
+(a::SpecialIndexedTerm,b::IndexedOperator) = QAdd([a,b])
+(a::IndexedOperator,b::SpecialIndexedTerm) = QAdd([a,b])
+(a::SpecialIndexedTerm,b::IndexedSingleSum) = QAdd([a,b])
+(a::IndexedSingleSum,b::SpecialIndexedTerm) = QAdd([a,b])

#Multiplications
#Sums
function *(sum::IndexedSingleSum,qmul::QMul)
    args_nc = qmul.args_nc
    arg_c = qmul.arg_c
    newSum = sum
    for i = 1:length(args_nc)
        newSum = newSum*args_nc[i]
    end
    return arg_c * newSum
end
function *(qmul::QMul,sum::IndexedSingleSum)
    args_nc = qmul.args_nc
    arg_c = qmul.arg_c
    newSum = sum
    len = length(args_nc)
    for i = 1:len
        newSum = args_nc[len+1-i]*newSum
    end
    return arg_c * newSum
end
function *(sum::IndexedSingleSum,elem::QNumber)
    NEIds = copy(sum.nonEqualIndices)
    if (typeof(elem) == IndexedOperator || typeof(elem) == IndexedVariable) && !(elem.ind == sum.sumIndex) && (elem.ind ∉ NEIds)
        qaddterm = nothing
        term = sum.term
        if length(NEIds) == 0
            extraterm = changeIndex(term,sum.sumIndex,elem.ind)
            qaddterm = extraterm*elem
        else
            specNEIs = Tuple{Index,Index}[]
            for ind in NEIds
                tuple = (elem.ind,ind) 
                push!(specNEIs,tuple)
            end
            extraterm_ = changeIndex(term,sum.sumIndex,elem.ind)
            extraterm = extraterm_*elem
            qaddterm = reorder(extraterm,specNEIs)
        end
        push!(NEIds,elem.ind)
        qmul = sum.term*elem
        if typeof(qmul) <: QMul
            sort!(qmul.args_nc, by=getIndName) #inside the sum everything always commutes
        end
        if (typeof(qmul) == QMul && (qmul.arg_c == 0 || qmul.args_nc == 0))
            return 0
        end
        newsum = IndexedSingleSum(qmul,sum.sumIndex,NEIds)
    
        if newsum == 0
            return qaddterm
        elseif qaddterm == 0
            return newsum
        end
    
        return QAdd([newsum,qaddterm])
    end
    qmul = sum.term*elem
    if (typeof(qmul) == QMul && (qmul.arg_c == 0 || qmul.args_nc == 0))
        return 0
    end
    return IndexedSingleSum(qmul,sum.sumIndex,NEIds)
end
function *(elem::QNumber,sum::IndexedSingleSum)
    NEIds = copy(sum.nonEqualIndices)
    if ((typeof(elem) == IndexedOperator) || (typeof(elem) == IndexedVariable)) && !(elem.ind == sum.sumIndex) && (elem.ind ∉ NEIds)
        push!(NEIds,elem.ind)
        term = sum.term
        extraterm = changeIndex(term,sum.sumIndex,elem.ind)
        qaddterm = elem*extraterm
        qmul = elem*sum.term
        indMap = Vector{Tuple{Index,Index}}()
        for ind in NEIds
            push!(indMap,(sum.sumIndex,ind))
        end
        push!(indMap,(sum.sumIndex,elem.ind))
        if typeof(qmul) <: QMul
            qmul = reorder(qmul,indMap)
        end
        if (typeof(qmul) == QMul && (qmul.arg_c == 0 || qmul.args_nc == 0))
            return 0
        end
        newsum = IndexedSingleSum(qmul,sum.sumIndex,NEIds)

        if newsum == 0
            return qaddterm
        elseif qaddterm == 0
            return newsum
        end
        
        return QAdd([newsum,qaddterm])
    end
    qmul = elem*sum.term
    if (typeof(qmul) == QMul && (qmul.arg_c == 0 || qmul.args_nc == 0))
        return 0
    end
    return IndexedSingleSum(qmul,sum.sumIndex,NEIds)
end
*(elem::SNuN,sum::IndexedSingleSum) = IndexedSingleSum(elem*sum.term,sum.sumIndex,sum.nonEqualIndices) #put elements from outside into sum
*(sum::IndexedSingleSum,elem::SNuN) = *(elem,sum)

-(sum::IndexedSingleSum,sum2::IndexedSingleSum) = sum + -1*sum2
-(sum::IndexedSingleSum,op::QNumber) = sum + -1*op
-(op::QNumber,sum::IndexedSingleSum) = -1*sum + op
-(op::Any,sum::IndexedSingleSum) = -1*sum + op
-(sum::IndexedSingleSum, op::Any) = -1*op + sum

-(op::IndexedOperator) = -1*op
-(op::IndexedVariable) = -1*op

*(a::Create,b::IndexedOperator) = QMul(1,[a,b])
*(b::IndexedOperator,a::Create) = QMul(1,[a,b])
*(a::Destroy,b::IndexedOperator) = QMul(1,[a,b])
*(b::IndexedOperator,a::Destroy) = QMul(1,[a,b])

*(a::Create,b::IndexedVariable) = QMul(1,[b,a])
*(b::IndexedVariable,a::Create) = QMul(1,[b,a])
*(a::Destroy,b::IndexedVariable) = QMul(1,[b,a])
*(b::IndexedVariable,a::Destroy) = QMul(1,[b,a])

*(a::IndexedOperator,b::IndexedVariable) = QMul(1,[b,a])
*(b::IndexedVariable,a::IndexedOperator) = QMul(1,[b,a])

function *(op1::IndexedOperator,op2::IndexedOperator)
    if op1.ind == op2.ind
        return IndexedOperator(op1.op*op2.op,op1.ind)
    else
        return QMul(1,[op1,op2])
    end
end
function *(a::IndexedOperator, b::SNuN)
    SymbolicUtils._iszero(b) && return b
    SymbolicUtils._isone(b) && return a
    return QMul(b,[a])
end
function *(a::IndexedOperator, b::QMul)
    check_hilbert(a, b)
    args_nc = vcat(a,b.args_nc)
    sort!(args_nc,by=acts_on)
    return merge_commutators(b.arg_c,args_nc)
end
function *(a::QMul, b::IndexedOperator)
    check_hilbert(a, b)
    args_nc = vcat(a.args_nc,b)
    sort!(args_nc, by=acts_on)
    return merge_commutators(a.arg_c,args_nc)
end
function *(a::IndexedVariable, b::QMul)
    check_hilbert(a, b)
    args_nc = vcat(a,b.args_nc)
    sort!(args_nc,by=acts_on)
    return merge_commutators(b.arg_c,args_nc)
end
function *(a::QMul, b::IndexedVariable)
    check_hilbert(a, b)
    args_nc = vcat(a.args_nc, b)
    sort!(args_nc,by=acts_on)
    return merge_commutators(a.arg_c,args_nc)
end
# Special terms
function *(x, term::SpecialIndexedTerm)
    return reorder(x*term.term,term.indexMapping)
end
function *(term::SpecialIndexedTerm,x)
    return reorder(term.term*x,term.indexMapping)
end

#acts on
acts_on(op::IndexedOperator) = acts_on(op.op)
acts_on(var::IndexedVariable) = 0
acts_on(var::SpecialIndexedTerm) = acts_on(var.term)

get_order(x::IndexedSingleSum) = get_order(x.term)
get_order(x::SpecialIndexedTerm) = get_order(x.term)

acts_on(indSum::IndexedSingleSum) = acts_on(indSum.term)

#extra commutators
#Indexed operators, evaluate the commutator directly, if 2 indexed ops have the same index
function commutator(op1::IndexedOperator,op2::IndexedOperator) 
    if (op1.ind == op2.ind)
        return IndexedOperator(commutator(op1.op,op2.op),op1.ind)
    else
        return op1*op2 - op2*op1
    end
end
commutator(op1::IndexedOperator,var::IndexedVariable) = 0
commutator(op1::IndexedOperator,b::SNuN) = 0
function commutator(a::IndexedOperator,b::QAdd)
    args = []
    for b_∈b.arguments
        c = commutator(a,b_)
        push_or_append_nz_args!(args, c)
    end
    isempty(args) && return 0
    return QAdd(args)
end

#adjoint
Base.adjoint(op::IndexedOperator) = IndexedOperator(Base.adjoint(op.op),op.ind)
Base.adjoint(op::IndexedVariable) = IndexedVariable(op.name,Base.adjoint.(op.values),op.ind)
Base.adjoint(op::IndexedSingleSum) = IndexedSingleSum(Base.adjoint(op.term),op.sumIndex,op.nonEqualIndices)

#Base Functionalities
#Hashing
function Base.hash(op::IndexedOperator, h::UInt)
    n = fieldcount(IndexedOperator)
    if n == 2
        # These three fields need to be defined for any QSym
        return hash(IndexedOperator, hash(op.ind, hash(op.op, h)))
    else
        # If there are more we'll need to iterate through
        h_ = copy(h)
        for k = n:-1:4
            if fieldname(typeof(op), k) !== :metadata
                h_ = hash(getfield(IndexedOperator, k), h_)
            end
        end
        return hash(IndexedOperator, hash(op.ind, hash(op.op, h_)))
    end
end
function Base.hash(var::IndexedVariable, h::UInt)
    n = fieldcount(IndexedVariable)
    if n == 2
        # These three fields need to be defined for any QSym
        return hash(IndexedVariable, hash(var.ind, hash(var.name, h)))
    else
        # If there are more we'll need to iterate through
        h_ = copy(h)
        for k = n:-1:4
            if fieldname(typeof(var), k) !== :metadata
                h_ = hash(getfield(IndexedVariable, k), h_)
            end
        end
        return hash(IndexedVariable, hash(var.ind, hash(var.name, h_)))
    end
end
function Base.hash(ind::Index, h::UInt)
    n = fieldcount(Index)
    if n == 3
        # These three fields need to be defined for any QSym
        return hash(Index, hash(ind.hilb, hash(ind.name, hash(ind.rangeN, h))))
    else
        # If there are more we'll need to iterate through
        h_ = copy(h)
        for k = n:-1:4
            if fieldname(typeof(ind), k) !== :metadata
                h_ = hash(getfield(Index, k), h_)
            end
        end
        return hash(Index, hash(ind.hilb, hash(ind.name, hash(ind.rangeN, h))))
    end
end

#Ordering of indices
Base.isless(a::IndexedOperator, b::IndexedOperator) = a.op.name < b.op.name
Base.isless(a::QMul, b::QMul) = isless(a.args_nc, b.args_nc)
Base.isless(a::IndexedOperator,b::QSym) = a.op.name < b.name
Base.isless(a::QSym,b::IndexedOperator) = a.name < b.op.name
Base.isless(a::IndexedVariable,b::IndexedVariable) = a.name < b.name
Base.isless(nothing,b::Symbol) = true
Base.isless(b::Symbol,nothing) = false

Base.isless(a::IndexedVariable,b::QSym) = true
Base.isless(b::QSym,a::IndexedVariable) = false 

Base.isless(a::Index,b::Index) = a.name < b.name
Base.isless(a::IndexedSingleSum,b::IndexedSingleSum) = Base.isless(a.sumIndex,b.sumIndex)

Base.isequal(ind1::Index,ind2::Index) = (ind1.name == ind2.name) && (ind1.rangeN == ind2.rangeN) && (ind1.hilb == ind2.hilb)
Base.:(==)(ind1::Index,ind2::Index) = (ind1.name == ind2.name) && (ind1.rangeN == ind2.rangeN) && (ind1.hilb == ind2.hilb)
Base.isequal(a::SpecialIndexedTerm,b::SpecialIndexedTerm) = isequal(a.term, b.term) && isequal(a.indexMapping, b.indexMapping)

#checks if two sums have opposite numeric values
check_sign(a::IndexedSingleSum,b::IndexedSingleSum) = isequal(a.term.arg_c, -1*b.term.arg_c)
function check_term(a::IndexedSingleSum,b::IndexedSingleSum)
    isequal(a.sumIndex,b.sumIndex) || return false
    isequal(a.term.arg_c, b.term.arg_c) || isequal(a.term.arg_c, -1*b.term.arg_c) || return false
    length(a.term.args_nc)==length(b.term.args_nc) || return false
    length(a.nonEqualIndices) == length(b.nonEqualIndices) || return false
    sort!(a.nonEqualIndices, by=getIndName) == sort!(b.nonEqualIndices, by=getIndName) || return false
    for (arg_a, arg_b) ∈ zip(sort!(a.term.args_nc,by=getIndName), sort!(b.term.args_nc,by=getIndName))
        isequal(arg_a,arg_b) || return false
    end
    return true
end
function Base.isequal(a::IndexedSingleSum, b::IndexedSingleSum)
    isequal(a.sumIndex,b.sumIndex) || return false
    isequal(a.term.arg_c, b.term.arg_c) || return false
    length(a.term.args_nc)==length(b.term.args_nc) || return false
    for (arg_a, arg_b) ∈ zip(orderByIndex(a.term.args_nc,[a.sumIndex]), orderByIndex(b.term.args_nc,[b.sumIndex]))
        isequal(arg_a,arg_b) || return false
    end
    return true
end
hilbert(a::SNuN) = 0

#Function that changes the index of a sum term into a different indexed term
#used for evaluating the extra terms when multiplying a sum with an operator with different index
#return a new QMul with indices swapped: from -> to index
function changeIndex(term::QMul, from::Index, to::Index)
    arg_c = term.arg_c
    arg_c_ = term.arg_c
    args_nc = copy(term.args_nc)
    for i = 1:length(args_nc)
        elem = args_nc[i]
        if typeof(elem) == IndexedOperator && elem.ind == from
            elem = IndexedOperator(elem.op,to)
        end
        args_nc[i] = elem #inplace exchange of element
    end
    if typeof(arg_c_) == IndexedVariable && arg_c_.ind == from
        arg_c = IndexedVariable(arg_c_.name,to)
    elseif typeof(arg_c_) <: SymbolicUtils.Mul
        args = copy(arguments(arg_c_))
        for i = 1:length(args)
            if typeof(args[i]) == SymbolicUtils.Sym{Parameter,IndexedVariable}
                if args[i].metadata.ind == from
                    var = args[i].metadata      #the actual indexed-value
                    args[i] = IndexedVariable(var.name,to) # inplace exchange
                end
            elseif typeof(args[i]) == SymbolicUtils.Sym{Parameter, DoubleIndexedVariable}
                if args[i].metadata.ind1 == args[i].metadata.ind2 && args[i].metadata.ind1 == from
                    var = args[i].metadata
                    args[i] = DoubleIndexedVariable(var.name,to,to,var.canHaveSame)
                elseif args[i].metadata.ind1 == from
                    var = args[i].metadata
                    args[i] = DoubleIndexedVariable(var.name,to,var.ind2,var.canHaveSame)
                elseif args[i].metadata.ind2 == from
                    var = args[i].metadata
                    args[i] = DoubleIndexedVariable(var.name,var.ind1,to,var.canHaveSame)
                end
            end
        end
        arg_c = *(args...)
    elseif typeof(arg_c_) == SymbolicUtils.Sym{Parameter, IndexedVariable} && arg_c_.metadata.ind == from
        arg_c = IndexedVariable(arg_c_.metadata.name,to)
    elseif  typeof(arg_c_) == SymbolicUtils.Sym{Parameter, DoubleIndexedVariable}
        DIndV = arg_c_.metadata
        if DIndV.ind1 == DIndV.ind2 && DIndV.ind1 == from
            arg_c = DoubleIndexedVariable(DIndV.name,to,to,DIndV.canHaveSame)
        elseif DIndV.ind1 == from
            arg_c = DoubleIndexedVariable(DIndV.name,to,DIndV.ind2,DIndV.canHaveSame)
        elseif DIndV.ind2 == from
            arg_c = DoubleIndexedVariable(DIndV.name,DIndV.ind1,to,DIndV.canHaveSame)
        end
    end
    if isempty(args_nc) || isequal(arg_c,0) || SymbolicUtils._iszero(args_nc)
        return 0
    end
    return merge_commutators(arg_c,args_nc)
end
function changeIndex(term::Term{AvgSym, Nothing}, from::Index,to::Index)
    qmul = arguments(term)[1]
    return average(changeIndex(qmul,from,to))
end
function changeIndex(op::IndexedOperator,from::Index,to::Index)
    if op.ind == from
        return IndexedOperator(op.op,to)
    else
        return op
    end
end
function changeIndex(ops::Vector,from::Index,to::Index)
    ops_ = copy(ops)
    for i = 1:length(ops_)
        ops_[i] = changeIndex(ops_[i],from,to)
    end
    return ops_
end
changeIndex(op::SymbolicUtils.Sym{Parameter,IndexedVariable},from::Index,to::Index) = op.metadata.ind == from ? IndexedVariable(op.metadata.name,to) : op
function changeIndex(op::SymbolicUtils.Sym{Parameter,DoubleIndexedVariable},from::Index,to::Index)
    if op.metadata.ind1 == from
        if op.metadata.ind1 == op.metadata.ind2 && op.metadata.canHaveSame
            return DoubleIndexedVariable(op.metadata.name,to,to,op.metadata.canHaveSame)
        elseif op.metadata.ind1 == op.metadata.ind2
            return 0
        else
            return DoubleIndexedVariable(op.metadata.name,to,op.metadata.ind2,op.metadata.canHaveSame)
        end
    elseif op.metadata.ind2 == from
        return DoubleIndexedVariable(op.metadata.name,op.metadata.ind1,to,op.metadata.canHaveSame)
    end
end
changeIndex(x,from::Index,to::Index) = x

ismergeable(a::IndexedOperator,b::IndexedOperator) = (a.ind == b.ind) ? ismergeable(a.op,b.op) : false

getIndName(op::IndexedOperator) = op.ind.name
getIndName(ind::Index) = ind.name
getIndName(x) = Symbol()

SymbolicUtils.istree(a::IndexedSingleSum) = false
SymbolicUtils.arguments(a::IndexedSingleSum) = SymbolicUtils.arguments(a.term)

get_order(::IndexedOperator) = 1
get_order(::IndexedVariable) = 0

#It is assumed that the term for which this operation is done already commutes with indices inside the indices-Vector
function orderByIndex(vec::Vector,indices::Vector{Index})
    vec_ = copy(vec)
    frontfront = []
    front = []
    back = []
    frontfront = filter(x -> !(typeof(x) == IndexedOperator),vec_)
    front = filter(x -> typeof(x) == IndexedOperator && x.ind in indices,vec_)
    back = filter(x -> typeof(x) == IndexedOperator && x.ind ∉ indices,vec_)
    sort!(front,by=getIndName)
    return vcat(frontfront,front,back)
end
function orderByIndex(qmul::QMul,inds::Vector{Index})
    return merge_commutators(qmul.arg_c,orderByIndex(qmul.args_nc,inds))
end
orderByIndex(qadd::QAdd) = +(orderByIndex.(qadd.arguments)...)
orderByIndex(x) = x
#Reorder function: given a tuple vector of indices meaning for each tuple: first ≠ second
#-> go through the term given and exchange 2 ops when the second has "lower" (i.e. its name is first in the alphabet) index than the first one
#-> results in a term, ordered by its commutating indices
function reorder(param::QMul,indexMapping::Vector{Tuple{Index,Index}})
    term = copy(param.args_nc)
    carg = param.arg_c
    indOps = []
    others = []
    for i = 1:length(term) #Split into indexed ops and non indexed ops
        if typeof(term[i]) == IndexedOperator
            push!(indOps,term[i])
        else
            push!(others,term[i])
        end
    end 
    while true #go over all ops ind indexed ops -> order by 
        finish = true
        for i = 1:(length(indOps)-1)
            if ((indOps[i].ind,indOps[i+1].ind) in indexMapping || (indOps[i+1].ind,indOps[i].ind) in indexMapping) && (indOps[i+1].ind < indOps[i].ind)
                temp = indOps[i+1]
                indOps[i+1] = indOps[i]
                indOps[i] = temp
                finish = false
            end
        end
        if finish
            break
        end
    end
    qmul = merge_commutators(carg,vcat(others,indOps))
    return SpecialIndexedTerm(qmul,indexMapping)
end
reorder(sum::IndexedSingleSum,indexMapping::Vector{Tuple{Index,Index}}) = IndexedSingleSum(reorder(sum.term,indexMapping),sum.sumIndex,sum.nonEqualIndices)
reorder(op::SpecialIndexedTerm) = reorder(op.term,op.indexMapping)
function reorder(term::SymbolicUtils.Add,indexMapping::Vector{Tuple{Index,Index}})
    args = []
    for arg in arguments(term)
        push!(args,reorder(arg,indexMapping))
    end
    if length(args) == 0
        return 0
    end
    if length(args) == 1
        return args[1]
    end
    return +(args...)
end
function reorder(term::SymbolicUtils.Mul,indexMapping::Vector{Tuple{Index,Index}})
    args = []
    for arg in arguments(term)
        push!(args,reorder(arg,indexMapping))
    end
    if length(args) == 0
        return 0
    end
    if length(args) == 1
        return args[1]
    end
    return *(args...)
end
function reorder(term::QAdd,indexMapping::Vector{Tuple{Index,Index}})
    args = []
    for arg in arguments(term)
        push!(args,reorder(arg,indexMapping))
    end
    if length(args) == 0
        return 0
    end
    if length(args) == 1
        return args[1]
    end
    return +(args...)
end
reorder(x,indMap) = x

#Show functions
function Base.show(io::IO,op::IndexedOperator)
    op_ = op.op
    if typeof(op_) <:Transition
        write(io,Symbol(op_.name,op_.i,op_.j))
    else
        write(io,op_.name)
    end
    write(io,op.ind.name)
end
function Base.show(io::IO,elem::IndexedVariable)
    write(io,elem.name,elem.ind.name)
end
function Base.show(io::IO,indSum::IndexedSingleSum) 
    write(io, "Σ", "($(indSum.sumIndex.name)", "=1:$(indSum.sumIndex.rangeN))",)
    if !(isempty(indSum.nonEqualIndices))
        write(io,"($(indSum.sumIndex.name)≠")
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
function Base.show(io::IO,op::SpecialIndexedTerm)
    if !isempty(op.indexMapping)
        Base.write(io,"(")
    end
    for i = 1:length(op.indexMapping)
        Base.write(io,first(op.indexMapping[i]).name)
        Base.write(io,"≠")
        Base.write(io,last(op.indexMapping[i]).name)
        if i != length(op.indexMapping)
            Base.write(io,";")
        else
            Base.write(io,")")
        end
    end
    Base.show(io,op.term)
end