#Main class for indexing, here indices, sums and indexed operators and variables are defined with the corresponding calculus.
#Many helping functions used in the different classes are also defined here.

const Ranges = Union{<:SymbolicUtils.Sym,<:Number,<:SymbolicUtils.Mul,<:SymbolicUtils.Div} #possible Types for the range of an index
const IndexableOps = Union{Transition,Create,Destroy} #every operator that can have an index

"""

    Index(hilb::HilbertSpace,name::Symbol,range::Union{Int64,Sym},aon::Int)

Defines an index, using a Symbol as a name, and a [`HilbertSpace`](@ref) for computation and
commutator-relations. Indices with all same fields will be considered equal.
See also: [`IndexedOperator`](@ref) and [`IndexedVariable`](@ref)

Fields:
======

* hilb: The whole [`HilbertSpace`](@ref), the index will be defined on.
* name: A Symbol, which defines the name of the index, and how product-terms of [`IndexedOperator`](@ref) are ordered (alphabetical)
* range: The upper bound limit of the index. This can be a SymbolicUitls.Symbolic or any Number.
* aon: Number specifying the specific [`HilbertSpace`](@ref), where the Index acts on.

"""
struct Index #main tool
    hilb::HilbertSpace
    name::Symbol
    range::Ranges
    aon::Int
end
function Index(hilb::HilbertSpace,name::Symbol,range::Ranges,specHilb::HilbertSpace)
    if hilb isa ProductSpace
        aon = findfirst(x -> isequal(x,specHilb), hilb.spaces)
        return Index(hilb,name,range,aon)
    else
        return Index(hilb,name,range,1)
    end
end

const IndexInt = Union{<:Index,<:Int64}
"""

    IndexedVariable <: CNumber
    IndexedVariable(name::Symbol,ind::Index)
    IndexedVariable(name::Symbol,ind1::Index,ind2:Index)

A indexed symbolic variable. The variable can (once equations are calculated) be easily exchanged for numerical values. Calling
a IndexedVariable using two different [`Index`](@ref) objects one can create [`DoubleIndexedVariable`](@ref) objects.
See also: [`value_map`](@ref)
"""
struct IndexedVariable <: CNumber #just a symbol, that can be manipulated via the metadata field
    name::Symbol
    ind::Index
    function IndexedVariable(name::Symbol,ind::Index)
        metadata = new(name,ind)
        return SymbolicUtils.Sym{Parameter, IndexedVariable}(Symbol("$(name)$(ind.name)"), metadata)
    end
end
"""

    DoubleIndexedVariable <: CNumber
    DoubleIndexedVariable(name::Symbol,ind1::Index,ind2::Index;identical::Bool)

A double-indexed symbolic variable. The variable can (once equations are calculated) be easily exchanged for numerical values.
See also: [`value_map`](@ref)

Fields:
======

* name: A Symbol, defining the name of the variable
* ind1: The first Index of the variable
* ind2: The second Index of the variable
* identical: A Bool, defining if the variable can have non-zero main-diagonal terms, e.g: Γᵢᵢ ≠ 0 would be specified with true.

"""
struct DoubleIndexedVariable <: CNumber #just a symbol, that can be manipulated via the metadata field
    name::Symbol
    ind1::Index
    ind2::Index
    identical::Bool
    function DoubleIndexedVariable(name,ind1,ind2;identical::Bool=true)
        if !(identical) && (ind1 == ind2)
            return 0
        end
        metadata = new(name,ind1,ind2,identical)
        return SymbolicUtils.Sym{Parameter, DoubleIndexedVariable}(Symbol("$(name)$(ind1.name)$(ind2.name)"), metadata)
    end
end
"""

    IndexedOperator <: QSym
    IndexedOperator(op::Union{Transition,Create,Destroy},ind::Index)

Operator, associated with an index.

Fields:
======

* op: Operator, either a [`Transition`](@ref), a [`Destroy`](@ref) or a [`Create`](@ref) can be defined.
* ind: The index, the operator will be associated with.

"""
struct IndexedOperator <: QSym
    op::IndexableOps
    ind::Index
    function IndexedOperator(op::IndexableOps,ind::Index)
        @assert isequal(ind.hilb,hilbert(op))
        isa(ind.hilb, ProductSpace) && (@assert isequal(acts_on(op),ind.aon))
        return new(op,ind)
    end
end

const Summable = Union{<:QNumber,<:CNumber,<:SymbolicUtils.Sym{Parameter,IndexedVariable},<:SymbolicUtils.Sym{Parameter,DoubleIndexedVariable}}

"""

    SingleSum <: QTerm

Defines a symbolic summation over a term, using one [`Index`](@ref) entity.

Fields:
======

* term: A multiplication of q-number terms. When the multiplication contains any [`IndexedOperator`](@ref) with the same index as the summation-index, a symbolic sum will be created.
* sum_index: The index, for which the summation will go over.
* non_equal_indices: (optional) A vector of indices, for which the summation-index can not be equal with.

"""
struct SingleSum{M} <: QTerm #Sum with an index, the term inside the sum must be a multiplication, either a QMul or a Symbolic one
    term::Summable
    sum_index::Index
    non_equal_indices::Vector{IndexInt}  #indices not equal to the summation index
    metadata::M
    function SingleSum(term::Summable,sum_index::Index,non_equal_indices::Vector,metadata)
        SymbolicUtils._iszero(term) ? 0 : new{typeof(metadata)}(term,sum_index,non_equal_indices,metadata)
    end
end
"""

    SpecialIndexedTerm <: QTerm

A multiplication of [`IndexedOperator`](@ref) entities, with special constraint on the index-values. For example σᵢ²² * σⱼ²² with the constraint i ≠ j

Fields:
======

* term: A multiplication of q-number terms.
* indexMapping: A Vector of [`Index`](@ref) tuples, specifying the contraints for the term. Each Tuple is considered to one constraint. e.g: (i,j) -> i ≠ j

"""
struct SpecialIndexedTerm <: QTerm    #A term, not in a sum, that has a condition on the indices, for example σⱼ*σₖ with condition j≠k
    term::Summable
    indexMapping::Vector{Tuple{Index,Index}}    #The conditions on indices are given via this tuple-vector, each tuple representing one condition (not to be confused with the numbered ones in averages)
    function SpecialIndexedTerm(term,indexMapping)
        if length(indexMapping) == 0
            return term
        elseif SymbolicUtils._iszero(term)
            return 0
        else
            return new(term,indexMapping)
        end
    end
end
const IndexedObSym = Union{IndexedOperator,SymbolicUtils.Sym{Parameter,IndexedVariable},SymbolicUtils.Sym{Parameter,DoubleIndexedVariable}}
const IndexedAdd = Union{QAdd, SymbolicUtils.Add}

SpecialIndexedTerm(term::IndexedObSym,indexMapping) = term
function SpecialIndexedTerm(add::IndexedAdd,indexMapping) 
    args = copy(arguments(add))
    return sum(reorder.(args,indexMapping))
end

function SingleSum(term::IndexedObSym, sum_index, non_equal_indices;metadata=NO_METADATA)
    term_indices = get_indices(term)
    if sum_index in term_indices
        return SingleSum(term,sum_index,non_equal_indices,metadata)
    else
        return (sum_index.range - length(non_equal_indices)) * term
    end
end

function SingleSum(term::IndexedAdd, sum_index, non_equal_indices;metadata=NO_METADATA)
    sum(SingleSum(arg,sum_index,non_equal_indices;metadata=metadata) for arg in arguments(term))
end
function SingleSum(term::QMul, sum_index, non_equal_indices;metadata=NO_METADATA)
    SymbolicUtils._iszero(term) && return 0
    NEI = Index[]
    NEI_ = copy(non_equal_indices)
    for arg in term.args_nc
        if typeof(arg) == IndexedOperator || typeof(arg) == IndexedVariable
            if arg.ind == sum_index || arg.ind in NEI_ || sum_index.aon != arg.ind.aon
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
        args_ = order_by_index(args,[sum_index]) #here all operators in the sum comute with operators indexed with the summation index -> push them in front
        term_ = 0
        if length(args_) == 1
            term_ = *(term.arg_c,args_[1])
        else
            term_ = *(term.arg_c,args_...) #merge operators again, since their order in the sum has changed
        end
        if term_ == 0 || SymbolicUtils._iszero(term_)
            return 0
        end
        return SingleSum(term_,sum_index,NEI_,metadata)
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
            push!(addTerms, reorder(change_index(term,sum_index,ind),indexMapping))
        else
            push!(addTerms,change_index(term,sum_index,ind))
        end
    end
    args = copy(term.args_nc)
    args_ = order_by_index(args,[sum_index]) #here all operators in the sum comute with operators indexed with the summation index -> push them in front
    if length(args_) == 1
        term_ = *(term.arg_c,args_[1])
    else
        term_ = *(term.arg_c,args_...) #merge operators again, since their order in the sum has changed
    end
    sort!(NEI_,by=getIndName)
    return +(SingleSum(term_,sum_index,NEI_;metadata=metadata),addTerms...)
end
function SingleSum(term::SpecialIndexedTerm,ind::Index,NEI;metadata=NO_METADATA)
    if length(term.indexMapping) == 0
        return SingleSum(term.term,ind,NEI)
    else
        NEI_ = copy(NEI)
        for tuple in term.indexMapping
            if first(tuple) == ind && last(tuple) ∉ NEI_
                push!(NEI_, last(tuple))
            elseif last(tuple) == ind && first(tuple) ∉ NEI_
                push!(NEI_, first(tuple))
            end
        end
        return SingleSum(term.term,ind,NEI_;metadata=metadata)
    end
end
SingleSum(ops::Vector{Any},ind::Index,NEI::Vector;metadata=NO_METADATA) = SingleSum(*(1,ops...),ind,NEI;metadata=metadata)
SingleSum(ops::QMul,ind::Index;metadata=NO_METADATA) = SingleSum(ops,ind,Index[];metadata=metadata)
SingleSum(ops::QAdd,ind::Index;metadata=NO_METADATA) = SingleSum(ops,ind,Index[];metadata=metadata)
SingleSum(op::QNumber,ind::Index;metadata=NO_METADATA) = SingleSum(op,ind,Index[];metadata=metadata)
SingleSum(ops::Number,ind::Index,NEI::Vector;metadata=NO_METADATA) = (ind.range - length(NEI))*ops
SingleSum(ops::SymbolicUtils.Mul,ind::Index,NEI::Vector;metadata=NO_METADATA) = (ind.range - length(NEI))*ops

SingleSum(term, sum_index, non_equal_indices;metadata=NO_METADATA) = (sum_index.range - length(non_equal_indices)) * term

function IndexedOperator(op::QMul,ind::Index)
    arg_c = op.arg_c
    ops_nc = []
    for op_ in op.args_nc
        op_ind = IndexedOperator(op_,ind)
        push!(ops_nc,op_ind)
    end
    return *(arg_c,ops_nc...)
end
function IndexedOperator(qadd::QAdd,ind::Index)
    terms = []
    for elem in qadd.arguments
        push!(terms,IndexedOperator(elem,ind))
    end
    return QAdd(terms)
end
IndexedOperator(op::SNuN,ind::Index) = op     #This is just declared, so one can ignore type-checking on numbers

#hilberts
hilbert(ind::Index) = ind.hilb
hilbert(op::IndexedOperator) = op.ind.hilb
hilbert(var::IndexedVariable) = var.ind.hilb
hilbert(indSum::SingleSum) = indSum.sum_index.hilb
hilbert(x::SpecialIndexedTerm) = hilbert(x.term)

#Basic functions for indexed Operators
import Base: *, +, -

#Multiplications
#Sums
*(sum::SingleSum,qmul::QMul) = qmul.arg_c*(*(sum,qmul.args_nc...))
function *(qmul::QMul,sum::SingleSum)
    sum_ = sum
    for i = length(qmul.args_nc):-1:1
        sum_ = qmul.args_nc[i] * sum_
    end
    return qmul.arg_c*sum_ 
end

function *(sum::SingleSum,elem::IndexedObSym)
    NEIds = copy(sum.non_equal_indices)
    if !(elem.ind == sum.sum_index) && (elem.ind ∉ NEIds) && (sum.sum_index.aon == elem.ind.aon)
        qaddterm = nothing
        term = sum.term
        if length(NEIds) == 0
            extraterm = change_index(term,sum.sum_index,elem.ind)
            qaddterm = extraterm*elem
        else
            specNEIs = Tuple{Index,Index}[]
            for ind in NEIds
                tuple = (elem.ind,ind)
                push!(specNEIs,tuple)
            end
            extraterm_ = change_index(term,sum.sum_index,elem.ind)
            qaddterm = reorder(extraterm_*elem,specNEIs)
        end
        push!(NEIds,elem.ind)
        qmul = sum.term*elem
        if qmul isa QMul
            qmul = order_by_index(qmul,[sum.sum_index])
        end
        if (qmul isa QMul && (isequal(qmul.arg_c,0) || SymbolicUtils._iszero(qmul.args_nc)))
            return 0
        end
        sort!(NEIds,by=getIndName)
        newsum = SingleSum(qmul,sum.sum_index,NEIds)

        if SymbolicUtils._iszero(newsum)
            return qaddterm
        elseif SymbolicUtils._iszero(qaddterm)
            return newsum
        end

        return QAdd([newsum,qaddterm])
    else
        qmul = sum.term*elem
        if qmul isa QMul 
            qmul = order_by_index(qmul,[sum.sum_index])
        end
        qmul isa QMul && (isequal(qmul.arg_c,0) || SymbolicUtils._iszero(qmul.args_nc)) && return 0
        return SingleSum(qmul,sum.sum_index,NEIds)
    end
end
function *(elem::IndexedObSym,sum::SingleSum)
    NEIds = copy(sum.non_equal_indices)
    if !(elem.ind == sum.sum_index) && (elem.ind ∉ NEIds) && (sum.sum_index.aon == elem.ind.aon)
        qaddterm = nothing
        term = sum.term
        if length(NEIds) == 0
            extraterm = change_index(term,sum.sum_index,elem.ind)
            qaddterm = elem*extraterm
        else
            specNEIs = Tuple{Index,Index}[]
            for ind in NEIds
                tuple = (elem.ind,ind)
                push!(specNEIs,tuple)
            end
            extraterm_ = change_index(term,sum.sum_index,elem.ind)
            qaddterm = reorder(elem*extraterm_,specNEIs)
        end
        push!(NEIds,elem.ind)
        qmul = elem*sum.term
        if qmul isa QMul
            qmul = order_by_index(qmul,[sum.sum_index])
        end
        if (qmul isa QMul && (isequal(qmul.arg_c,0) || SymbolicUtils._iszero(qmul.args_nc)))
            return 0
        end
        sort!(NEIds,by=getIndName)
        newsum = SingleSum(qmul,sum.sum_index,NEIds)

        if SymbolicUtils._iszero(newsum)
            return qaddterm
        elseif SymbolicUtils._iszero(qaddterm)
            return newsum
        end

        return QAdd([newsum,qaddterm])
    else
        qmul = elem*sum.term
        if qmul isa QMul 
            qmul = order_by_index(qmul,[sum.sum_index])
        end
        qmul isa QMul && (isequal(qmul.arg_c,0) || SymbolicUtils._iszero(qmul.args_nc)) && return 0
        return SingleSum(qmul,sum.sum_index,NEIds)
    end
end
function *(sum::SingleSum,elem::QSym)
    NEIds = copy(sum.non_equal_indices)
    qmul = sum.term*elem
    if qmul isa QMul 
        qmul = order_by_index(qmul,[sum.sum_index])
    end
    qmul isa QMul && (isequal(qmul.arg_c,0) || SymbolicUtils._iszero(qmul.args_nc)) && return 0
    return SingleSum(qmul,sum.sum_index,NEIds)
end
function *(elem::QSym,sum::SingleSum)
    NEIds = copy(sum.non_equal_indices)
    qmul = elem*sum.term
    if qmul isa QMul 
        qmul = order_by_index(qmul,[sum.sum_index])
    end
    qmul isa QMul && (isequal(qmul.arg_c,0) || SymbolicUtils._iszero(qmul.args_nc)) && return 0
    return SingleSum(qmul,sum.sum_index,NEIds)
end


*(elem::SNuN,sum::SingleSum) = SingleSum(elem*sum.term,sum.sum_index,sum.non_equal_indices) #put elements from outside into sum
*(sum::SingleSum,elem::SNuN) = *(elem,sum)

-(sum::SingleSum,sum2::SingleSum) = sum + -1*sum2
-(sum::SingleSum,op::QNumber) = sum + -1*op
-(op::QNumber,sum::SingleSum) = -1*sum + op
-(op::Any,sum::SingleSum) = -1*sum + op
-(sum::SingleSum, op::Any) = -1*op + sum

function *(op1::IndexedOperator,op2::IndexedOperator)
    if isequal(op1.ind, op2.ind)
        if op1.op isa Transition
            return IndexedOperator(op1.op*op2.op,op1.ind)
        end
        if op1.op isa Destroy && op2.op isa Create
            return op2*op1 + 1
        else
            return QMul(1,[op1,op2])
        end
    else
        return QMul(1,[op1,op2])
    end
end

# Special terms
*(a::SNuN,b::SpecialIndexedTerm) = reorder(a*b.term,b.indexMapping)
*(b::SpecialIndexedTerm,a::SNuN) = a*b
*(a::SpecialIndexedTerm,b) = reorder(a.term*b,a.indexMapping)
*(a,b::SpecialIndexedTerm) = reorder(a*b.term,b.indexMapping)

function *(a::QAdd,b::SpecialIndexedTerm)
    check_hilbert(a, b)
    args = Any[a_ * b for a_ ∈ a.arguments]
    flatten_adds!(args)
    isempty(args) && return 0
    q = QAdd(args)
    return q
end
function *(a::SpecialIndexedTerm, b::QAdd)
    check_hilbert(a, b)
    args = Any[a * b_ for b_ ∈ b.arguments]
    flatten_adds!(args)
    isempty(args) && return 0
    q = QAdd(args)
    return q
end

function *(x::IndexedObSym, term::SpecialIndexedTerm)
    map = term.indexMapping
    if x.ind ∉ first.(map) && x.ind ∉ last.(map)
        indices = get_indices(term.term)
        for ind in indices
            push!(map,(ind,x.ind))
        end
    end
    return reorder(x*term.term,map)
end
function *(term::SpecialIndexedTerm,x::IndexedObSym)
    map = term.indexMapping
    if x.ind ∉ first.(map) && x.ind ∉ last.(map)
        indices = get_indices(term.term)
        for ind in indices
            push!(map,(ind,x.ind))
        end
    end
    return reorder(term.term*x,map)
end

#acts on
acts_on(op::IndexedOperator) = acts_on(op.op)
acts_on(var::SpecialIndexedTerm) = acts_on(var.term)

get_order(x::SingleSum) = get_order(x.term)
get_order(x::SpecialIndexedTerm) = get_order(x.term)

acts_on(indSum::SingleSum) = acts_on(indSum.term)

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
    length(args) == 1 && return args[1]
    return +(args...)
end

#adjoint
Base.adjoint(op::IndexedOperator) = IndexedOperator(Base.adjoint(op.op),op.ind)
Base.adjoint(op::SingleSum) = SingleSum(Base.adjoint(op.term),op.sum_index,op.non_equal_indices)

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
function Base.hash(ind::Index, h::UInt)
    n = fieldcount(Index)
    if n == 3
        # These three fields need to be defined for any QSym
        return hash(Index, hash(ind.hilb, hash(ind.name, hash(ind.range, h))))
    else
        # If there are more we'll need to iterate through
        h_ = copy(h)
        for k = n:-1:4
            if fieldname(typeof(ind), k) !== :metadata
                h_ = hash(getfield(Index, k), h_)
            end
        end
        return hash(Index, hash(ind.hilb, hash(ind.name, hash(ind.range, h))))
    end
end

#Ordering of indices
Base.isless(a::IndexedOperator, b::IndexedOperator) = a.op.name < b.op.name
Base.isless(a::QMul, b::QMul) = isless(a.args_nc, b.args_nc)
Base.isless(a::IndexedOperator,b::QSym) = a.op.name < b.name
Base.isless(a::QSym,b::IndexedOperator) = a.name < b.op.name
Base.isless(nothing,b::Symbol) = true
Base.isless(b::Symbol,nothing) = false

Base.isless(a::Index,b::Index) = a.name < b.name
Base.isless(a::SingleSum,b::SingleSum) = Base.isless(a.sum_index,b.sum_index)

Base.isequal(ind1::Index,ind2::Index) = (ind1.name == ind2.name) && isequal(ind1.range,ind2.range) && (ind1.hilb == ind2.hilb) && isequal(ind1.aon,ind2.aon)

Base.isequal(op1::IndexedOperator,op2::IndexedOperator) = isequal(op1.op,op2.op) && isequal(op1.ind,op2.ind)
Base.:(==)(ind1::Index,ind2::Index) = isequal(ind1,ind2)
function Base.isequal(a::SpecialIndexedTerm,b::SpecialIndexedTerm)
    isequal(a.term, b.term) || return false
    isequal(length(a.indexMapping),length(b.indexMapping)) || return false
    for tuple in a.indexMapping
        if tuple ∉ b.indexMapping && (last(tuple),first(tuple)) ∉ b.indexMapping
            return false
        end
    end
    return true
end

function Base.isequal(a::SingleSum, b::SingleSum)
    if isequal(a.sum_index,b.sum_index)
        typeof(a.term) == typeof(b.term) || return false
        isequal(a.term,b.term) || return false    
    elseif isequal(a.sum_index.aon,b.sum_index.aon)
        typeof(a.term) == typeof(b.term) || return false
        isequal(a.term,change_index(b.term,b.sum_index,a.sum_index)) || return false 
    end
    return isequal(a.non_equal_indices,b.non_equal_indices)
end

#Function that changes the index of a sum term into a different indexed term
#used for evaluating the extra terms when multiplying a sum with an operator with different index
#return a new QMul with indices swapped: from -> to index
"""
    change_index(term,from::Index,to::Index)

Exchanges all occuring indices inside the given term, that are equal to the `from` to the `to` index.

Examples
========

    change_index(σⱼ²¹,j,i) = σᵢ²¹

    change_index(σⱼ²¹ * σᵢ¹²,j,i) = σᵢ²²


"""
function change_index(term::QMul, from::Index, to::Index)
    arg_c = change_index(term.arg_c,from,to)
    args_nc = [change_index(arg,from,to) for arg in copy(term.args_nc)]
    return arg_c*prod(args_nc)

end
change_index(term::Average, from::Index,to::Index) = average(change_index(arguments(term)[1],from,to))
change_index(op::IndexedOperator,from::Index,to::Index) = isequal(op.ind,from) ? IndexedOperator(op.op,to) : op

change_index(op::SymbolicUtils.Sym{Parameter,IndexedVariable},from::Index,to::Index) = isequal(op.metadata.ind,from) ? IndexedVariable(op.metadata.name,to) : op
function change_index(op::SymbolicUtils.Sym{Parameter,DoubleIndexedVariable},from::Index,to::Index)
    if op.metadata.ind1 == from
        if op.metadata.ind1 == op.metadata.ind2 && op.metadata.identical
            return DoubleIndexedVariable(op.metadata.name,to,to;identical=op.metadata.identical)
        elseif op.metadata.ind1 == op.metadata.ind2
            return 0
        else
            return DoubleIndexedVariable(op.metadata.name,to,op.metadata.ind2;identical=op.metadata.identical)
        end
    elseif op.metadata.ind2 == from
        return DoubleIndexedVariable(op.metadata.name,op.metadata.ind1,to;identical=op.metadata.identical)
    end
end
change_index(mul::SymbolicUtils.Mul,from::Index,to::Index) = prod([change_index(arg,from,to) for arg in arguments(mul)])
change_index(x,from::Index,to::Index) = x

ismergeable(a::IndexedOperator,b::IndexedOperator) = isequal(a.ind,b.ind) ? ismergeable(a.op,b.op) : false

getIndName(op::IndexedOperator) = op.ind.name
getIndName(ind::Index) = ind.name
getIndName(x) = Symbol()

SymbolicUtils.istree(a::SingleSum) = false
SymbolicUtils.arguments(a::SingleSum) = SymbolicUtils.arguments(a.term)
SymbolicUtils.arguments(a::IndexedOperator) = [a]

get_order(::IndexedOperator) = 1
#It is assumed that the term for which this operation is done already commutes with indices inside the indices-Vector
function order_by_index(vec::Vector,indices::Vector{Index})
    vec_ = copy(vec)
    frontfront = filter(x -> !(typeof(x) == IndexedOperator),vec_)
    front = filter(x -> typeof(x) == IndexedOperator && x.ind in indices,vec_)
    back = filter(x -> typeof(x) == IndexedOperator && x.ind ∉ indices,vec_)
    sort!(front,by=getIndName)
    return vcat(frontfront,front,back)
end
order_by_index(qmul::QMul,inds::Vector{Index}) = qmul.arg_c*prod(order_by_index(qmul.args_nc,inds))
order_by_index(avrg::Average,inds::Vector{Index}) = order_by_index(arguments(avrg)[1],inds)
order_by_index(x,inds) = x
#Reorder function: given a tuple vector of indices meaning for each tuple: first ≠ second
#-> go through the term given and exchange 2 ops when the second has "lower" (i.e. its name is first in the alphabet) index than the first one
#-> results in a term, ordered by its commutating indices
"""
    reorder(param,indexMapping)

Reorders a given term (param) regarding to a given indexMapping, which specifies, which [`Index`](@ref) entities can not be equal
inside the given term. reorder() creates a [`SpecialIndexedTerm`](@ref) as a result.

Examples
========

    reorder(σⱼ²¹ * σᵢ²¹,[(i,j)]) = σᵢ²¹ * σⱼ²¹

    reorder(σⱼ²¹ * σᵢ²¹ * σⱼ¹²,[(i,j)]) = σᵢ²¹ * σⱼ²²

"""
function reorder(param::QMul,indexMapping::Vector{Tuple{Index,Index}})
    term = copy(param.args_nc)
    carg = param.arg_c
    indOps = []
    others = []
    for i = 1:length(term) #Split into indexed ops and non indexed ops
        if term[i] isa IndexedOperator
            push!(indOps,term[i])
        else
            push!(others,term[i])
        end
    end
    if isequal(carg,0) || (0 in term)
        return 0
    end
    finish = false
    while !(finish) #go over all ops ind indexed ops -> order by
        finish = true
        for i = 1:(length(indOps)-1)
            if ((indOps[i].ind,indOps[i+1].ind) in indexMapping || (indOps[i+1].ind,indOps[i].ind) in indexMapping) && (indOps[i+1].ind < indOps[i].ind)
                temp = indOps[i+1]
                indOps[i+1] = indOps[i]
                indOps[i] = temp
                finish = false
            end
        end
    end
    args = vcat(others,indOps)
    qmul = carg*prod(args)

    if qmul isa QMul
        mapping_ = orderMapping(indexMapping)
        return SpecialIndexedTerm(qmul,mapping_)
    else
        return reorder(qmul,indexMapping)
    end
end
reorder(sum::SingleSum,indexMapping::Vector{Tuple{Index,Index}}) = SingleSum(reorder(sum.term,indexMapping),sum.sum_index,sum.non_equal_indices)
reorder(op::SpecialIndexedTerm) = reorder(op.term,op.indexMapping)
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
reorder(x::IndexedOperator,indexMapping::Vector{Tuple{Index,Index}}) = SpecialIndexedTerm(x,indexMapping)
reorder(x::SpecialIndexedTerm,indexMapping::Vector) = reorder(x.term,indexMapping)
reorder(x,indMap) = x

function orderMapping(mapping::Vector{Tuple{Index,Index}})
    mapping_ = Vector{Union{Missing,Tuple{Index,Index}}}(missing,length(mapping))
    for i = 1:length(mapping)
        sort_ = sort([first(mapping[i]),last(mapping[i])],by=getIndName)
        mapping_[i] = (sort_[1],sort_[2])
    end
    return mapping_
end

#Show functions
function Base.show(io::IO,op::IndexedOperator)
    op_ = op.op
    if typeof(op_) <:Transition
        write(io,Symbol(op_.name,op_.i,op_.j,op.ind.name))
    elseif op_ isa Destroy
        write(io,Symbol(op_.name,op.ind.name))
    elseif op_ isa Create
        write(io,Symbol(op_.name,op.ind.name,"'"))
    else
        write(io,op_.name)
    end
end
function Base.show(io::IO,indSum::SingleSum)
    write(io, "Σ", "($(indSum.sum_index.name)", "=1:$(indSum.sum_index.range))",)
    if !(isempty(indSum.non_equal_indices))
        write(io,"($(indSum.sum_index.name)≠")
        for i = 1:length(indSum.non_equal_indices)
            write(io, "$(indSum.non_equal_indices[i].name)")
            if i == length(indSum.non_equal_indices)
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
#Functions for easier symbol creation in Constructor
function writeNEIs(neis::Vector{IndexInt})
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

_to_expression(ind::Index) = ind.name
function _to_expression(x::IndexedOperator)
    x.op isa Transition && return :( IndexedOperator($(x.op.name),$(x.ind.name),$(x.op.i),$(x.op.j)) )
    x.op isa Destroy && return :(IndexedDestroy($(x.op.name),$(x.ind.name)))
    x.op isa Create && return :(dagger(IndexedDestroy($(x.op.name),$(x.ind.name))))
end
_to_expression(s::SingleSum) = :( SingleSum($(_to_expression(s.term)),$(s.sum_index.name),$(s.sum_index.range),$(writeNEIs(s.non_equal_indices))))
_to_expression(a::SymbolicUtils.Sym{Parameter,IndexedVariable}) = :(IndexedVariable($(a.metadata.name),$(a.metadata.ind.name)))
_to_expression(a::SymbolicUtils.Sym{Parameter,DoubleIndexedVariable}) = :(DoubleIndexedVariable($(a.metadata.name),$(a.metadata.ind1.name),$(a.metadata.ind2.name)))

@latexrecipe function f(s::SingleSum)
    neis = writeNEIs(s.non_equal_indices)

    ex = latexify(s.term)
    sumString = nothing
    if neis != ""
        sumString = L"$\underset{%$(s.sum_index.name) ≠%$(neis) }{\overset{%$(s.sum_index.range)}{\sum}}$ %$(ex)"
    else
        sumString = L"$\underset{%$(s.sum_index.name)}{\overset{%$(s.sum_index.range)}{\sum}}$ %$(ex)"
    end
    return sumString
end
SymbolicUtils._iszero(x::SpecialIndexedTerm) = SymbolicUtils._iszero(x.term)
get_range(i::Index) = i.range
get_aon(i::Index) = i.aon

