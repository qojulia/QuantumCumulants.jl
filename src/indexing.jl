#Main class for indexing, here indices, sums and indexed operators and variables are defined with the corresponding calculus.
#Many helping functions used in the different classes are also defined here.

const Ranges = Union{<:SymbolicUtils.Sym,<:Number,<:SymbolicUtils.Mul,<:SymbolicUtils.Div,<:BasicSymbolic} #possible Types for the range of an index
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
* range: The upper bound limit of the index. This can be a SymbolicUtils.Symbolic or any Number.
* aon: Number specifying the specific [`HilbertSpace`](@ref), where the Index acts on.

"""
struct Index <: SymbolicUtils.Symbolic{Int}
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

acts_on(i::Index) = i.aon

function Base.:(==)(i::Index, j::Index)
    if acts_on(i) != acts_on(j) || hilbert(i) != hilbert(j)
        return false
    end

    return TermInterface.maketerm(SymbolicUtils.BasicSymbolic{Bool}, ==, [i, j], nothing)
end

# TermInterface for Index
TermInterface.iscall(::Index) = false
TermInterface.metadata(::Index) = nothing

# const IndexInt = Union{<:Index,<:Int64}
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
        sym = SymbolicUtils.Sym{IndexedVariable}(Symbol("$(name)$(ind.name)"))
        sym = SymbolicUtils.setmetadata(sym,typeof(metadata),metadata)
        return sym
    end
end

# TODO: deprecate IndexedVariable and DoubelIndexedVariable

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
        if !(identical) && isequal(ind1, ind2)
            return 0
        end
        metadata = new(name,ind1,ind2,identical)
        sym = SymbolicUtils.Sym{DoubleIndexedVariable}(Symbol("$(name)$(ind1.name)$(ind2.name)"))
        sym = SymbolicUtils.setmetadata(sym,typeof(metadata),metadata)
        return sym
    end
end

struct IndexedParameterSym <: CNumber end
const SymbolicIndexedParameter = SymbolicUtils.BasicSymbolic{<:IndexedParameterSym}
const sym_idx_parameter = begin
    T = SymbolicUtils.FnType{Tuple{Complex{Real}, Vector{<:Index}}, IndexedParameterSym}
    SymbolicUtils.Sym{T}(:getindex_parameter)
end
SymbolicUtils.symtype(::T) where T <: SymbolicIndexedParameter = IndexedParameterSym

function IndexedParameter(name::Symbol, indices::Vector{Index})
    return SymbolicUtils.Term{IndexedParameterSym}(sym_idx_parameter, [Parameter(name), indices])
end
IndexedParameter(name::Symbol, index::Index) = IndexedParameter(name, [index])
IndexedParameter(name::Symbol, indices::Index...) = IndexedParameter(name, [indices...])
IndexedParameter(name::Symbol) = (args...) -> IndexedParameter(name, args...)

# TermInterface
# TermInterface.iscall(::IndexedParameter) = true
# TermInterface.arguments(x::IndexedParameter) = x.indices
# TermInterface.operation(x::IndexedParameter) = IndexedParameter(x.name)
# TermInterface.maketerm(::Type{<:IndexedParameter}, ::typeof(IndexedParameter), args, metadata) = IndexedParameter(args..., metadata)

# function Base.isequal(pi::IndexedParameter, pj::IndexedParameter)
#     if !isequal(pi.name, pj.name)
#         return false
#     end

#     if !isequal(length(pi.indices), length(pj.indices))
#         return false
#     end

#     for (i, j) in zip(pi.indices, pj.indices)
#         isequal(i, j) || return false
#     end

#     return true
# end

"""
    IndexedOperator <: QSym
    IndexedOperator(op::Union{Transition,Create,Destroy},ind::Index)

Operator, associated with an index.

Fields:
======

* op: Operator, either a [`Transition`](@ref), a [`Destroy`](@ref) or a [`Create`](@ref) can be defined.
* ind: The index, the operator will be associated with.

"""
struct IndexedOperator{T,I} <: QSym
    op::T
    ind::I  # TODO: why not spell out index here?
    merge_events::Vector{UUID}
    function IndexedOperator(op::T,ind::I) where {T<:QSym,I}
        if I <: Index
            @assert isequal(ind.hilb,hilbert(op))
        end
        isa(ind.hilb, ProductSpace) && (@assert isequal(acts_on(op),ind.aon))
        return new{T,I}(op,ind,UUID[])
    end
end

#hilberts
hilbert(ind::Index) = ind.hilb
hilbert(op::IndexedOperator) = op.ind.hilb
hilbert(var::IndexedVariable) = var.ind.hilb

#Basic functions for indexed Operators
import Base: *, +, -

#Multiplications
function Base.:*(a::IndexedOperator{<:Destroy}, b::IndexedOperator{<:Create})
    check_hilbert(a,b)
    aon_a = acts_on(a)
    aon_b = acts_on(b)
    if aon_a == aon_b
        if isequal(a.ind, b.ind)
            # shortcut simplification here
            return b*a + 1
        else
            return b*a + a.ind == b.ind
        end
    elseif aon_a < aon_b
        return QMul(1, [a,b])
    else
        return QMul(1, [b,a])
    end
end

function Base.:*(a::IndexedOperator{<:Transition}, b::IndexedOperator{<:Transition})
    if was_merged(a, b)
        # case when they have been merged, but weren't sorted
        # due to the check in ismergeable, this can only mean that a has more merge events
        # and hence should be moved to the right
        return QMul(1, [b,a])
    end

    check_hilbert(a, b)
    aon_a = acts_on(a)
    aon_b = acts_on(b)
    if aon_a == aon_b
        i = a.ind
        j = b.ind
        op_ = a.op * b.op
        t1 = iszero(op_) ? 0 : IndexedOperator(a.op * b.op, i)
        if isequal(i, j)
            return t1
        else
            a_copy = deepcopy(a)
            b_copy = deepcopy(b)

            # assign a unique id to signal that these have been merged already in the same "merge event"
            id = uuid4()
            push!(a_copy.merge_events, id)
            push!(b_copy.merge_events, id)

            # the operator with more merge events goes right
            if length(a_copy.merge_events) > length(b_copy.merge_events)
                t2 = QMul(1, [b_copy, a_copy])
            else
                t2 = QMul(1, [a_copy, b_copy])
            end
            return t1 * (i == j) + (1 - (i == j)) * t2
        end
    elseif aon_a < aon_b
        return QMul(1, [a,b])
    else
        return QMul(1, [b,a])
    end
end

ismergeable(a::IndexedOperator{<:Destroy}, b::IndexedOperator{<:Create}) = true
function ismergeable(a::IndexedOperator{<:Transition}, b::IndexedOperator{<:Transition})
    if length(a.merge_events) > length(b.merge_events)
        # hop into the merge function if the expression is not yet sorted by length of merge events
        # this is so that operators with fewer merge events go left
        # that's required to ensure every operator will be combined with any other operator in longer products
        return true
    end

    return !was_merged(a, b)
end

function was_merged(a::IndexedOperator, b::IndexedOperator)
    # if they were merged, then a unique id (UUID4) was added to the merge events vector of each of the operators
    # hence, we can just check if they share any common element
    return !isdisjoint(a.merge_events, b.merge_events)
end

#acts on
acts_on(op::IndexedOperator) = acts_on(op.op)

#extra commutators
#Indexed operators, evaluate the commutator directly, if 2 indexed ops have the same index
function commutator(op1::IndexedOperator,op2::IndexedOperator)
    commutated_op = IndexedOperator(commutator(op1.op,op2.op),op1.ind)
    if isequal(op1.ind, op2.ind)
        return commutated_op
    else
        return (op1.ind == op2.ind) * commutated_op
    end
end
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
# Base.isless(a::IndexedOperator, b::IndexedOperator) = a.op.name < b.op.name
# Base.isless(a::QMul, b::QMul) = isless(a.args_nc, b.args_nc)
# Base.isless(a::IndexedOperator,b::QSym) = a.op.name < b.name
# Base.isless(a::QSym,b::IndexedOperator) = a.name < b.op.name
# Base.isless(nothing,b::Symbol) = true
# Base.isless(b::Symbol,nothing) = false

# Base.isless(a::Index,b::Index) = a.name < b.name
# Base.isless(a::SingleSum,b::SingleSum) = Base.isless(a.sum_index,b.sum_index)

Base.isequal(ind1::Index,ind2::Index) = (ind1.name === ind2.name) && isequal(ind1.range,ind2.range) && (ind1.hilb == ind2.hilb) && isequal(ind1.aon,ind2.aon)

Base.isequal(op1::IndexedOperator,op2::IndexedOperator) = isequal(op1.op,op2.op) && isequal(op1.ind,op2.ind)

#Function that changes the index of a sum term into a different indexed term
#used for evaluating the extra terms when multiplying a sum with an operator with different index
#return a new QMul with indices swapped: from -> to index
"""
    change_index(term,from::Index,to::Index)

Exchanges all occurring indices inside the given term, that are equal to the `from` to the `to` index.

Examples
========

    change_index(σⱼ²¹,j,i) = σᵢ²¹

    change_index(σⱼ²¹ * σᵢ¹²,j,i) = σᵢ²²


"""
change_index(x, from::Index, to::Index) = x

function change_index(t::SymbolicUtils.Symbolic, from::Index, to::Index)
    isequal(from, to) || !TermInterface.iscall(t) && return t

    f = SymbolicUtils.operation(t)
    args = SymbolicUtils.arguments(t)
    return f(change_index(args, from, to)...)
end

function change_index(a::IndexedOperator, from::Index, to::Index)
    if !isequal(a.ind, from)
        return a
    end

    return IndexedOperator(a.op, to)
end

function change_index(v::IndexedVariable, from::Index, to::Index)
    if !isequal(v.ind, from)
        return v
    end
    return IndexedVariable(v.name, to)
end

function change_index(v::DoubleIndexedVariable, from::Index, to::Index)
    if isequal(v.ind1, from) && isequal(v.ind2, from)
        return DoubleIndexedVariable(v.name, to, to; identical=v.identical)
    elseif isequal(v.ind1, from)
        return DoubleIndexedVariable(v.name, to, v.ind2; identical=v.identical)
    elseif isequal(v.ind2, from)
        return DoubleIndexedVariable(v.name, v.ind1, to; identical=v.identical)
    end
    return v
end

# function change_index(p::IndexedParameter, from::Index, to::Index)
#     p_ = deepcopy(p)
#     for i=1:length(p.indices)
#         if isequal(p_.indices[i], from)
#             p_.indices[i] = to
#         end
#     end
#     return p_
# end


function change_index(args::Vector, from::Index, to::Index)
    isequal(from, to) && return args

    return [change_index(arg, from, to) for arg in args]
end

function change_index(i::Index, from::Index, to::Index)
    if isequal(i, from)
        return to
    else
        return i
    end
end


# getIndName(op::IndexedOperator) = op.ind.name
# getIndName(ind::Index) = ind.name
# getIndName(x) = Symbol()

# # SymbolicUtils.iscall(a::SingleSum) = false
# # SymbolicUtils.arguments(a::SingleSum) = SymbolicUtils.arguments(a.term)
# # SymbolicUtils.arguments(a::IndexedOperator) = [a]

# get_order(::IndexedOperator) = 1
# #It is assumed that the term for which this operation is done already commutes with indices inside the indices-Vector
# function order_by_index(vec::Vector,indices::Vector{Index})
#     vec_ = copy(vec)
#     frontfront = filter(x -> !(typeof(x) == IndexedOperator),vec_)
#     front = filter(x -> typeof(x) == IndexedOperator && x.ind in indices,vec_)
#     back = filter(x -> typeof(x) == IndexedOperator && x.ind ∉ indices,vec_)
#     sort!(front,by=getIndName)
#     return vcat(frontfront,front,back)
# end
# order_by_index(qmul::QMul,inds::Vector{Index}) = qmul.arg_c*prod(order_by_index(qmul.args_nc,inds))
# order_by_index(avrg::Average,inds::Vector{Index}) = order_by_index(arguments(avrg)[1],inds)
# order_by_index(x,inds) = x

# #Reorder function: given a tuple vector of indices meaning for each tuple: first ≠ second
# #-> go through the term given and exchange 2 ops when the second has "lower" (i.e. its name is first in the alphabet) index than the first one
# #-> results in a term, ordered by its commutating indices
# """
#     reorder(param,indexMapping)

# Reorders a given term (param) regarding to a given indexMapping, which specifies, which [`Index`](@ref) entities can not be equal
# inside the given term. reorder() creates a [`SpecialIndexedTerm`](@ref) as a result.

# Examples
# ========

#     reorder(σⱼ²¹ * σᵢ²¹,[(i,j)]) = σᵢ²¹ * σⱼ²¹

#     reorder(σⱼ²¹ * σᵢ²¹ * σⱼ¹²,[(i,j)]) = σᵢ²¹ * σⱼ²²

# """
# function reorder(param::QMul,indexMapping::Vector{Tuple{Index,Index}})
#     term = copy(param.args_nc)
#     carg = param.arg_c
#     indOps = []
#     others = []
#     for i = 1:length(term) #Split into indexed ops and non indexed ops
#         if term[i] isa IndexedOperator
#             push!(indOps,term[i])
#         else
#             push!(others,term[i])
#         end
#     end
#     if isequal(carg,0) || (0 in term)
#         return 0
#     end
#     finish = false
#     while !(finish) #go over all ops ind indexed ops -> order by
#         finish = true
#         for i = 1:(length(indOps)-1)
#             if ((indOps[i].ind,indOps[i+1].ind) in indexMapping || (indOps[i+1].ind,indOps[i].ind) in indexMapping) && (indOps[i+1].ind < indOps[i].ind)
#                 temp = indOps[i+1]
#                 indOps[i+1] = indOps[i]
#                 indOps[i] = temp
#                 finish = false
#             end
#         end
#     end
#     args = vcat(others,indOps)
#     qmul = carg*prod(args)

#     if qmul isa QMul
#         mapping_ = orderMapping(indexMapping)
#         return SpecialIndexedTerm(qmul,mapping_)
#     else
#         return reorder(qmul,indexMapping)
#     end
# end
# # reorder(sum::SingleSum,indexMapping::Vector{Tuple{Index,Index}}) = SingleSum(reorder(sum.term,indexMapping),sum.sum_index,sum.non_equal_indices)
# function reorder(term::QAdd,indexMapping::Vector{Tuple{Index,Index}})
#     args = []
#     for arg in arguments(term)
#         push!(args,reorder(arg,indexMapping))
#     end
#     if length(args) == 0
#         return 0
#     end
#     if length(args) == 1
#         return args[1]
#     end
#     return +(args...)
# end
# reorder(x::IndexedOperator,indexMapping::Vector{Tuple{Index,Index}}) = SpecialIndexedTerm(x,indexMapping)
# reorder(x,indMap) = x

# function orderMapping(mapping::Vector{Tuple{Index,Index}})
#     mapping_ = Vector{Union{Missing,Tuple{Index,Index}}}(missing,length(mapping))
#     for i = 1:length(mapping)
#         sort_ = sort([first(mapping[i]),last(mapping[i])],by=getIndName)
#         mapping_[i] = (sort_[1],sort_[2])
#     end
#     return mapping_
# end

#Show functions
Base.show(io::IO, index::Index) = write(io, index.name)

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
# function Base.show(io::IO,indSum::SingleSum)
#     write(io, "Σ", "($(indSum.sum_index.name)", "=1:$(indSum.sum_index.range))",)
#     if !(isempty(indSum.non_equal_indices))
#         write(io,"($(indSum.sum_index.name)≠")
#         for i = 1:length(indSum.non_equal_indices)
#             write(io, "$(indSum.non_equal_indices[i].name)")
#             if i == length(indSum.non_equal_indices)
#                 write(io,")")
#             else
#                 write(io,",")
#             end
#         end
#     end
#     Base.show(io,indSum.term)
# end
# function Base.show(io::IO,op::SpecialIndexedTerm)
#     if !isempty(op.indexMapping)
#         Base.write(io,"(")
#     end
#     for i = 1:length(op.indexMapping)
#         Base.write(io,first(op.indexMapping[i]).name)
#         Base.write(io,"≠")
#         Base.write(io,last(op.indexMapping[i]).name)
#         if i != length(op.indexMapping)
#             Base.write(io,";")
#         else
#             Base.write(io,")")
#         end
#     end
#     Base.show(io,op.term)
# end
# #Functions for easier symbol creation in Constructor
# function writeNEIs(neis::Vector{IndexInt})
#     syms = ""
#     for i = 1:length(neis)
#         syms = typeof(neis[i]) == Index ? join([syms,neis[i].name]) : join([syms,neis[i]])
#         if i != length(neis)
#             syms = join([syms,","])
#         end
#     end
#     return syms
# end
# function writeNEIs(neis::Vector{Index})
#     syms = ""
#     for i = 1:length(neis)
#         syms = join([syms,neis[i].name])
#         if i != length(neis)
#             syms = join([syms,","])
#         end
#     end
#     return syms
# end

_to_expression(ind::Index) = ind.name
function _to_expression(x::IndexedOperator)
    x.op isa Transition && return :( IndexedOperator($(x.op.name),$(x.ind.name),$(x.op.i),$(x.op.j)) )
    x.op isa Destroy && return :(IndexedDestroy($(x.op.name),$(x.ind.name)))
    x.op isa Create && return :(dagger(IndexedDestroy($(x.op.name),$(x.ind.name))))
end
# _to_expression(s::SingleSum) = :( SingleSum($(_to_expression(s.term)),$(s.sum_index.name),$(s.sum_index.range),$(writeNEIs(s.non_equal_indices))))
_to_expression(a::IndexedVariable) = :(IndexedVariable($(a.name),$(a.ind.name)))
_to_expression(a::DoubleIndexedVariable) = :(DoubleIndexedVariable($(a.name),$(a.ind1.name),$(a.ind2.name)))
function _to_expression(a::BasicSymbolic{IndexedVariable})
    if SymbolicUtils.hasmetadata(a,IndexedVariable)
        meta = SymbolicUtils.metadata(a)[IndexedVariable]
        return _to_expression(meta)
    end
end
function _to_expression(a::BasicSymbolic{DoubleIndexedVariable})
    if SymbolicUtils.hasmetadata(a,DoubleIndexedVariable)
        meta = SymbolicUtils.metadata(a)[DoubleIndexedVariable]
        return _to_expression(meta)
    end
end

# @latexrecipe function f(s::SingleSum)
#     neis = writeNEIs(s.non_equal_indices)

#     ex = latexify(s.term)
#     sumString = nothing
#     if neis != ""
#         sumString = L"$\underset{%$(s.sum_index.name) ≠%$(neis) }{\overset{%$(s.sum_index.range)}{\sum}}$ %$(ex)"
#     else
#         sumString = L"$\underset{%$(s.sum_index.name)}{\overset{%$(s.sum_index.range)}{\sum}}$ %$(ex)"
#     end
#     return sumString
# end
get_range(i::Index) = i.range
get_aon(i::Index) = i.aon

-(a::BasicSymbolic{IndexedVariable}) = -1*a
