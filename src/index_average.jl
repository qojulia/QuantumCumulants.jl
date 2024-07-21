#Main file for manipulating indexed averages and sums over averages.
using ModelingToolkit

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


const symbolics_terms = Union{<:Average,<:BasicSymbolic{<:CNumber}}
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
* sum_index: The index, for which the summation will go over.
* non_equal_indices: (optional) A vector of indices, for which the summation-index can not be equal with.

"""
struct IndexedAverageSum <: CNumber
    term::symbolics_terms
    sum_index::Index
    non_equal_indices::Vector{IndexInt}
    metadata
    function IndexedAverageSum(term::symbolics_terms,sum_index::Index,non_equal_indices::Vector,metadata)
        neis_sym = ""
        if !(isempty(non_equal_indices))
            neis_sym = string("(",neis_sym)
            neis_sym = string(neis_sym, "$(sum_index.name)≠")
            neis_sym = string(neis_sym, writeNEIs(non_equal_indices))
            neis_sym = string(neis_sym,")")
        end
        _metadata = new(term,sum_index,non_equal_indices,metadata)
        sym = SymbolicUtils.Sym{IndexedAverageSum}(Symbol("∑($(sum_index.name)=1:$(sum_index.range))$(neis_sym)$(term)"))
        sym = SymbolicUtils.setmetadata(sym,typeof(_metadata),_metadata)
        sym = SymbolicUtils.setmetadata(sym,typeof(metadata),metadata)
        return sym
    end
end
function IndexedAverageSum(term::symbolics_terms,sum_index::Index,non_equal_indices::Vector;metadata=NO_METADATA)
    if sum_index ∉ get_indices(term)
        return (sum_index.range - length(non_equal_indices)) * term
    end
    prefact = 1.0
    if iscall(term)
        op = operation(term)
        args = arguments(term)
        # move numbers outside of sum
        if op === *
            if args[1] isa Number
                prefact = args[1]
                args_nc = args[2:end]
                if length(args_nc) == 1
                    term = args_nc[1]
                else
                    term = *(args_nc...)
                end
            end
        end
        if op === +
            return sum(IndexedAverageSum(arg,sum_index,non_equal_indices;metadata=metadata) for arg in arguments(term))
        end
    end
    return prefact*IndexedAverageSum(term,sum_index,non_equal_indices,metadata)
end
IndexedAverageSum(x,args...;kwargs...) = average(SingleSum(x,args...;kwargs...))
IndexedAverageSum(x::Number) = x

"""

    IndexedAverageDoubleSum <: CNumber

Defines a symbolic summation over an [`IndexedAverageSum`](@ref), using a [`Index`](@ref) entity. This schematically represent a double-sum over a multiplication of averages.

Fields:
======

* innerSum: An [`IndexedAverageSum`](@ref) entity.
* sum_index: The index, for which the (outer) summation will go over.
* non_equal_indices: (optional) A vector of indices, for which the (outer) summation-index can not be equal with.

"""
struct IndexedAverageDoubleSum <: CNumber
    innerSum::BasicSymbolic{IndexedAverageSum}
    sum_index::Index
    non_equal_indices::Vector{IndexInt}
    function IndexedAverageDoubleSum(term::BasicSymbolic{IndexedAverageSum},sum_index::Index,non_equal_indices)
        _metadata = new(term,sum_index,non_equal_indices)
        neis_sym = ""
        if !(isempty(non_equal_indices))
            neis_sym = string("(",neis_sym)
            neis_sym = string(neis_sym, "$(sum_index.name)≠")
            neis_sym = string(neis_sym, writeNEIs(non_equal_indices))
            neis_sym = string(neis_sym,")")
        end
        sym = SymbolicUtils.Sym{IndexedAverageDoubleSum}(Symbol("∑($(sum_index.name):=1:$(sum_index.range))$(neis_sym)$(String(term.name))"))
        sym = SymbolicUtils.setmetadata(sym,typeof(_metadata),_metadata)
        return sym
    end
end
function IndexedAverageDoubleSum(term::symbolics_terms,sum_index::Index,non_equal_indices)
    if iscall(term)
        op = operation(term)
        args = arguments(term)
        param = 1.0
        if op === *
            if args[1] isa Number #put numbers out in front
                param = args[1]
                deleteat!(args,1)
            end
            if length(args) == 1 && args[1] isa BasicSymbolic{IndexedAverageSum}
                return param*IndexedAverageDoubleSum(args[1],sum_index,non_equal_indices)
            end
        elseif op === +
            return sum(IndexedAverageDoubleSum(arg,sum_index,non_equal_indices) for arg in arguments(term))
        end
    end
    return IndexedAverageSum(term,sum_index,non_equal_indices)
end
IndexedAverageDoubleSum(x,y,z) = IndexedAverageSum(x,y,z)

#For representing in average terms
"""

    NumberedOperator <: QSym

Defines an operator, associated with a Number. Commutator-relations are calculated using these numbers, as a sort of a specific index-value.

Fields:
======

* op: An Operator, either a [`Transition`](@ref), a [`Destroy`](@ref) or a [`Create`](@ref) can be defined.
* numb: An Integer Number.

"""
struct NumberedOperator <:QSym
    op::IndexableOps
    numb::Int64
    function NumberedOperator(op::IndexableOps,numb::Int64)
        (numb <= 0) && error("can not create numbered-operator with negative or 0 number")
        return new(op,numb)
    end
end
function NumberedOperator(op,numb::Int64)
    (op isa SNuN) && return op
    if SymbolicUtils.iscall(op)
        f = SymbolicUtils.operation(op)
        args = [NumberedOperator(arg,numb) for arg in SymbolicUtils.arguments(op)]
        isempty(args) && return 0
        isequal(length(args),1) && return args[1]
        return f(args...)
    end
end

#Variables
"""

    SingleNumberedVariable <: numberedVariable

Defines a variable, associated with a Number. Used by [`value_map`](@ref)

Fields:
======

* name: The name of the variable.
* numb: An Integer Number.

"""
struct SingleNumberedVariable <: numberedVariable
    name::Symbol
    numb::Int64
    function SingleNumberedVariable(name,numb)
        sym_name = Symbol("$(name)_$(numb)")
        return Parameter(sym_name)
    end
end
"""

    DoubleNumberedVariable <: numberedVariable

Defines a variable, associated with two Numbers. Used by [`value_map`](@ref)

Fields:
======

* name: The name of the variable.
* numb1: An Integer Number.
* numb2: Another Integer Number.

"""
struct DoubleNumberedVariable <: numberedVariable
    name::Symbol
    numb1::IndexInt
    numb2::IndexInt
    function DoubleNumberedVariable(name,numb1,numb2;identical::Bool=true)
        if !(identical) && (numb1 == numb2)
            return 0
        end
        if typeof(numb1) == typeof(numb2) && numb1 isa Int64
            sym_name = Symbol("$(name)_{$(numb1),$(numb2)}")
        return Parameter(sym_name)
        else
            metadata = new(name,numb1,numb2)
            sym = SymbolicUtils.Sym{DoubleNumberedVariable}(Symbol("$(name)_{$(numb1),$(numb2)}"))
            sym = SymbolicUtils.setmetadata(sym,typeof(metadata),metadata)
            return sym
        end
    end
end
struct SpecialIndexedAverage <: CNumber #An average-Term with special condition, for example l ≠ k; needed for correct calculus of Double indexed Sums
    term::symbolics_terms
    indexMapping::Vector{Tuple{IndexInt,IndexInt}}
    function SpecialIndexedAverage(term::Average,indexMapping)
        if isempty(indexMapping)
            return term
        end
        if SymbolicUtils._iszero(arguments(term)[1])
            return 0
        end
        metadata = new(term,indexMapping)
        neis = writeIndexNEIs(indexMapping)
        sym = SymbolicUtils.Sym{SpecialIndexedAverage}(Symbol("$(neis)$(term)"))
        sym = SymbolicUtils.setmetadata(sym,typeof(metadata),metadata)
        return sym
    end
end
function SpecialIndexedAverage(term::symbolics_terms,indexMapping)
    if iscall(term)
        op = operation(term)
        args = arguments(term)
        if op === *
            prefac = 1
            if args[1] isa Number
                prefac = args[1]
                deleteat!(args,1)
            end
            if length(args) == 1
                return prefac * SpecialIndexedAverage(args[1],indexMapping)
            end
            specInds = [SpecialIndexedAverage(arg,indexMapping) for arg in args]
            return prefac * prod(specInds)
        elseif op === +
            return sum(SpecialIndexedAverage(arg,indexMapping) for arg in args)
        end
    end
    return term
end
SpecialIndexedAverage(x,args...) = x

const AvgSums = Union{BasicSymbolic{IndexedAverageSum},BasicSymbolic{IndexedAverageDoubleSum},BasicSymbolic{SpecialIndexedAverage},IndexedAverageSum,IndexedAverageDoubleSum,SpecialIndexedTerm}
const AvgS = Union{IndexedAverageSum,IndexedAverageDoubleSum,SpecialIndexedTerm}

average(indOp::IndexedOperator) = SymbolicUtils._iszero(indOp) ? 0 : _average(indOp)
average(x::SpecialIndexedTerm) = SpecialIndexedAverage(average(x.term),x.indexMapping)
average(indSum::SingleSum; kwargs...) = IndexedAverageSum(average(indSum.term),indSum.sum_index,indSum.non_equal_indices)
average(indDSum::DoubleSum) = IndexedAverageDoubleSum(average(indDSum.innerSum),indDSum.sum_index,indDSum.NEI)

undo_average(a::IndexedAverageSum) = SingleSum(undo_average(a.term),a.sum_index,a.non_equal_indices)
undo_average(a::IndexedAverageDoubleSum) = DoubleSum(undo_average(a.innerSum),a.sum_index,a.non_equal_indices)
undo_average(a::SpecialIndexedAverage) = reorder(undo_average(a.term),a.indexMapping)
function undo_average(a::BasicSymbolic{IndexedAverageSum})
    if SymbolicUtils.hasmetadata(a,IndexedAverageSum)
        meta = SymbolicUtils.metadata(a)[IndexedAverageSum]
        return undo_average(meta)
    end
end
function undo_average(a::BasicSymbolic{IndexedAverageDoubleSum})
    if SymbolicUtils.hasmetadata(a,IndexedAverageDoubleSum)
        meta = SymbolicUtils.metadata(a)[IndexedAverageDoubleSum]
        return undo_average(meta)
    end
end
function undo_average(a::BasicSymbolic{SpecialIndexedAverage})
    if SymbolicUtils.hasmetadata(a,SpecialIndexedAverage)
        meta = SymbolicUtils.metadata(a)[SpecialIndexedAverage]
        return undo_average(meta)
    end
end


#define calculus for numbered operators -> break it down into QNumber multiplication

ismergeable(a::NumberedOperator,b::NumberedOperator) = isequal(a.numb,b.numb) ? ismergeable(a.op,b.op) : false

*(numOp::NumberedOperator, qmul::QMul) = inorder!(QMul(qmul.arg_c,vcat(numOp,qmul.args_nc)))
*(qmul::QMul, numOp::NumberedOperator) = inorder!(QMul(qmul.arg_c,vcat(qmul.args_nc,numOp)))

function *(numOp1::NumberedOperator,numOp2::NumberedOperator)
    if numOp1.op isa Create || numOp1.op isa Destroy || numOp2.op isa Create || numOp2.op isa Destroy
        return merge_commutators(1,[numOp1,numOp2])
    end
    return (numOp1.numb == numOp2.numb && isequal(acts_on(numOp1.op),acts_on(numOp2.op))) ? NumberedOperator(numOp1.op*numOp2.op,numOp1.numb) : QMul(1,sort([numOp1,numOp2],by=get_numbers))
end
#Symbolics functions
get_order(a::IndexedAverageSum) = get_order(a.term)
get_order(a::IndexedAverageDoubleSum) = get_order(a.innerSum)
get_order(a::SpecialIndexedAverage) = get_order(a.term)
function get_order(a::BasicSymbolic{IndexedAverageSum})
    if SymbolicUtils.hasmetadata(a,IndexedAverageSum)
        meta = SymbolicUtils.metadata(a)[IndexedAverageSum]
        return get_order(meta)
    end
end
function get_order(a::BasicSymbolic{IndexedAverageDoubleSum})
    if SymbolicUtils.hasmetadata(a,IndexedAverageDoubleSum)
        meta = SymbolicUtils.metadata(a)[IndexedAverageDoubleSum]
        return get_order(meta)
    end
end
function get_order(a::BasicSymbolic{SpecialIndexedAverage})
    if SymbolicUtils.hasmetadata(a,SpecialIndexedAverage)
        meta = SymbolicUtils.metadata(a)[SpecialIndexedAverage]
        return get_order(meta)
    end
end

SymbolicUtils._iszero(a::IndexedAverageSum) = SymbolicUtils._iszero(a.term)
SymbolicUtils._isone(a::IndexedAverageSum) = SymbolicUtils._isone(a.term)

SymbolicUtils.iscall(a::IndexedAverageSum) = false
SymbolicUtils.iscall(a::BasicSymbolic{SpecialIndexedAverage}) = false
SymbolicUtils.iscall(a::IndexedAverageDoubleSum) = false
SymbolicUtils.iscall(a::BasicSymbolic{IndexedAverageDoubleSum}) = false


average(x::NumberedOperator) = _average(x)
hilbert(x::NumberedOperator) = hilbert(x.op)
Base.adjoint(x::NumberedOperator) = NumberedOperator(Base.adjoint(x.op),x.numb)
has_cluster(x::NumberedOperator) = has_cluster(x.op)
acts_on(x::NumberedOperator) = acts_on(x.op)
get_order(x::NumberedOperator) = get_order(x.op)

function writeIndexNEIs(neis::Vector{Tuple{IndexInt,IndexInt}})
    syms = ""
    syms = join([syms,"("])
    for i = 1:length(neis)
        if first(neis[i]) isa Index
            syms = join([syms,first(neis[i]).name])
        else
            syms = join([syms,first(neis[i])])
        end
        syms = join([syms,"≠"])
        if last(neis[i]) isa Index
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
writeIndexNEIs(neis::Vector{Tuple{Index,Index}}) = writeIndexNEIs(convert(Vector{Tuple{IndexInt,IndexInt}},neis))
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
    return hash(IndexedAverageSum, hash(a.term, hash(a.sum_index, hash(a.non_equal_indices,h))))
end
function Base.hash(a::NumberedOperator,h::UInt)
    return hash(NumberedOperator, hash(a.op, hash(a.numb, h)))
end
Base.isequal(x::NumberedOperator,y::NumberedOperator) = isequal(x.op,y.op) && isequal(x.numb,y.numb)
Base.isless(a::IndexedAverageSum,b::IndexedAverageSum) = a.sum_index < b.sum_index
function Base.isequal(a::BasicSymbolic{IndexedAverageSum},b::BasicSymbolic{IndexedAverageSum})
    a_meta = SymbolicUtils.metadata(a)[IndexedAverageSum]
    b_meta = SymbolicUtils.metadata(b)[IndexedAverageSum]
    return isequal(a_meta,b_meta)
end
function Base.isequal(a::IndexedAverageSum, b::IndexedAverageSum)
    isequal(a.sum_index,b.sum_index) || return false
    isequal(a.term, b.term) || return false
    isequal(a.non_equal_indices,b.non_equal_indices) || return false
    return true
end
function Base.isequal(a::BasicSymbolic{SpecialIndexedAverage},b::BasicSymbolic{SpecialIndexedAverage})
    a_meta = SymbolicUtils.metadata(a)[SpecialIndexedAverage]
    b_meta = SymbolicUtils.metadata(b)[SpecialIndexedAverage]
    return isequal(a_meta.term,b_meta.term) && isequal(a_meta.indexMapping,b_meta.indexMapping)
end

function _cumulant_expansion(x::IndexedAverageSum,order;kwargs...)
    return IndexedAverageSum(simplifyMultiplication(cumulant_expansion(x.term,order;kwargs...)),x.sum_index,x.non_equal_indices)
end
function _cumulant_expansion(x::IndexedAverageDoubleSum,order;kwargs...)
    inner = _cumulant_expansion(x.innerSum,order;kwargs...)
    return IndexedAverageDoubleSum(inner,x.sum_index,x.non_equal_indices)
end
function _cumulant_expansion(a::BasicSymbolic{IndexedAverageSum},order;kwargs...)
    if SymbolicUtils.hasmetadata(a,IndexedAverageSum)
        meta = SymbolicUtils.metadata(a)[IndexedAverageSum]
        return _cumulant_expansion(meta,order;kwargs...)
    end
end
function _cumulant_expansion(a::BasicSymbolic{IndexedAverageDoubleSum},order;kwargs...)
    if SymbolicUtils.hasmetadata(a,IndexedAverageDoubleSum)
        meta = SymbolicUtils.metadata(a)[IndexedAverageDoubleSum]
        return _cumulant_expansion(meta,order;kwargs...)
    end
end
function _cumulant_expansion(a::BasicSymbolic{SpecialIndexedAverage},order;kwargs...)
    if SymbolicUtils.hasmetadata(a,SpecialIndexedAverage)
        meta = SymbolicUtils.metadata(a)[SpecialIndexedAverage]
        return SpecialIndexedAverage(cumulant_expansion(meta.term,order;kwargs...),meta.indexMapping)
    end
end

SymbolicUtils.arguments(op::BasicSymbolic{IndexedAverageSum}) = arguments(SymbolicUtils.metadata(op)[IndexedAverageSum])
SymbolicUtils.arguments(op::IndexedAverageSum) = arguments(op.term)
SymbolicUtils.arguments(op::BasicSymbolic{IndexedAverageDoubleSum}) = arguments(SymbolicUtils.metadata(op)[IndexedAverageDoubleSum])
SymbolicUtils.arguments(op::IndexedAverageDoubleSum) = op.innerSum
SymbolicUtils.arguments(op::BasicSymbolic{SpecialIndexedAverage}) = arguments(SymbolicUtils.metadata(op)[SpecialIndexedAverage])
SymbolicUtils.arguments(op::SpecialIndexedAverage) = arguments(op.term)

#this is the new method, insert values directly into the average before calculating anything, simplifies evaluation afterwards extremely
#function for inserting index, k -> 1,2,...,N
"""
    insert_index(term,ind::Index,value::Int)

Function, that inserts an integer value for a index in a specified term.
This function creates Numbered- Variables/Operators/Sums upon calls.

Examples
========

    insert_index(σⱼ²¹,j,1) = σ₁²¹

"""
function insert_index(sum::BasicSymbolic{IndexedAverageSum}, ind::Index, value::Int64)
    meta = SymbolicUtils.metadata(sum)[IndexedAverageSum]
    if ind == meta.sum_index
        error("cannot exchange summation index with number!")
    end
    if ind in meta.non_equal_indices
        newNEI = filter(x-> !isequal(x,ind),meta.non_equal_indices)
        push!(newNEI,value)
        return IndexedAverageSum(insert_index(meta.term,ind,value),meta.sum_index,newNEI)
    else
        return IndexedAverageSum(insert_index(meta.term,ind,value),meta.sum_index,meta.non_equal_indices)
    end
end
function insert_index(sum::BasicSymbolic{IndexedAverageDoubleSum}, ind::Index,value::Int64)
    meta = SymbolicUtils.metadata(sum)[IndexedAverageDoubleSum]
    inner = insert_index(meta.innerSum,ind,value)
    return IndexedAverageDoubleSum(inner,meta.sum_index,meta.non_equal_indices)
end
function insert_index(term::BasicSymbolic{<:CNumber},ind::Index,value::Int64)
    if iscall(term)
        op = operation(term)
        if op === *
            return prod(insert_index(arg,ind,value) for arg in arguments(term))
        elseif op === +
            return sum(insert_index(arg,ind,value) for arg in arguments(term))
        elseif op === ^
            return insert_index(arguments(term)[1],ind,value)^insert_index(arguments(term)[2],ind,value)
        # issue 198
        elseif op === /
            return insert_index(arguments(term)[1],ind,value)/insert_index(arguments(term)[2],ind,value)
        elseif length(arguments(term)) == 1 # exp, sin, cos, ln, ... #TODO: write tests
            return op(insert_index(arguments(term)[1],ind,value))
        end
    end
    return term
end
function insert_index(term::Average,ind::Index,value::Int64)
    f = operation(term)
    if f == conj
        return conj(insert_index(arguments(term)[1],ind,value))
    end
    return average(inorder!(insert_index(arguments(term)[1],ind,value)))
end
function insert_index(term_::BasicSymbolic{DoubleIndexedVariable},ind::Index,value::Int64)
    term = SymbolicUtils.metadata(term_)[DoubleIndexedVariable]
    if term.ind1 == ind && term.ind2 == ind
        return DoubleNumberedVariable(term.name,value,value)
    elseif term.ind1 == ind
        return DoubleNumberedVariable(term.name,value,term.ind2)
    elseif term.ind2 == ind
        return DoubleNumberedVariable(term.name,term.ind1,value)
    end
    return term_
end
function insert_index(term::BasicSymbolic{DoubleNumberedVariable},ind::Index,value::Int64)
    if iscall(term)
        op = operation(term)
        if op === *
            return prod(insert_index(arg,ind,value) for arg in arguments(term))
        elseif op === +
            return sum(insert_index(arg,ind,value) for arg in arguments(term))
        elseif op === ^
            return insert_index(arguments(term)[1],ind,value)^(arguments(term)[2])
        end
    end
    data = SymbolicUtils.metadata(term)[DoubleNumberedVariable]
    if data.numb1 isa Index && data.numb1 == ind
        return DoubleNumberedVariable(data.name,value,data.numb2)
    elseif data.numb2 isa Index && data.numb2 == ind
        return DoubleNumberedVariable(data.name,data.numb1,value)
    end
    return term
end
function insert_index(term::BasicSymbolic{SpecialIndexedAverage},ind::Index,value::Int64)
    meta = SymbolicUtils.metadata(term)[SpecialIndexedAverage]
    newterm = insert_index(meta.term,ind,value)
    newMapping = Tuple{IndexInt,IndexInt}[]
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
    filter!(x -> !(first(x) isa Int64 && last(x) isa Int64),newMapping)
    return SpecialIndexedAverage(newterm,newMapping)
end
insert_index(qmul::QMul,ind::Index,value::Int64) = qmul.arg_c*prod(insert_index(arg,ind,value) for arg in qmul.args_nc)
insert_index(eq::Symbolics.Equation,ind::Index,value::Int64) = Symbolics.Equation(insert_index(eq.lhs,ind,value),insert_index(eq.rhs,ind,value))
insert_index(term::IndexedOperator,ind::Index,value::Int64) = term.ind == ind ? NumberedOperator(term.op,value) : term
function insert_index(term::BasicSymbolic{IndexedVariable},ind::Index,value::Int64)
    meta = SymbolicUtils.metadata(term)[IndexedVariable]
    meta.ind == ind ? SingleNumberedVariable(meta.name,value) : term
end
insert_index(x,args...) = x
"""
    insert_indices(eq::Symbolics.Equation,map::Dict{Index,Int64};limits=Dict{SymbolicUtils.BasicSymbolic,Int64}())

Function, that inserts an integer value for a index in a specified Equation. This function creates Numbered- Variables/Operators/Sums upon calls.
Mainly used by [`evalEquation`](@ref).

# Arguments
*`eq::Symbolics.Equation`: The equation for which the indices will be inserted.
*`map::Dict{Index,Int64}`: A dictionary, which contains specifications for the insertions
    the entry (i => 5) would result in all `i` indices being replaced with the number 5.

# Optional argumentes
*`limits::Dict{SymbolicUtils.BasicSymbolic,Int64}=Dict{Symbol,Int64}()`: A seperate dictionary, to
    specify any symbolic limits used when [`Index`](@ref) entities were defined. This needs
    to be specified, when the equation contains summations, for which the upper bound is given
    by a Symbolic.

"""
function insert_indices(eq::Symbolics.Equation,map::Dict{Index,Int64};limits=Dict{SymbolicUtils.BasicSymbolic,Int64}(),kwargs...)
    eq_rhs = eq.rhs
    while !isempty(map)
        pair = first(map)
        eq_rhs = insert_index(eq_rhs,first(pair),last(pair))
        delete!(map,first(pair))
    end
    return eval_term(eq_rhs;limits,kwargs...) #return finished equation
end
function insert_indices_lhs(term::Average,map::Dict{Index,Int64};kwargs...)
    lhs = term
    map_ = copy(map)
    while !isempty(map_)
        pair = first(map_)
        lhs = insert_index(lhs,first(pair),last(pair))
        delete!(map_,first(pair))
        inorder!(lhs)
    end
    return lhs
end
"""
    evalME(me::MeanfieldEquations;limits::Dict{SymbolicUtils.BasicSymbolic,Int64}=Dict{SymbolicUtils.BasicSymbolic,Int64}())

Function, that evaluates a given [`MeanfieldEquations`](@ref) entity and returns again equations,
where indices have been inserted and sums evaluated.

# Arguments
*`me::MeanfieldEquations`: A [`MeanfieldEquations`](@ref) entity, which shall be evaluated.

# Optional argumentes
*`limits=Dict{SymbolicUtils.BasicSymbolic,Int64}()`: A seperate dictionary, to
    specify any symbolic limits used when [`Index`](@ref) entities were defined. This needs
    to be specified, when the equations contain summations, for which the upper bound is given
    by a Symbolic.

"""
function evalME(me::AbstractMeanfieldEquations;limits=Dict{SymbolicUtils.BasicSymbolic,Int64}(),h=nothing,kwargs...)#this is still pretty slow
    vs = me.states
    maxRange = count_eq_number(vs;limits=limits,h=h,kwargs...)
    if !(maxRange isa Int)
        error("Not all upper limits of indices are set as a Number! You can do this by using the \"limits\" keyword argument.")
    end
    newEqs = Vector{Any}(nothing,maxRange)
    states = Vector{Any}(nothing,maxRange)
    counter = 1
    for i = 1:length(vs)
        inds = get_indices(vs[i])
        eq = me.equations[i]
        if !=(h,nothing)
            filter!(x->x.aon in h,inds)
        end
        if isempty(inds)
            eval = evalEq(eq;limits=limits,h=h,kwargs...)
            if (eval.lhs ∉ states) && (_inconj(eval.lhs) ∉ states)
                states[counter] = eval.lhs
                newEqs[counter] = eval
                counter = counter + 1
            end
        else
            if !=(h,nothing)
                filter!(x->x.aon in h,inds)
            end
            ranges_ = Vector{Any}(nothing,length(inds))
            for i=1:length(inds)
                ranges_[i] = (inds[i].range in keys(limits)) ? (1:limits[inds[i].range]) : (1:inds[i].range)
            end
            arr = create_index_arrays(inds,ranges_)
            for j=1:length(arr)
                if !isempty(get_numbers(eq.lhs)) && !(check_arr(eq.lhs,arr[j]))
                    continue
                end
                dict = Dict{Index,Int}(inds .=> arr[j])
                eq_lhs = insert_indices_lhs(eq.lhs,dict)
                if (eq_lhs ∉ states) && (_inconj(eq_lhs) ∉ states)
                    eq_rhs = insert_indices(eq,dict;limits=limits,h=h,kwargs...)
                    states[counter] = eq_lhs
                    if SymbolicUtils._iszero(eq_rhs)
                        newEqs[counter] = Symbolics.Equation(eq_lhs,0)
                    else
                        newEqs[counter] = Symbolics.Equation(eq_lhs,eq_rhs)
                    end
                    counter = counter + 1
                end
            end
        end
    end
    states = states[1:(counter-1)]
    operats = undo_average.(states)
    newEqs = newEqs[1:(counter-1)]
    varmap = make_varmap(states, me.iv)
    return EvaledMeanfieldEquations(newEqs,me.operator_equations,states,operats,me.hamiltonian,me.jumps,me.jumps_dagger,me.rates,me.iv,varmap,me.order)
 end
 # function that counts how many equations are needed for a given set of states
 function count_eq_number(vs;limits=Dict(),h=nothing,kwargs...)
    if !=(h,nothing) && !(h isa Vector)
        h = [h]
    end
    counter = 0
    for state in vs
        inds = get_indices(state)
        if !=(h,nothing)
            filter!(x->x.aon in h, inds)
        end
        if isempty(inds)
            counter = counter + 1
        else
            ranges = get_range.(inds)
            counter = counter + prod(ranges)
        end
    end
    return substitute(counter,limits)
end
function eval_term(sum_::BasicSymbolic{IndexedAverageSum};limits=Dict{SymbolicUtils.BasicSymbolic,Int64}(), h=nothing, kwargs...)
    meta = SymbolicUtils.metadata(sum_)[IndexedAverageSum]
    if !=(h,nothing)
        if !(h isa Vector)
            h = [h]
        end
        meta.sum_index.aon ∉ h && return sum_
    end
    rangeEval = 0
    if meta.sum_index.range in keys(limits)
        rangeEval = limits[meta.sum_index.range]
    else
        if meta.sum_index.range isa BasicSymbolic{<:CNumber}
            if iscall(meta.sum_index.range)
                args = arguments(meta.sum_index.range)
                args_ = Vector{Any}(nothing,length(args))
                for i=1:length(args)
                    if args[i] in keys(limits)
                        args_[i] = limits[args[i]]
                    else
                        args_[i] = args[i]
                    end
                end
                if operation(meta.sum_index.range) === *
                    rangeEval = prod(args_)
                end
            else
                rangeEval = meta.sum_index.range
            end
        else
            rangeEval = meta.sum_index.range
        end
    end
    adds = Vector{Any}(nothing,rangeEval)
    for i = 1:rangeEval
        if i in meta.non_equal_indices
            adds[i] = 0
        else
            temp = insert_index(meta.term,meta.sum_index,i)
            inorder!(temp)
            adds[i]=temp
        end
    end
    if isempty(adds)
        return 0
    end
    return sum(adds)
end
function eval_term(sum::BasicSymbolic{IndexedAverageDoubleSum};kwargs...)
    meta = SymbolicUtils.metadata(sum)[IndexedAverageDoubleSum]
    return eval_term(IndexedAverageDoubleSum(eval_term(meta.innerSum;kwargs...),meta.sum_index,meta.non_equal_indices);kwargs...)
end
function eval_term(term::BasicSymbolic{<:CNumber};kwargs...)
    if iscall(term)
        op = operation(term)
        if op === +
            return sum(eval_term(arg;kwargs...) for arg in arguments(term))
        end
        if op === *
            return prod(eval_term(arg;kwargs...) for arg in arguments(term))
        end
        # issue 198 #TODO: tests
        if op === ^
            args = arguments(term)
            return eval_term(args[1];kwargs...)^eval_term(args[2];kwargs...)
        end
        if op === /
            args = arguments(term)
            return eval_term(args[1];kwargs...)/eval_term(args[2];kwargs...)
        end

        if length(arguments(term)) == 1 # exp, sin, cos, ln, ...
            return op(eval_term(arguments(term)[1];kwargs...))
        end
    end
    return term
end

function eval_term(x;kwargs...)
    inorder!(x)
    return x
end
function evalEq(eq::Symbolics.Equation;kwargs...)
    rhs_ = eval_term(eq.rhs;kwargs...)
    if SymbolicUtils._iszero(rhs_)
        return Symbolics.Equation(eq.lhs,0)
    else
        return Symbolics.Equation(eq.lhs,rhs_)
    end
end

getLHS(eq::Symbolics.Equation) = eq.lhs
getLHS(x) = []

Base.:(==)(term1::Average,term2::Average) = isequal(arguments(term1), arguments(term2))

#Value map creation, for easier inserting into the ODEProblem
"""
    create_value_map(sym::BasicSymbolic{IndexedVariable}, values::Vector;limits::Dict{SymbolicUtils.BasicSymbolic,Int64}=Dict{SymbolicUtils.BasicSymbolic,Int64}())
    create_value_map(sym::BasicSymbolic{IndexedVariable}, value::Number)
    create_value_map(sym::BasicSymbolic{DoubleIndexedVariable},values::Matrix;limits::Dict{SymbolicUtils.BasicSymbolic,Int64}=Dict{SymbolicUtils.BasicSymbolic,Int64}())

Function, that creates a Dictionary, for which a indexedVariable is associated with a series of (number) values. The dictionary contains Symbols of either [`SingleNumberedVariable`](@ref)
or [`DoubleNumberedVariables`](@ref) as keys and the values as values. For a Single-indexed variable, one can
create such a limits by giving a Vector of values, and for double-indexed variables by giving a Matrix. One can also create such a limits, by using
only a single value, then all possible numbered-Variables are set to the same values.

# Arguments
*`sym`: Either a [`IndexedVariable`](@ref) or a [`DoubleIndexedVariable`](@ref)
*`values`: For a [`IndexedVariable`](@ref) either a vector or a single number, and for [`DoubleIndexedVariable`](@ref) a matrix.

# Optional argumentes
*`limits::Dict{SymbolicUtils.BasicSymbolic,Int64}=Dict{BasicSymbolic,Int64}()`: A seperate dictionary, to
    specify any symbolic limits used when [`Index`](@ref) entities were defined. This needs
    to be specified, when the equations contain summations, for which the upper bound is given
    by a Symbolic.

"""
function create_value_map(sym::BasicSymbolic{IndexedVariable}, values::Vector;limits=Dict{BasicSymbolic,Int64}(),kwargs...)
    iVar = SymbolicUtils.metadata(sym)[IndexedVariable]
    if iVar.ind.range isa SymbolicUtils.BasicSymbolic
        if iVar.ind.range in keys(limits)
            range1 = limits[iVar.ind.range]
        else
            error("Can not evaluate without a limits")
        end
    else
        range1 = iVar.ind.range
    end
    if range1 != length(values)
        error("different length of index-range and given values!")
    end
    dict = Dict{BasicSymbolic,ComplexF64}()
    for i = 1:range1
        push!(dict,(SingleNumberedVariable(iVar.name,i) => values[i]))
    end
    return dict
end
function create_value_map(sym::BasicSymbolic{IndexedVariable}, value::Number;limits=Dict{BasicSymbolic,Int64}(),kwargs...)
    iVar = SymbolicUtils.metadata(sym)[IndexedVariable]
    dict = Dict{BasicSymbolic,ComplexF64}()
    if iVar.ind.range isa SymbolicUtils.BasicSymbolic
        if iVar.ind.range in keys(limits)
            range1 = limits[iVar.ind.range]
        else
            error("Can not evaluate without a limits")
        end
    else
        range1 = iVar.ind.range
    end
    for i = 1:range1
        push!(dict,(SingleNumberedVariable(iVar.name,i) => value))
    end
    return dict
end
function create_value_map(sym::BasicSymbolic{DoubleIndexedVariable},values::Matrix;limits=Dict{BasicSymbolic,Int64}(),kwargs...)
    dict = Dict{BasicSymbolic,ComplexF64}()
    var = SymbolicUtils.metadata(sym)[DoubleIndexedVariable]
    if var.ind1.range isa SymbolicUtils.BasicSymbolic
        if var.ind1.range in keys(limits)
            range1 = limits[var.ind1.range]
        else
            error("Can not evaluate without a limits")
        end
    else
        range1 = var.ind1.range
    end
    if var.ind2.range isa SymbolicUtils.BasicSymbolic
        if var.ind2.range in keys(limits)
            range2 = limits[var.ind2.range]
        else
            error("Can not evaluate without a limits")
        end
    else
        range2 = var.ind2.range
    end
    for i = 1:range1
        for j = 1:range2
            push!(dict,(DoubleNumberedVariable(var.name,i,j) => values[i,j]))
        end
    end
    return dict
end

#functions for checking if indices occure in specific terms
function containsIndexedOps(term::Average)
    arg_ = arguments(term)
    if arg_[1] isa QMul
        for arg in arg_[1].args_nc
            if arg isa IndexedOperator
                return true
            end
        end
    else
        return arg_[1] isa IndexedOperator
    end
    return false
end
containsIndex(term::Average,ind::Index) = ind ∈ get_indices(term)

function SymbolicUtils.simplify(sym::BasicSymbolic{SpecialIndexedAverage})
    meta = SymbolicUtils.metadata(sym)[SpecialIndexedAverage]
    SpecialIndexedAverage(SymbolicUtils.simplify(meta.term),meta.indexMapping)
end

#function that creates an array consisting of all possible number values for each index given
#ind_vec should be sorted beforehand
function create_index_arrays(ind_vec,ranges)
    if length(ind_vec) == 1
        return ranges[1]
    end
    @assert length(ind_vec) == length(ranges)
    array = unique(collect(Iterators.product(ranges...)))
    length(ind_vec) == 1 && return array
    length(unique(get_spec_hilb.(ind_vec))) == length(ind_vec) && return array #every index has its own specHilb
    for vec in get_not_allowed(ind_vec)
        array = array[Bool[all_different(array[i],vec) for i=1:length(array)]]
    end
    return collect(array)
end
all_different(x,vec) = length(unique(getindex(x,vec))) == length(getindex(x,vec))
function get_not_allowed(ind_vec)
    spec_hilbs = get_spec_hilb.(ind_vec)
    not_allowed = []
    for ind in ind_vec
        indices = findall(x -> isequal(x,ind.aon) ,spec_hilbs)
        length(indices) == 1 && continue
        if indices ∉ not_allowed
            push!(not_allowed,indices)
        end
    end
    return not_allowed
end
get_spec_hilb(ind::Index) = ind.aon

function check_arr(lhs,arr)
    numbs = get_numbers(lhs)
    inds = get_indices(lhs)
    D = Dict(inds.=>arr)
    args_ = arguments(lhs)[1]
    if args_ isa QMul
        args = args_.args_nc
    else
        args = [args_]
    end
    for i = 1:length(hilbert(args[1]).spaces)
        as = filter(x->isequal(acts_on(x),i),args)
        isempty(get_numbers(as)) && continue
        isempty(get_indices(as)) && continue
        inds_ = get_indices(as)
        numbs = get_numbers(as)
        for i in inds_
            if D[i] in numbs
                return false
            end
        end
    end
    return true
end

getAvrgs(sum::BasicSymbolic{SpecialIndexedAverage}) = getAvrgs(SymbolicUtils.metadata(sum)[SpecialIndexedAverage].term)
getAvrgs(sum::BasicSymbolic{IndexedAverageSum}) = getAvrgs(SymbolicUtils.metadata(sum)[IndexedAverageSum].term)
getAvrgs(Dsum::BasicSymbolic{IndexedAverageDoubleSum}) = getAvrgs(SymbolicUtils.metadata(Dsum)[IndexedAverageDoubleSum].innerSum)
function getAvrgs(term::BasicSymbolic{<:CNumber})
    if iscall(term)
        return  vcat(filter(x->!=(x,nothing),[getAvrgs(arg) for arg in arguments(term)])...)
    else
        return nothing
    end
end
getAvrgs(avrg::Average) = avrg
getAvrgs(x) = nothing

getNumber(x::NumberedOperator) = [acts_on(x) + x.numb]
getNumber(x::QMul) = acts_on(x)
getNumber(x) = [acts_on(x)] # this is so that, any other operator still behaves the same as before

function Base.show(io::IO,indSum::IndexedAverageSum)
    write(io, "Σ", "($(indSum.sum_index.name)", "=1:$(indSum.sum_index.range))",)
    if !(isempty(indSum.non_equal_indices))
        write(io,"($(indSum.sum_index.name) ≠ ")
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
function Base.show(io::IO,indSum::IndexedAverageDoubleSum)
    write(io, "Σ", "($(indSum.sum_index.name)", "=1:$(indSum.sum_index.range))",)
    if !(isempty(indSum.non_equal_indices))
        write(io,"($(indSum.sum_index.name) ≠ ")
        for i = 1:length(indSum.non_equal_indices)
            write(io, "$(indSum.non_equal_indices[i].name)")
            if i == length(indSum.non_equal_indices)
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

function _to_expression(x::NumberedOperator)
    x.op isa Transition && return :( NumberedOperator($(x.op.name),$(x.numb),$(x.op.i),$(x.op.j)) )
    x.op isa Destroy && return :(NumberedDestroy($(x.op.name),$(x.numb)))
    x.op isa Create && return :(dagger(NumberedDestroy($(x.op.name),$(x.numb))))
end
function _to_expression(x::BasicSymbolic{IndexedAverageSum})
    meta = SymbolicUtils.metadata(x)[IndexedAverageSum]
    return :( IndexedAverageSum($(_to_expression(meta.term)),$(meta.sum_index.name),$(meta.sum_index.range),$(writeNEIs(meta.non_equal_indices))) )
end
function _to_expression(x::BasicSymbolic{SpecialIndexedAverage})
    meta = SymbolicUtils.metadata(x)[SpecialIndexedAverage]
    return _to_expression(meta.term)
end
function _to_expression(x::BasicSymbolic{IndexedAverageDoubleSum})
    meta = SymbolicUtils.metadata(x)[IndexedAverageDoubleSum]
    return :( IndexedAverageDoubleSum($(_to_expression(meta.innerSum)),$(meta.sum_index.name),$(meta.sum_index.range),$(writeNEIs(meta.non_equal_indices))) )
end

@latexrecipe function f(s_::BasicSymbolic{IndexedAverageSum})
    s = SymbolicUtils.metadata(s_)[IndexedAverageSum]
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

#simplify functions not "really" needed, they are nice to have, since equation creation of Symbolics sometimes does not simplify certain terms
#function to reduce multiplication of numbers with a sum into just a sum of multiplication
function simplifyMultiplication(term::BasicSymbolic{<:CNumber})
    if iscall(term) && operation(term) === *
        args = arguments(term)
        ind = findfirst(x-> (iscall(x) && operation(x) === +),args)
        (ind === nothing) && return term #no add-terms were found inside the multiplication

        args_ = arguments(args[ind]) # arguments of the addition
        lefts = isempty(args[1:(ind-1)]) ? 1 : args[1:(ind-1)]
        rights = isempty(args[(ind+1):end]) ? 1 : args[(ind+1):end]
        adds = [simplifyMultiplication(prod(lefts)*arg*prod(rights)) for arg in args_]

        return sum(adds)
    end
    return term
end
simplifyMultiplication(x) = x


function +(a::BasicSymbolic{SpecialIndexedAverage},b::BasicSymbolic{SpecialIndexedAverage})
    if isequal(a,b)
        return SymbolicUtils.Add(CNumber,0,Dict(a=>2))
    end
    return SymbolicUtils.Add(CNumber,0,Dict(a=>1,b=>1))
end
function +(a::BasicSymbolic{IndexedAverageDoubleSum},b::BasicSymbolic{IndexedAverageDoubleSum})
    if isequal(a,b)
        return SymbolicUtils.Add(CNumber,0,Dict(a=>2))
    end
    return SymbolicUtils.Add(CNumber,0,Dict(a=>1,b=>1))
end
function +(a::BasicSymbolic{IndexedAverageSum},b::BasicSymbolic{IndexedAverageSum})
    if isequal(a,b)
        return SymbolicUtils.Add(CNumber,0,Dict(a=>2))
    end
    return SymbolicUtils.Add(CNumber,0,Dict(a=>1,b=>1))
end
function *(a::BasicSymbolic{SpecialIndexedAverage},b::BasicSymbolic{SpecialIndexedAverage})
    if isequal(a,b)
        return SymbolicUtils.Mul(CNumber,1,Dict(a=>2))
    end
    return SymbolicUtils.Mul(CNumber,1,Dict(a=>1,b=>1))
end
function +(a::BasicSymbolic{IndexedVariable},b::BasicSymbolic{IndexedVariable})
    if isequal(a,b)
        return SymbolicUtils.Add(CNumber,0,Dict(a=>2))
    end
    return SymbolicUtils.Add(CNumber,0,Dict(a=>1,b=>1))
end
function +(a::BasicSymbolic{DoubleIndexedVariable},b::BasicSymbolic{DoubleIndexedVariable})
    if isequal(a,b)
        return SymbolicUtils.Add(CNumber,0,Dict(a=>2))
    end
    return SymbolicUtils.Add(CNumber,0,Dict(a=>1,b=>1))
end
