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
* sum_index: The index, for which the summation will go over.
* non_equal_indices: (optional) A vector of indices, for which the summation-index can not be equal with.

"""
struct IndexedAverageSum{M} <: CNumber
    term::symbolics_terms
    sum_index::Index
    non_equal_indices::Vector{IndexInt}
    metadata::M
    function IndexedAverageSum(term::symbolics_terms,sum_index::Index,non_equal_indices::Vector,metadata)
        neis_sym = ""
        if !(isempty(non_equal_indices))
            neis_sym = string("(",neis_sym)
            neis_sym = string(neis_sym, "$(sum_index.name)≠")
            neis_sym = string(neis_sym, writeNEIs(non_equal_indices))
            neis_sym = string(neis_sym,")")
        end
        _metadata = new{typeof(metadata)}(term,sum_index,non_equal_indices,metadata)
        return SymbolicUtils.Sym{Parameter, IndexedAverageSum}(Symbol("∑($(sum_index.name)=1:$(sum_index.range))$(neis_sym)$(term)"), _metadata) #Symbol("∑($(sum_index.name)=1:$(sum_index.range))$(neis_sym)$(term)")
    end
end
function IndexedAverageSum(term::IndexedAdd,sum_index::Index,non_equal_indices::Vector;metadata=NO_METADATA)
    return sum(IndexedAverageSum(arg,sum_index,non_equal_indices;metadata=metadata) for arg in arguments(term))
end
function IndexedAverageSum(term::symbolics_terms,sum_index::Index,non_equal_indices::Vector;metadata=NO_METADATA)
    if sum_index ∉ get_indices(term)
        return (sum_index.range - length(non_equal_indices)) * term
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
    innerSum::Sym{Parameter, IndexedAverageSum}
    sum_index::Index
    non_equal_indices::Vector{IndexInt}
    function IndexedAverageDoubleSum(term::Sym{Parameter, IndexedAverageSum},sum_index::Index,non_equal_indices)
        _metadata = new(term,sum_index,non_equal_indices)
        neis_sym = ""
        if !(isempty(non_equal_indices))
            neis_sym = string("(",neis_sym)
            neis_sym = string(neis_sym, "$(sum_index.name)≠")
            neis_sym = string(neis_sym, writeNEIs(non_equal_indices))
            neis_sym = string(neis_sym,")")
        end
        return SymbolicUtils.Sym{Parameter, IndexedAverageDoubleSum}(Symbol("∑($(sum_index.name):=1:$(sum_index.range))$(neis_sym)$(String(term.name))"), _metadata)
    end
end
function IndexedAverageDoubleSum(term::IndexedAdd,sum_index::Index,non_equal_indices)
    return sum(IndexedAverageDoubleSum(arg,sum_index,non_equal_indices) for arg in arguments(term))
end
function IndexedAverageDoubleSum(term::SymbolicUtils.Mul,sum_index::Index,non_equal_indices)
    args = arguments(term)
    param = 1.0
    if args[1] isa Number #put numbers out in front
        param = args[1]
        deleteat!(args,1)
    end
    if length(args) == 1 && args[1] isa Sym{Parameter, IndexedAverageSum}
        return param*IndexedAverageDoubleSum(args[1],sum_index,non_equal_indices)
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
    if SymbolicUtils.istree(op)
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
        metadata=source_metadata(:Parameter, name)
        s = SymbolicUtils.Sym{Parameter, typeof(metadata)}(Symbol("$(name)_$(numb)"), metadata)
        return SymbolicUtils.setmetadata(s, MTK.MTKParameterCtx, true)
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
            metadata = source_metadata(:Parameter, name)
            s = SymbolicUtils.Sym{Parameter, typeof(metadata)}(Symbol("$(name)_{$(numb1)$(numb2)}"), metadata)
            return SymbolicUtils.setmetadata(s, MTK.MTKParameterCtx, true)
        else
            return SymbolicUtils.Sym{Parameter, numberedVariable}(Symbol("$(name)_{$(numb1)$(numb2)}"), new(name,numb1,numb2))
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
        return SymbolicUtils.Sym{Parameter, SpecialIndexedAverage}(Symbol("$(neis)$(term)"), metadata)
    end
end
function SpecialIndexedAverage(term::SymbolicUtils.Add,indexMapping) 
    sum(SpecialIndexedAverage(arg,indexMapping) for arg in arguments(term))
end
function SpecialIndexedAverage(term::SymbolicUtils.Mul,indexMapping)
    args = arguments(term)
    prefac = 1
    if args[1] isa Number
        prefac = args[1]
        deleteat!(args,1)
    end
    return prefac * prod(SpecialIndexedAverage(arg,indexMapping) for arg in args)
end
SpecialIndexedAverage(x,args...) = x

const AvgSums = Union{SymbolicUtils.Sym{Parameter,IndexedAverageSum},SymbolicUtils.Sym{Parameter,IndexedAverageDoubleSum},SymbolicUtils.Sym{Parameter,SpecialIndexedAverage},IndexedAverageSum,IndexedAverageDoubleSum,SpecialIndexedTerm}

average(indOp::IndexedOperator) = SymbolicUtils._iszero(indOp) ? 0 : _average(indOp)
average(x::SpecialIndexedTerm) = SpecialIndexedAverage(average(x.term),x.indexMapping)
average(indSum::SingleSum; kwargs...) = IndexedAverageSum(average(indSum.term),indSum.sum_index,indSum.non_equal_indices)
average(indDSum::DoubleSum) = IndexedAverageDoubleSum(average(indDSum.innerSum),indDSum.sum_index,indDSum.NEI)

undo_average(a::IndexedAverageSum) = SingleSum(undo_average(a.term),a.sum_index,a.non_equal_indices)
undo_average(a::Sym{Parameter,IndexedAverageSum}) = undo_average(a.metadata)
undo_average(a::Sym{Parameter,IndexedAverageDoubleSum}) = undo_average(a.metadata)
undo_average(a::IndexedAverageDoubleSum) = DoubleSum(undo_average(a.innerSum),a.sum_index,a.non_equal_indices)
undo_average(a::SymbolicUtils.Sym{Parameter,SpecialIndexedAverage}) = reorder(undo_average(a.metadata.term),a.metadata.indexMapping)

#define calculus for numbered operators -> break it down into QNuber multiplication

*(numOp::NumberedOperator, qmul::QMul) = merge_commutators(qmul.arg_c,inorder!(vcat(numOp,qmul.args_nc)))
*(qmul::QMul, numOp::NumberedOperator) = merge_commutators(qmul.arg_c,inorder!(vcat(qmul.args_nc,numOp)))

function *(numOp1::NumberedOperator,numOp2::NumberedOperator)
    if numOp1.op isa Create || numOp1.op isa Destroy || numOp2.op isa Create || numOp2.op isa Destroy
        return merge_commutators(1,[numOp1,numOp2])
    end
    return (numOp1.numb == numOp2.numb && isequal(acts_on(numOp1.op),acts_on(numOp2.op))) ? NumberedOperator(numOp1.op*numOp2.op,numOp1.numb) : inorder!(QMul(1,[numOp1,numOp2]))
end
#Symbolics functions
get_order(a::Sym{Parameter,IndexedAverageSum}) = get_order(a.metadata.term)
get_order(a::Sym{Parameter,IndexedAverageDoubleSum}) = get_order(a.metadata.innerSum)
get_order(a::SymbolicUtils.Sym{Parameter,SpecialIndexedAverage}) = get_order(a.metadata.term)

SymbolicUtils._iszero(a::IndexedAverageSum) = SymbolicUtils._iszero(a.term)
SymbolicUtils._isone(a::IndexedAverageSum) = SymbolicUtils._isone(a.term)

SymbolicUtils.istree(a::IndexedAverageSum) = false
SymbolicUtils.istree(a::SymbolicUtils.Sym{Parameter,SpecialIndexedAverage}) = false
SymbolicUtils.istree(a::IndexedAverageDoubleSum) = false
SymbolicUtils.istree(a::SymbolicUtils.Sym{Parameter,IndexedAverageDoubleSum}) = false


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
Base.isequal(a::SymbolicUtils.Sym{Parameter,IndexedAverageSum},b::SymbolicUtils.Sym{Parameter,IndexedAverageSum}) = isequal(a.metadata,b.metadata)
function Base.isequal(a::IndexedAverageSum, b::IndexedAverageSum)
    isequal(a.sum_index,b.sum_index) || return false
    isequal(a.term, b.term) || return false
    isequal(a.non_equal_indices,b.non_equal_indices) || return false
    return true
end
Base.isequal(a::SymbolicUtils.Sym{Parameter,SpecialIndexedAverage},b::SymbolicUtils.Sym{Parameter,SpecialIndexedAverage}) = isequal(a.metadata.term,b.metadata.term) && isequal(a.metadata.indexMapping,b.metadata.indexMapping)
function Base.isequal(nVal1::Sym{Parameter,numberedVariable},nVal2::Sym{Parameter,numberedVariable})
    if typeof(nVal1) == typeof(nVal2) && nVal1.metadata isa SingleNumberedVariable
        return (nVal1.metadata.name == nVal2.metadata.name) && (nVal1.metadata.numb == nVal2.metadata.numb)
    elseif typeof(nVal1) == typeof(nVal2) && nVal1.metadata isa DoubleNumberedVariable
        return (nVal1.metadata.name == nVal2.metadata.name) && (nVal1.metadata.numb1 == nVal2.metadata.numb1) && (nVal1.metadata.numb2 == nVal2.metadata.numb2)
    end
    return false
end
Base.:(==)(nVal1::Sym{Parameter,numberedVariable},nVal2::Sym{Parameter,numberedVariable}) = (nVal1.name == nVal2.name) && (nVal1.numb == nVal2.numb)

function cumulant_expansion(x::SymbolicUtils.Sym{Parameter,IndexedAverageSum},order::Integer;kwargs...)
    sum = x.metadata
    return IndexedAverageSum(simplifyMultiplication(cumulant_expansion(sum.term,order;kwargs...)),sum.sum_index,sum.non_equal_indices)
end
function cumulant_expansion(x::SymbolicUtils.Sym{Parameter,IndexedAverageDoubleSum},order::Integer;kwargs...)
    inner = cumulant_expansion(x.metadata.innerSum,order;kwargs...)
    return IndexedAverageDoubleSum(inner,x.metadata.sum_index,x.metadata.non_equal_indices)
end
cumulant_expansion(a::SymbolicUtils.Sym{Parameter,SpecialIndexedAverage},order::Int;kwargs...) = SpecialIndexedAverage(cumulant_expansion(a.metadata.term,order;kwargs...),a.metadata.indexMapping)
cumulant_expansion(a::IndexedAverageSum,order::Int) = IndexedAverageSum(simplifyMultiplication(cumulant_expansion(a.term,order)),a.sum_index,a.non_equal_indices)

SymbolicUtils.arguments(op::Sym{Parameter,IndexedAverageSum}) = arguments(op.metadata)
SymbolicUtils.arguments(op::IndexedAverageSum) = arguments(op.term)
SymbolicUtils.arguments(op::Sym{Parameter, IndexedAverageDoubleSum}) = arguments(op.metadata)
SymbolicUtils.arguments(op::IndexedAverageDoubleSum) = op.innerSum
SymbolicUtils.arguments(op::Sym{Parameter,SpecialIndexedAverage}) = arguments(op.metadata)
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
function insert_index(sum::SymbolicUtils.Sym{Parameter,IndexedAverageSum}, ind::Index, value::Int64)
    if ind == sum.metadata.sum_index
        error("cannot exchange summation index with number!")
    end
    if ind in sum.metadata.non_equal_indices
        newNEI = filter(x-> x != ind,sum.metadata.non_equal_indices)
        push!(newNEI,value)
        return IndexedAverageSum(insert_index(sum.metadata.term,ind,value),sum.metadata.sum_index,newNEI)
    else
        return IndexedAverageSum(insert_index(sum.metadata.term,ind,value),sum.metadata.sum_index,sum.metadata.non_equal_indices)
    end
end
function insert_index(sum::SymbolicUtils.Sym{Parameter,IndexedAverageDoubleSum}, ind::Index,value::Int64)
    inner = insert_index(sum.metadata.innerSum,ind,value)
    return IndexedAverageDoubleSum(inner,sum.metadata.sum_index,sum.metadata.non_equal_indices)
end
insert_index(term::SymbolicUtils.Mul, ind::Index, value::Int64) = prod(insert_index(arg,ind,value) for arg in arguments(term))
insert_index(term::SymbolicUtils.Add,ind::Index,value::Int64) = sum(insert_index(arg,ind,value) for arg in arguments(term))
insert_index(term::SymbolicUtils.Pow,ind::Index,value::Int64) = insert_index(arguments(term)[1],ind,value)^(arguments(term)[2])
function insert_index(term::SymbolicUtils.Term{AvgSym,Nothing},ind::Index,value::Int64)
    f = operation(term)
    arg = inorder!(insert_index(arguments(term)[1],ind,value))
    return f(arg)
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
    if term.metadata isa SingleNumberedVariable
        return term
    end
    data = term.metadata
    if data.numb1 isa Index && data.numb1 == ind
        return DoubleNumberedVariable(data.name,value,data.numb2)
    elseif data.numb2 isa Index && data.numb2 == ind
        return DoubleNumberedVariable(data.name,data.numb1,value)
    end
    return term
end
function insert_index(term::SymbolicUtils.Sym{Parameter,SpecialIndexedAverage},ind::Index,value::Int64)
    meta = term.metadata
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
insert_index(term::SymbolicUtils.Sym{Parameter,IndexedVariable},ind::Index,value::Int64) = term.metadata.ind == ind ? SingleNumberedVariable(term.metadata.name,value) : term
insert_index(x,args...) = x
"""
    insert_indices(eq::Symbolics.Equation,map::Dict{Index,Int64};limits=Dict{SymbolicUtils.Sym,Int64}())

Function, that inserts an integer value for a index in a specified Equation. This function creates Numbered- Variables/Operators/Sums upon calls.
Mainly used by [`evalEquation`](@ref).

# Arguments
*`eq::Symbolics.Equation`: The equation for which the indices will be inserted.
*`map::Dict{Index,Int64}`: A dictionary, which contains specifications for the insertions
    the entry (i => 5) would result in all `i` indices being replaced with the number 5.

# Optional argumentes
*`limits::Dict{SymbolicUtils.Sym,Int64}=Dict{Symbol,Int64}()`: A seperate dictionary, to
    specify any symbolic limits used when [`Index`](@ref) entities were defined. This needs
    to be specified, when the equation contains summations, for which the upper bound is given
    by a Symbolic.

"""
function insert_indices(eq::Symbolics.Equation,map::Dict{Index,Int64};limits=Dict{SymbolicUtils.Sym,Int64}(),kwargs...)
    eq_rhs = eq.rhs
    while !isempty(map)
        pair = first(map)
        eq_rhs = insert_index(eq_rhs,first(pair),last(pair))
        delete!(map,first(pair))
    end
    return eval_term(eq_rhs;limits,kwargs...) #return finished equation
end
function insert_indices_lhs(term::Term{AvgSym, Nothing},map::Dict{Index,Int64};kwargs...)
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
    evalME(me::MeanfieldEquations;limits::Dict{SymbolicUtils.Sym,Int64}=Dict{SymbolicUtils.Sym,Int64}())

Function, that evaluates a given [`MeanfieldEquations`](@ref) entity and returns again equations,
where indices have been inserted and sums evaluated.

# Arguments
*`me::MeanfieldEquations`: A [`MeanfieldEquations`](@ref) entity, which shall be evaluated.

# Optional argumentes
*`limits=Dict{SymbolicUtils.Sym,Int64}()`: A seperate dictionary, to
    specify any symbolic limits used when [`Index`](@ref) entities were defined. This needs
    to be specified, when the equations contain summations, for which the upper bound is given
    by a Symbolic.

"""
function evalME(me::AbstractMeanfieldEquations;limits=Dict{SymbolicUtils.Sym,Int64}(),h=nothing,kwargs...)#this is still pretty slow
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
                dict = Dict{Index,Int}(inds .=> arr[j])
                eq_lhs = insert_indices_lhs(eq.lhs,dict)
                if !(eq_lhs in states || _inconj(eq_lhs) in states)
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
    newEqs = newEqs[1:(counter-1)]
    varmap = make_varmap(states, me.iv)
    return EvaledMeanfieldEquations(newEqs,me.operator_equations,states,me.operators,me.hamiltonian,me.jumps,me.jumps_dagger,me.rates,me.iv,varmap,me.order)
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
function eval_term(sum_::SymbolicUtils.Sym{Parameter,IndexedAverageSum};limits=Dict{SymbolicUtils.Sym,Int64}(), h=nothing, kwargs...)
    if !=(h,nothing)
        if !(h isa Vector)
            h = [h]
        end
        sum_.metadata.sum_index.aon ∉ h && return sum_
    end
    rangeEval = 0
    if sum_.metadata.sum_index.range in keys(limits)
        rangeEval = limits[sum_.metadata.sum_index.range]
    else
        if sum_.metadata.sum_index.range isa SymbolicUtils.Mul
            args = arguments(sum_.metadata.sum_index.range)
            args_ = Vector{Any}(nothing,length(args))
            for i=1:length(args)
                if args[i] in keys(limits)
                    args_[i] = limits[args[i]]
                else
                    args_[i] = args[i]
                end
            end
            rangeEval = prod(args_)
        else
            rangeEval = sum_.metadata.sum_index.range
        end
    end
    adds = Vector{Any}(nothing,rangeEval)
    for i = 1:rangeEval
        if i in sum_.metadata.non_equal_indices
            adds[i] = 0
            continue
        end
        temp = insert_index(sum_.metadata.term,sum_.metadata.sum_index,i)
        inorder!(temp)
        adds[i]=temp
    end
    if isempty(adds)
        return 0
    end
    return sum(adds)
end
function eval_term(sum::SymbolicUtils.Sym{Parameter,IndexedAverageDoubleSum};kwargs...)
    return eval_term(IndexedAverageDoubleSum(eval_term(sum.metadata.innerSum;kwargs...),sum.metadata.sum_index,sum.metadata.non_equal_indices);kwargs...)
end
eval_term(term::SymbolicUtils.Mul;kwargs...) = prod(eval_term(arg;kwargs...) for arg in arguments(term))
eval_term(term::SymbolicUtils.Add;kwargs...) = sum(eval_term(arg;kwargs...) for arg in arguments(term))

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

Base.:(==)(term1::Term{AvgSym, Nothing},term2::Term{AvgSym, Nothing}) = isequal(arguments(term1), arguments(term2))

#Value map creation, for easier inserting into the ODEProblem
"""
    create_value_map(sym::Sym{Parameter, IndexedVariable}, values::Vector;limits::Dict{SymbolicUtils.Sym,Int64}=Dict{SymbolicUtils.Sym,Int64}())
    create_value_map(sym::Sym{Parameter, IndexedVariable}, value::Number)
    create_value_map(sym::Sym{Parameter,DoubleIndexedVariable},values::Matrix;limits::Dict{SymbolicUtils.Sym,Int64}=Dict{SymbolicUtils.Sym,Int64}())

Function, that creates a Dictionary, for which a indexedVariable is associated with a series of (number) values. The dictionary contains Symbols of either [`SingleNumberedVariable`](@ref)
or [`DoubleNumberedVariables`](@ref) as keys and the values as values. For a Single-indexed variable, one can
create such a limits by giving a Vector of values, and for double-indexed variables by giving a Matrix. One can also create such a limits, by using
only a single value, then all possible numbered-Variables are set to the same values.

# Arguments
*`sym`: Either a [`IndexedVariable`](@ref) or a [`DoubleIndexedVariable`](@ref)
*`values`: For a [`IndexedVariable`](@ref) either a vector or a single number, and for [`DoubleIndexedVariable`](@ref) a matrix.

# Optional argumentes
*`limits::Dict{SymbolicUtils.Sym,Int64}=Dict{Symbol,Int64}()`: A seperate dictionary, to
    specify any symbolic limits used when [`Index`](@ref) entities were defined. This needs
    to be specified, when the equations contain summations, for which the upper bound is given
    by a Symbolic.

"""
function create_value_map(sym::Sym{Parameter, IndexedVariable}, values::Vector;limits=Dict{SymbolicUtils.Sym,Int64}(),kwargs...)
    iVar = sym.metadata
    if iVar.ind.range isa SymbolicUtils.Sym
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
    dict = Dict{Sym{Parameter, Base.ImmutableDict{DataType, Any}},ComplexF64}()
    for i = 1:range1
        push!(dict,(SingleNumberedVariable(iVar.name,i) => values[i]))
    end
    return dict
end
function create_value_map(sym::Sym{Parameter, IndexedVariable}, value::Number)
    iVar = sym.metadata
    dict = Dict{Sym{Parameter, Base.ImmutableDict{DataType, Any}},ComplexF64}()
    if iVar.ind.range isa SymbolicUtils.Sym
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
function create_value_map(sym::Sym{Parameter,DoubleIndexedVariable},values::Matrix;limits=Dict{SymbolicUtils.Sym,Int64}(),kwargs...)
    dict = Dict{Sym{Parameter, Base.ImmutableDict{DataType, Any}},ComplexF64}()
    var = sym.metadata
    if var.ind1.range isa SymbolicUtils.Sym
        if var.ind1.range in keys(limits)
            range1 = limits[var.ind1.range]
        else
            error("Can not evaluate without a limits")
        end
    else
        range1 = var.ind1.range
    end
    if var.ind2.range isa SymbolicUtils.Sym
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
function containsIndexedOps(term::SymbolicUtils.Term{AvgSym, Nothing})
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
containsIndex(term::SymbolicUtils.Term{AvgSym,Nothing},ind::Index) = ind ∈ get_indices(term)

SymbolicUtils.simplify(sym::SymbolicUtils.Sym{Parameter,SpecialIndexedAverage}) = SpecialIndexedAverage(SymbolicUtils.simplify(sym.metadata.term),sym.metadata.indexMapping)

#function that creates an array consisting of all possible number values for each index given
#ind_vec should be sorted beforehand
function create_index_arrays(ind_vec,ranges)
    if length(ind_vec) == 1
        return ranges[1]
    end
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

getAvrgs(sum::SymbolicUtils.Sym{Parameter,SpecialIndexedAverage}) = getAvrgs(sum.metadata.term)
getAvrgs(sum::SymbolicUtils.Sym{Parameter,IndexedAverageSum}) = getAvrgs(sum.metadata.term)
getAvrgs(term::SymbolicUtils.Mul) = vcat(filter(x -> !=(x,nothing),[getAvrgs(arg) for arg in arguments(term)]))
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
_to_expression(x::SymbolicUtils.Sym{Parameter,IndexedAverageSum}) = :( IndexedAverageSum($(_to_expression(x.metadata.term)),$(x.metadata.sum_index.name),$(x.metadata.sum_index.range),$(writeNEIs(x.metadata.non_equal_indices))) )
_to_expression(x::SymbolicUtils.Sym{Parameter,SpecialIndexedAverage}) = :($(x.metadata.term))
_to_expression(x::SymbolicUtils.Sym{Parameter,IndexedAverageDoubleSum}) = :( IndexedAverageDoubleSum($(_to_expression(x.metadata.innerSum)),$(x.metadata.sum_index.name),$(x.metadata.sum_index.range),$(writeNEIs(x.metadata.non_equal_indices))) )

@latexrecipe function f(s_::SymbolicUtils.Sym{Parameter,IndexedAverageSum})
    s = s_.metadata
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
function simplifyMultiplication(term::SymbolicUtils.Mul)
    args = arguments(term)
    ind = findfirst(x-> x isa SymbolicUtils.Add,args)
    (ind === nothing) && return term #no add-terms were found inside the multiplication

    args_ = arguments(args[ind]) # arguments of the addition
    lefts = isempty(args[1:(ind-1)]) ? 1 : args[1:(ind-1)]
    rights = isempty(args[(ind+1):end]) ? 1 : args[(ind+1):end]
    adds = [simplifyMultiplication(prod(lefts)*arg*prod(rights)) for arg in args_]

    return sum(adds)
end
simplifyMultiplication(x) = x
