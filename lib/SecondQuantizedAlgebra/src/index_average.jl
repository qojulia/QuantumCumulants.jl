#Main file for manipulating indexed averages and sums over averages.

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
        args = copy(arguments(term))
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
        args = copy(arguments(term))
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

# variable
IndexedVariable(x,numb::Int64) = SingleNumberedVariable(x,numb)
IndexedVariable(x,num1::Int64,num2::Int64;kwargs...) = DoubleNumberedVariable(x,num1,num2;kwargs...)
IndexedVariable(name::Symbol,ind1::Index,ind2::Index;kwargs...) = DoubleIndexedVariable(name,ind1,ind2;kwargs...)
DoubleIndexedVariable(x,num1::Int64,num2::Int64;kwargs...) = DoubleNumberedVariable(x,num1,num2;kwargs...)
function get_indices(a::BasicSymbolic{DoubleIndexedVariable})
    meta = SymbolicUtils.metadata(a)[DoubleIndexedVariable]
    return unique([meta.ind1,meta.ind2])
end
get_indices(a::BasicSymbolic{IndexedVariable}) = [SymbolicUtils.metadata(a)[IndexedVariable].ind]

#Symbolics functions
SymbolicUtils._iszero(a::IndexedAverageSum) = SymbolicUtils._iszero(a.term)
SymbolicUtils._isone(a::IndexedAverageSum) = SymbolicUtils._isone(a.term)

SymbolicUtils.iscall(a::IndexedAverageSum) = false
SymbolicUtils.iscall(a::BasicSymbolic{SpecialIndexedAverage}) = false
SymbolicUtils.iscall(a::IndexedAverageDoubleSum) = false
SymbolicUtils.iscall(a::BasicSymbolic{IndexedAverageDoubleSum}) = false


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
end # TODO: only used in tests, remove?

#Base functions
function Base.hash(a::IndexedAverageSum, h::UInt)
    return hash(IndexedAverageSum, hash(a.term, hash(a.sum_index, hash(a.non_equal_indices,h))))
end

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

SymbolicUtils.arguments(op::BasicSymbolic{IndexedAverageSum}) = arguments(SymbolicUtils.metadata(op)[IndexedAverageSum])
SymbolicUtils.arguments(op::IndexedAverageSum) = arguments(op.term)
SymbolicUtils.arguments(op::BasicSymbolic{IndexedAverageDoubleSum}) = arguments(SymbolicUtils.metadata(op)[IndexedAverageDoubleSum])
SymbolicUtils.arguments(op::IndexedAverageDoubleSum) = op.innerSum
SymbolicUtils.arguments(op::BasicSymbolic{SpecialIndexedAverage}) = arguments(SymbolicUtils.metadata(op)[SpecialIndexedAverage])
SymbolicUtils.arguments(op::SpecialIndexedAverage) = arguments(op.term)

function SymbolicUtils.simplify(sym::BasicSymbolic{SpecialIndexedAverage})
    meta = SymbolicUtils.metadata(sym)[SpecialIndexedAverage]
    SpecialIndexedAverage(SymbolicUtils.simplify(meta.term),meta.indexMapping)
end




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
