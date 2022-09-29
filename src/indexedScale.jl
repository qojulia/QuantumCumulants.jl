#file for functions regarding scaling
"""
    scaleME(me::IndexedMeanfieldEquations)

Function, that evaluates a given [`MeanfieldEquations`](@ref) entity and returns again equations,
where indices have been inserted and sums evaluated, regarding the same relations, as done when calculating
with oparators using a [`ClusterSpace`](@ref).

# Arguments
*`me::MeanfieldEquations`: A [`MeanfieldEquations`](@ref) entity, which shall be evaluated.

"""
function scaleME(me::IndexedMeanfieldEquations)
    newEqs = []
    for eq in me.equations
        tempEq = scaleEq(eq)
        if tempEq.lhs in getLHS.(newEqs)
            continue
        elseif isNotIn(getOps(tempEq.lhs;scaling=true),getOps.(getLHS.(newEqs);scaling=true),true) && isNotIn(getOps(_conj(tempEq.lhs);scaling=true),getOps.(getLHS.(newEqs);scaling=true),true)
            push!(newEqs,tempEq)
        end
    end
    filter!(x -> !isequal(x,missing), newEqs)
    vs = getLHS.(newEqs)
    varmap = make_varmap(vs, me.iv)
    ops = undo_average.(vs)
    return IndexedMeanfieldEquations(newEqs,me.operator_equations,vs,ops,me.hamiltonian,me.jumps,me.jumps_dagger,me.rates,me.iv,varmap,me.order)
end
scaleTerm(sym::SymbolicUtils.Sym{Parameter,IndexedAverageSum}) = scaleTerm(sym.metadata)
function scaleTerm(sum::IndexedAverageSum)
    prefact = sum.sumIndex.rangeN - length(sum.nonEqualIndices)
    term_ = scaleTerm(sum.term)
    return prefact*term_
end
function scaleTerm(pow::SymbolicUtils.Pow)
    args = arguments(pow)
    return scaleTerm(args[1])^(args[2])
end
function scaleTerm(add::SymbolicUtils.Add)
    args = arguments(add)
    adds = Vector{Any}(nothing,length(args))
    for i=1:length(args)
        adds[i] = scaleTerm(args[i])
    end
    filter!(x -> !=(x,nothing),adds)
    return sum(adds)
end
function scaleTerm(mul::SymbolicUtils.Mul)
    mults = []
    for arg in arguments(mul)
        push!(mults,scaleTerm(arg))
    end
    if length(mults) == 1
        return mults[1]
    elseif isempty(mults)
        return 0
    else
        return *(mults...)
    end
end
function scaleTerm(x::Average)
    indices = getIndices(x)
    newterm = x
    for i=1:length(indices)
        newterm = insertIndex(newterm,indices[i],i)
    end
    return newterm
end
scaleTerm(x::IndexedOperator) = NumberedOperator(x.op,1)
function scaleTerm(x::QMul)
    indices = getIndices(x)
    term = x
    for i = 1:length(indices)
        term = insertIndex(term,indices[i],i)
    end
    return term
end
scaleTerm(x::SymbolicUtils.Sym{Parameter,SpecialIndexedAverage}) = scaleTerm(x.metadata.term) #this is fine, since the intrinsic conditions on the indices go away automatically from the scaling 
scaleEq(eq::Symbolics.Equation) = Symbolics.Equation(scaleTerm(eq.lhs),scaleTerm(eq.rhs))
scaleTerm(x) = x

SymbolicUtils.substitute(avrg::SymbolicUtils.Sym{Parameter,SpecialIndexedAverage},subs;fold=false) = SymbolicUtils.substitute(avrg.metadata.term,subs;fold=fold)#SpecialIndexedAverage(SymbolicUtils.substitute(term.metadata.term,subs;fold=fold),term.metadata.indexMapping)
SymbolicUtils.substitute(sum::SymbolicUtils.Sym{Parameter,IndexedAverageSum},subs;fold=false) = IndexedAverageSum(SymbolicUtils.substitute(sum.metadata.term,subs;fold=fold),sum.metadata.sumIndex,sum.metadata.nonEqualIndices)
#function to split sums into different clusters, keeping 
#assume: all the indices that are not equal to the summation index are in the same Sum
#for anything else than the assumption, there needs an extra argument, to specify where or how the extra indices are handled
"""
    splitSums(term::SymbolicUtils.Symbolic,amount::Union{<:SymbolicUtils.Sym,<:Int64})
    splitSums(me::MeanfieldEquations,amount)

Function, that splits sums inside a given term. The sums are split into a number of equal sums, specified in the `amount` argument,
where in only one of the sums the dependencies for the indices (non equal indices) is considered.

# Arguments
*`me::MeanfieldEquations`: A [`MeanfieldEquations`](@ref) entity, which shall be evaluated, can also be any symbolic expression.
*`amount::Union{<:SymbolicUtils.Sym,<:Int64}`: A Number or Symbolic determining, in how many terms a sum is split

"""
function splitSums(term::SymbolicUtils.Symbolic,ind::Index,amount::Union{<:SymbolicUtils.Sym,<:Int64})
    if term isa Average
        return term
    end
    if SymbolicUtils.istree(term)
        args = []
        op = SymbolicUtils.operation(term)
        for i = 1:length(arguments(term))
            push!(args,splitSums(arguments(term)[i],ind,amount))
        end
        if isempty(args)
            return 0
        elseif length(args) == 1
            return args[1]
        end
        return op(args...)
    end
    if typeof(term) == SymbolicUtils.Sym{Parameter,IndexedAverageSum}
        term_ = term.metadata.term
        sumInd = term.metadata.sumIndex
        if isequal(ind,sumInd)
            ind2 = Index(sumInd.hilb,sumInd.name,(term.metadata.sumIndex.rangeN/amount),sumInd.specHilb)
            extrasum = IndexedAverageSum(term_,ind,term.metadata.nonEqualIndices)
            return extrasum + (amount-1)*IndexedAverageSum(term_,ind,Index[])
        end
    end
    return term
end
splitSums(x::Symbolics.Equation,ind::Index,amount) = Symbolics.Equation(x.lhs,splitSums(x.rhs,ind,amount))
function splitSums(me::AbstractMeanfieldEquations,ind::Index,amount)
    newEqs = Vector{Union{Missing,Symbolics.Equation}}(missing,length(me.equations))
    for i=1:length(me.equations)
        newEqs[i] = splitSums(me.equations[i],ind,amount)
    end
    newEqs = filter(x -> !isequal(x,missing), newEqs)
    vs = getLHS.(newEqs)
    varmap = make_varmap(vs, me.iv)
    return IndexedMeanfieldEquations(newEqs,me.operator_equations,vs,me.operators,me.hamiltonian,me.jumps,me.jumps_dagger,me.rates,me.iv,varmap,me.order)
end
splitSums(x,ind,amount) = x
#TODO: specify also for double sums (maybe ?)

scale(eqs::IndexedMeanfieldEquations;kwargs...) = scaleME(eqs;kwargs...)

#Some utility functions -> implemented here, since this file is the last to get imported
"""
    createMap(ps::Vector,p0::Vector)

A Function to create parameter values for indexed Variables more convenient.

# Arguments
*`ps::Vector`: A vector of parameters, that have no value assigned to them.
*`p0::Vector`: A vector for numeric values, that should get assigned to the corresponding
    entry in the `ps` vector. For Single-Indexed Variables the entry in the vector can also be again
    a Vector, that has an amount of entries as the index of the variables has range. For Double-Indexed
    Variables, this can also be a Matrix of a dimension, that corresponds to the ranges of the indices
    of the given variable.

"""
function createMap(ps::Vector,p0::Vector)
    length(ps) != length(p0) && error("Vectors given have non-equal length!")

    dict = Dict{Sym{Parameter, Base.ImmutableDict{DataType, Any}},ComplexF64}()
    for i=1:length(ps)
        dicVal = nothing
        if ps[i] isa SymbolicUtils.Sym{Parameter, IndexedVariable}
            if p0[i] isa Vector || p0[i] isa Number
                dicVal = createValueMap(ps[i],p0[i])
            else
                error("cannot resolve entry at $i-th position in values-vector")
            end
        elseif ps[i] isa SymbolicUtils.Sym{Parameter, DoubleIndexedVariable}
            if p0[i] isa Matrix || p0[i] isa Number
                dicVal = createValueMap(ps[i],p0[i])
            end
        else
            push!(dict,ps[i]=>p0[i])
            continue
        end
        dict = merge(dict,dicVal)
    end
    return collect(dict)
end

