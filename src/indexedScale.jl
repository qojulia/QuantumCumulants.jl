#file for functions regarding scaling
#include("indexedMeanfield.jl")
"""
    scaleME(me::MeanfieldEquations)

Function, that evaluates a given [`MeanfieldEquations`](@ref) entity and returns again equations,
where indices have been inserted and sums evaluated, regarding the same relations, as done when calculating
with oparators using a [`ClusterSpace`](@ref).

# Arguments
*`me::MeanfieldEquations`: A [`MeanfieldEquations`](@ref) entity, which shall be evaluated.

"""
function scaleME(me::MeanfieldEquations)
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
    return MeanfieldEquations(newEqs,me.operator_equations,vs,ops,me.hamiltonian,me.jumps,me.jumps_dagger,me.rates,me.iv,varmap,me.order)
end
getLHS(eq::Symbolics.Equation) = eq.lhs
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
function splitSums(term::SymbolicUtils.Symbolic,amount::Union{<:SymbolicUtils.Sym,<:Int64})
    if term isa Average
        return term
    end
    if SymbolicUtils.istree(term)
        args = []
        op = SymbolicUtils.operation(term)
        for i = 1:length(arguments(term))
            push!(args,splitSums(arguments(term)[i],amount))
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
        ind = Index(term.metadata.sumIndex.hilb,term.metadata.sumIndex.name,(term.metadata.sumIndex.rangeN/amount))
        extrasum = IndexedAverageSum(term_,ind,term.metadata.nonEqualIndices)
        return extrasum + (amount-1)*IndexedAverageSum(term_,ind,Index[])
    end
    return term
end
splitSums(x::Symbolics.Equation,amount) = Symbolics.Equation(x.lhs,splitSums(x.rhs,amount))
function splitSums(me::MeanfieldEquations,amount)
    newEqs = Vector{Union{Missing,Symbolics.Equation}}(missing,length(me.equations))
    for i=1:length(me.equations)
        newEqs[i] = splitSums(me.equations[i],amount)
    end
    newEqs = filter(x -> !isequal(x,missing), newEqs)
    vs = getLHS.(newEqs)
    varmap = make_varmap(vs, me.iv)
    return MeanfieldEquations(newEqs,me.operator_equations,vs,me.operators,me.hamiltonian,me.jumps,me.jumps_dagger,me.rates,me.iv,varmap,me.order)
end
splitSums(x,amount) = x