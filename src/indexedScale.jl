#file for functions regarding scaling
include("indexedMeanfield.jl")

function scaleME(me::MeanfieldEquations)
    newEqs = Vector{Union{Missing,Symbolics.Equation}}(missing,length(me.equations))
    for i = 1:length(me.equations)
        newEqs[i] = scaleEq(me.equations[i])
    end
    filter!(x -> !isequal(x,missing), newEqs)
    vs = getLHS.(newEqs)
    varmap = make_varmap(vs, me.iv)
    ops = scaleTerm.(me.operators)
    return MeanfieldEquations(newEqs,me.operator_equations,vs,ops,me.hamiltonian,me.jumps,me.jumps_dagger,me.rates,me.iv,varmap,me.order)
end
scaleTerm(sym::SymbolicUtils.Sym{Parameter,IndexedAverageSum}) = scaleTerm(sym.metadata)
function scaleTerm(sum::IndexedAverageSum)
    prefact = sum.sumIndex.rangeN - length(sum.nonEqualIndices)
    term_ = scaleTerm(sum.term)
    return prefact*term_
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

substitute(term::SymbolicUtils.Sym{Parameter,SpecialIndexedAverage},subs;fold=false) = SpecialIndexedAverage(substitute(term.metadata.term,subs;fold=fold),term.metadata.indexMapping)
function substitute(sum::SymbolicUtils.Sym{Parameter,IndexedAverageSum},subs;fold=false)
    mults = []
    for arg in arguments(sum.metadata.term)
        push!(mults,substitute(arg,subs;fold))
    end
    return IndexedAverageSum(*(mults...),sum.metadata.sumIndex,sum.metadata.nonEqualIndices)
end
#function to split sums into different clusters, keeping 
#assume: all the indices that are not equal to the summation index are in the same Sum
#for anything else than the assumption, there needs an extra argument, to specify where or how the extra indices are handled
function splitSums(term::SymbolicUtils.Symbolic,amount::Union{<:SymbolicUtils.Sym,<:Int64})
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
    for i=1:length(eqs)
        newEqs[i] = splitSums(me.equations[i],amount)
    end
    newEqs = filter(x -> !isequal(x,missing), newEqs)
    vs = getLHS.(newEqs)
    varmap = make_varmap(vs, me.iv)
    return MeanfieldEquations(newEqs,me.operator_equations,vs,me.operators,me.hamiltonian,me.jumps,me.jumps_dagger,me.rates,me.iv,varmap,me.order)
end
splitSums(x,amount) = x