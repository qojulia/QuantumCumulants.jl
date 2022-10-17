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
    NEI = filter(x->x in getIndices(sum.term),sum.nonEqualIndices)
    prefact = sum.sumIndex.rangeN - length(NEI)
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
        if arg isa SymbolicUtils.Sym{Parameter,DoubleIndexedVariable}
            # Double numbered variables need extra computation, since they depend on the number of indices within an average
            # for example Γij * <σi> * <σj> -> Γ11, since both averages only have 1 index
            # however Γij <σi * σj> -> Γ12 = Γ21 for scaled terms
            meta = arg.metadata
            inds = []
            for arg in arguments(mul)
                push!(inds,getIndices(arg))
            end
            order = maximum(length.(inds))
            if order <= 1
                push!(mults,DoubleNumberedVariable(meta.name,1,1))
            else
                push!(mults,DoubleNumberedVariable(meta.name,1,2))
            end
        else
            push!(mults,scaleTerm(arg))
        end
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
        newterm = insert_index(newterm,indices[i],i)
    end
    return newterm
end
scaleTerm(x::IndexedOperator) = NumberedOperator(x.op,1)
function scaleTerm(x::QMul)
    indices = getIndices(x)
    term = x
    for i = 1:length(indices)
        term = insert_index(term,indices[i],i)
    end
    return term
end
scaleTerm(x::SymbolicUtils.Sym{Parameter,SpecialIndexedAverage}) = scaleTerm(x.metadata.term) #this is fine, since the intrinsic conditions on the indices go away automatically from the scaling 
scaleEq(eq::Symbolics.Equation) = Symbolics.Equation(scaleTerm(eq.lhs),scaleTerm(eq.rhs))
scaleTerm(sym::SymbolicUtils.Sym{Parameter,IndexedVariable}) = SingleNumberedVariable(sym.metadata.name,1)
scaleTerm(x) = x

SymbolicUtils.substitute(avrg::SymbolicUtils.Sym{Parameter,SpecialIndexedAverage},subs;fold=false) = SymbolicUtils.substitute(avrg.metadata.term,subs;fold=fold)#SpecialIndexedAverage(SymbolicUtils.substitute(term.metadata.term,subs;fold=fold),term.metadata.indexMapping)
function SymbolicUtils.substitute(sum::SymbolicUtils.Sym{Parameter,IndexedAverageSum},subs;fold=false) 
    subTerm = SymbolicUtils.substitute(sum.metadata.term,subs;fold=fold)
    if SymbolicUtils._iszero(subTerm)
        return 0
    elseif subTerm isa symbolics_terms
        return IndexedAverageSum(subTerm,sum.metadata.sumIndex,sum.metadata.nonEqualIndices)
    else
        return (sum.metadata.sumIndex - length(sum.metadata.nonEqualIndices)) * subTerm
    end
end
#function to split sums into different clusters, keeping 
#assume: all the indices that are not equal to the summation index are in the same Sum
#for anything else than the assumption, there needs an extra argument, to specify where or how the extra indices are handled
"""
    split_sums(term::SymbolicUtils.Symbolic,amount::Union{<:SymbolicUtils.Sym,<:Int64})
    split_sums(me::MeanfieldEquations,amount)

Function, that splits sums inside a given term. The sums are split into a number of equal sums, specified in the `amount` argument,
where in only one of the sums the dependencies for the indices (non equal indices) is considered.

# Arguments
*`me::MeanfieldEquations`: A [`MeanfieldEquations`](@ref) entity, which shall be evaluated, can also be any symbolic expression.
*`amount::Union{<:SymbolicUtils.Sym,<:Int64}`: A Number or Symbolic determining, in how many terms a sum is split

"""
function split_sums(term::SymbolicUtils.Symbolic,ind::Index,amount::Union{<:SymbolicUtils.Sym,<:Int64})
    if term isa Average
        return term
    end
    if SymbolicUtils.istree(term)
        args = []
        op = SymbolicUtils.operation(term)
        for i = 1:length(arguments(term))
            push!(args,split_sums(arguments(term)[i],ind,amount))
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
            term_ = change_index(term_,ind,ind2)
            if isempty(term.metadata.nonEqualIndices)
                return (amount)*IndexedAverageSum(term_,ind2,Index[])
            end
            extrasum = IndexedAverageSum(term_,ind2,term.metadata.nonEqualIndices)
            return extrasum + (amount-1)*IndexedAverageSum(term_,ind2,Index[])
        end
    end
    if typeof(term) == SymbolicUtils.Sym{Parameter,IndexedAverageDoubleSum}
        Dsum = term.metadata
        if isequal(ind,Dsum.sumIndex)
            #create multiple doubleSums with the same innerSum
            ind2 = Index(Dsum.sumIndex.hilb,Dsum.sumIndex.name,(Dsum.sumIndex.rangeN/amount),Dsum.sumIndex.specHilb)
            if isempty(Dsum.nonEqualIndices)
                return amount*IndexedAverageDoubleSum(Dsum.innerSum,ind2,Index[])
            end
            extraSum = IndexedAverageDoubleSum(Dsum.innerSum,ind2,Dsum.nonEqualIndices)
            return extraSum + (amount-1)*Σ(Dsum.innerSum,ind2,Index[])
        elseif isequal(ind,Dsum.innerSum.metadata.sumIndex)
            #create doublesum of split innersums
            innerSums = split_sums(Dsum.innerSum,ind,amount)
            return IndexedAverageDoubleSum(innerSums,Dsum.sumIndex,Dsum.nonEqualIndices)
        end
    end
    return term
end
split_sums(x::Symbolics.Equation,ind::Index,amount) = Symbolics.Equation(x.lhs,split_sums(x.rhs,ind,amount))
function split_sums(me::AbstractMeanfieldEquations,ind::Index,amount)
    newEqs = Vector{Union{Missing,Symbolics.Equation}}(missing,length(me.equations))
    for i=1:length(me.equations)
        newEqs[i] = split_sums(me.equations[i],ind,amount)
    end
    newEqs = filter(x -> !isequal(x,missing), newEqs)
    vs = getLHS.(newEqs)
    varmap = make_varmap(vs, me.iv)
    return IndexedMeanfieldEquations(newEqs,me.operator_equations,vs,me.operators,me.hamiltonian,me.jumps,me.jumps_dagger,me.rates,me.iv,varmap,me.order)
end
split_sums(x,ind,amount) = x

scale(eqs::IndexedMeanfieldEquations;kwargs...) = subst_reds(scaleME(eqs;kwargs...);scaling=true)

"""
    value_map(ps::Vector,p0::Vector)

A Function to create parameter values for indexed Variables more convenient.

# Arguments
*`ps::Vector`: A vector of parameters, that have no value assigned to them.
*`p0::Vector`: A vector for numeric values, that should get assigned to the corresponding
    entry in the `ps` vector. For Single-Indexed Variables the entry in the vector can also be again
    a Vector, that has an amount of entries as the index of the variables has range. For Double-Indexed
    Variables, this can also be a Matrix of a dimension, that corresponds to the ranges of the indices
    of the given variable.

"""
function value_map(ps::Vector,p0::Vector;mapping=nothing)
    length(ps) != length(p0) && error("Vectors given have non-equal length!")

    if !=(mapping,nothing) && mapping isa Pair
        mapping_ = Dict{SymbolicUtils.Sym,Int64}(first(mapping)=>last(mapping))
        mapping = mapping_
    end
    if mapping === nothing
        mapping = Dict{SymbolicUtils.Sym,Int64}()
    end

    dict = Dict{Sym{Parameter, Base.ImmutableDict{DataType, Any}},ComplexF64}()
    for i=1:length(ps)
        dicVal = nothing
        if ps[i] isa SymbolicUtils.Sym{Parameter, IndexedVariable}
            if p0[i] isa Vector || p0[i] isa Number
                dicVal = create_value_map(ps[i],p0[i];mapping)
            else
                error("cannot resolve entry at $i-th position in values-vector")
            end
        elseif ps[i] isa SymbolicUtils.Sym{Parameter, DoubleIndexedVariable}
            if p0[i] isa Matrix || p0[i] isa Number
                dicVal = create_value_map(ps[i],p0[i];mapping)
            end
        else
            push!(dict,ps[i]=>p0[i])
            continue
        end
        dict = merge(dict,dicVal)
    end
    return collect(dict)
end

