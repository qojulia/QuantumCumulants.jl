 #Main file for manipulating indexed averages and sums over averages.

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
            ranges = SQA.get_range.(inds)
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
        if meta.sum_index.range isa BasicSymbolic{<:Complex{<:Real}}
        # if meta.sum_index.range isa BasicSymbolic{<:CNumber}
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
