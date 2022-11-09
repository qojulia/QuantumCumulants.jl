#File for adapting the meanfield, complete,... algorithms to indices. could be included in the already existing meanfield file.
#I did not want to change anything on the files in the already existing package, so I just copied and renamed a bunch of stuff here.


#function that takes indexed operators and double indexed varaibles to calculate the meanfield equations
#the jump operators have to have same indices as the indices specified by the double indexed variable
struct IndexedMeanfieldEquations <: AbstractMeanfieldEquations
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
"""
    indexed_meanfield(ops::Vector,H::QNumber,J::Vector;
        Jdagger::Vector=adjoint.(J),rates=ones(length(J)))

Compute the set of equations for the indexed-operators [`IndexedOperator`](@ref) in `ops` under the Hamiltonian
`H` and with loss operators contained in `J`. The resulting equation is
equivalent to the Quantum-Langevin equation where noise is neglected.
This is a modified version of the [`meanfield`](@ref) function, that can now also take [`IndexedOperator`](@ref)
entities for both the `ops` argument aswell as for the `J` arguments.
See also: [`meanfield`](@ref).

# Arguments
*`ops::Vector`: The operators of which the equations are to be computed.
*`H::QNumber`: The Hamiltonian describing the reversible dynamics of the
    system.
*`J::Vector{<:QNumber}`: A vector containing the collapse operators of
    the system. A term of the form
    ``\\sum_i J_i^\\dagger O J_i - \\frac{1}{2}\\left(J_i^\\dagger J_i O + OJ_i^\\dagger J_i\\right)``
    is added to the Heisenberg equation.

# Optional argumentes
*`Jdagger::Vector=adjoint.(J)`: Vector containing the hermitian conjugates of
    the collapse operators.
*`rates=ones(length(J))`: Decay rates corresponding to the collapse operators in `J`.
*`multithread=false`: Specify whether the derivation of equations for all operators in `ops`
    should be multithreaded using `Threads.@threads`.
*`simplify=true`: Specify whether the derived equations should be simplified.
*`order=nothing`: Specify to which `order` a [`cumulant_expansion`](@ref) is performed.
    If `nothing`, this step is skipped.
*`mix_choice=maximum`: If the provided `order` is a `Vector`, `mix_choice` determines
    which `order` to prefer on terms that act on multiple Hilbert spaces.
*`iv=SymbolicUtils.Sym{Real}(:t)`: The independent variable (time parameter) of the system.

"""
function indexed_meanfield(a::Vector,H,J;Jdagger::Vector=adjoint.(J),rates=ones(Int,length(J)),
    multithread=false,
    simplify::Bool=true,
    order=nothing,
    mix_choice=maximum,
    iv=SymbolicUtils.Sym{Real}(:t))

    for ind in get_indices(a)
        if ind in get_indices(H)
            error("Index $(ind.name) in operator-vector is already used in H!")
        end
    end

    # Derive operator equations
    rhs = Vector{Any}(undef, length(a))
    imH = im*H
    for i=1:length(a)
        try
            rhs_ = commutator(imH,a[i])
            rhs_diss = indexed_master_lindblad(a[i],J,Jdagger,rates)
            indices = get_indices(a[i])
            if length(indices) <= 1
                rhs[i] = rhs_ + rhs_diss
            else #everything on lhs commutes -> reorder corresponding terms on rhs
                mapping = Tuple{Index,Index}[]
                for j = 1:length(indices)
                    for k = 1:j
                        if k != j
                            push!(mapping,(indices[k],indices[j]))
                        end
                    end
                end
                rhs[i] = reorder((rhs_+rhs_diss),mapping)
            end    
        catch err
            println("could not calculate meanfield-equations for operator $(a[i])")
            rethrow(err)
        end
    end

    # Average
    vs = map(average, a)
    rhs_avg = map(average, rhs)
    if simplify
        rhs_avg = map(SymbolicUtils.simplify, rhs_avg)
    end
    rhs = map(undo_average, rhs_avg)

    if order !== nothing
        rhs_avg = [cumulant_expansion(r, order; simplify=simplify, mix_choice=mix_choice) for r∈rhs_avg]
    end

    eqs_avg = [Symbolics.Equation(l,r) for (l,r)=zip(vs,rhs_avg)]
    eqs = [Symbolics.Equation(l,r) for (l,r)=zip(a,rhs)]
    varmap = make_varmap(vs, iv)

    me = IndexedMeanfieldEquations(eqs_avg,eqs,vs,a,H,J,Jdagger,rates,iv,varmap,order)
    if has_cluster(H)
        return scale(me;simplify=simplify,order=order,mix_choice=mix_choice)
    else
        return me
    end
end
indexed_meanfield(a::QNumber,args...;kwargs...) = indexed_meanfield([a],args...;kwargs...)
indexed_meanfield(a::Vector,H;kwargs...) = indexed_meanfield(a,H,[];Jdagger=[],kwargs...)

function indexed_master_lindblad(a_,J,Jdagger,rates)
    args = Any[]
    for k=1:length(J)
        if typeof(J[k]) == IndexedOperator
            c1 = 0.5*rates[k]*Jdagger[k]*commutator(a_,J[k])
            c2 = 0.5*rates[k]*commutator(Jdagger[k],a_)*J[k]
            c = nothing
            args_ = []
            push_or_append_nz_args!(args_,c1)
            push_or_append_nz_args!(args_,c2)
            if isempty(args_)
                c = 0
            else
                c = sum(args_)
            end
            if !SymbolicUtils._iszero(c)
                push!(args, SingleSum(c,J[k].ind,Index[]))
            end
        else
            if typeof(rates[k]) == SymbolicUtils.Sym{Parameter,DoubleIndexedVariable}
                if J[k][1].ind != rates[k].metadata.ind1
                    error("unequal index of first jump operator and variable")
                end
                if J[k][2].ind != rates[k].metadata.ind2
                    error("unequal index of second jump operator and variable")
                end
                c1 = Jdagger[k][1]*commutator(a_,J[k][2])
                c2 = commutator(Jdagger[k][1],a_)*J[k][2]
                c = 0.5*rates[k]*(c1+c2)
                push!(args,IndexedDoubleSum(SingleSum((c),J[k][1].ind,Index[]),J[k][2].ind,Index[]))
            elseif isa(rates[k],SymbolicUtils.Symbolic) || isa(rates[k],Number) || isa(rates[k],Function)
                c1 = 0.5*rates[k]*Jdagger[k]*commutator(a_,J[k])
                c2 = 0.5*rates[k]*commutator(Jdagger[k],a_)*J[k]
                push_or_append_nz_args!(args, c1)
                push_or_append_nz_args!(args, c2)
            elseif isa(rates[k],Matrix)
                for i=1:length(J[k]), j=1:length(J[k])
                    c1 = 0.5*rates[k][i,j]*Jdagger[k][i]*commutator(a_,J[k][j])
                    c2 = 0.5*rates[k][i,j]*commutator(Jdagger[k][i],a_)*J[k][j]
                    push_or_append_nz_args!(args, c1)
                    push_or_append_nz_args!(args, c2)
                end
            else 
                error("Unknown rates type!")
            end
        end
    end
    isempty(args) && return 0
    return QAdd(args)
end

"""
    indexed_complete(de::MeanfieldEquations)

From a set of differential equation of averages, find all averages that are missing
and derive the corresponding equations of motion. Uses [`find_missing`](@ref)
and [`indexed_meanfield`](@ref) to do so.

Optional arguments
==================

*`order=de.order`: The order at which the [`cumulant_expansion`](@ref) is performed
    on the newly derived equations. If `nothing`, the order is inferred from the
    existing equations.
*`filter_func=nothing`: Custom function that specifies whether some averages should
    be ignored when completing a system. This works by calling `filter!(filter_func, missed)`
    where `missed` is the vector resulting from [`find_missing`](@ref). Occurrences
    of averages for which `filter_func` returns `false` are substituted to 0.
*`extra_indices`: A Vector of symbols, representing extra [`Index`](@ref) entities, which are
    needed and created in the process of finding missing terms.
*`kwargs...`: Further keyword arguments are passed on to [`indexed_meanfield`](@ref) and
    simplification.

see also: [`find_missing`](@ref), [`indexed_meanfield`](@ref), [`meanfield`](@ref), [`find_missing_sums`](@ref)
"""
function indexed_complete(de::AbstractMeanfieldEquations;kwargs...)
    de_ = deepcopy(de)
    indexed_complete!(de_;kwargs...)
    return de_
end

"""
    indexed_complete!(de::MeanfieldEquations)

In-place version of [`indexed_complete`](@ref)
"""
function indexed_complete!(de::AbstractMeanfieldEquations;
    order=de.order,
    multithread=false,
    filter_func=nothing,
    mix_choice=maximum,
    simplify=true,
    extra_indices::Vector=[:i,:j,:k,:l,:m,:n,:p,:q,:r,:s,:t],
    scaling::Bool=false,
    kwargs...)

    maxNumb = maximum(length.(get_indices.(de.operators)))

    if isempty(extra_indices)
        error("can not complete equations with empty extra_indices!")
    end

    for ind in extra_indices
        if typeof(ind) != typeof(extra_indices[1])
            error("Cannot use extra_indices of different types. Use either only Symbols or Index-Objects!")
        end
    end

    if containsMultiple(get_all_indices(de))
        if extra_indices[1] isa Symbol
            error("It is not possible to complete equations, containing indices, that act on different hilbertspaces using Symbols as
            extra_indices. For this case use specific Indices.")
        end
        #maybe write also a check that checks for the indices being correct/enough
    end

    allInds = get_all_indices(de)
    if extra_indices[1] isa Symbol
        filter!(x -> x ∉ getIndName.(allInds),extra_indices)
    else
        filter!(x -> x ∉ allInds,extra_indices)
    end

    if de.order > maxNumb && de.order - maxNumb > length(extra_indices)
        error("Too few extra_indices provided! Please make sure that for higher orders of cumulant expansion, 
            you also use the extra_indices argument to provide additional indices for calculation. The Number of
            extra_indices provided should be at least $(de.order - maxNumb).
        ")
    end

    vs = de.states
    order_lhs = maximum(get_order.(vs))
    order_rhs = 0
    for i=1:length(de.equations)
        k = get_order(de.equations[i].rhs)
        k > order_rhs && (order_rhs = k)
    end
    if order === nothing
        order_ = max(order_lhs, order_rhs)
    else
        order_ = order
    end
    maximum(order_) >= order_lhs || error("Cannot form cumulant expansion of derivative; you may want to use a higher order!")

    if order_ != de.order
        for i=1:length(de.equations)
            lhs = de.equations[i].lhs
            rhs = cumulant_expansion(de.equations[i].rhs,order_;
                        mix_choice=mix_choice,
                        simplify=simplify)
            de.equations[i] = Symbolics.Equation(lhs, rhs)
        end
    end

    indices_lhs = get_all_indices(vs) #indices that are used in the beginning -> have priority

    #if there are no indices on the LHS -> use the first index of extra_indices as starting point
    #this first index is (if it is a symbol) created similar to the first occuring index in the Equations
    if isempty(indices_lhs) && extra_indices[1] isa Symbol
        for ind in allInds
            indices_lhs = [Index(ind.hilb,extra_indices[1],ind.range,ind.specHilb)]
            deleteat!(extra_indices,1)
            break
        end
    elseif isempty(indices_lhs) && extra_indices[1] isa Index
        indices_lhs = [extra_indices[1]]
        deleteat!(extra_indices,1)
    end

    #create (or append to the extras vector) all other extra_indices using the indices on the LHS
    extras = indices_lhs
    if !isempty(extra_indices)
        if extra_indices[1] isa Symbol
            first = extras[1]
            for name in extra_indices
                push!(extras,Index(first.hilb,name,first.range,first.specHilb))
                if length(extras) >= order_
                    break
                end
            end
        end
        if extra_indices[1] isa Index
            extras = [extras;extra_indices]
        end
    end
    if length(extras) < order_
        error("More extra_indices are needed for order $(order_)")
    end

    #at this point extras is a list of extra_indices, sorted by their priority 
    # (meaning that indices that were used in the ops of in indexed_meanfield come first)


    vhash = map(hash, vs)
    vs′ = map(_conj, vs)
    vs′hash = map(hash, vs′)
    filter!(!in(vhash), vs′hash)
    missed = find_missing(de.equations, vhash, vs′hash; get_adjoints=false)

    missed = find_missing_sums(missed,de;extra_indices=extras,scaling=scaling,kwargs...)
    missed = findMissingSpecialTerms(missed,de;scaling=scaling,kwargs...)
    
    missed = sortByIndex.(missed)
    filter!(x -> (x isa Average),missed)

    isnothing(filter_func) || filter!(filter_func, missed) # User-defined filter
    
    filter!(x -> filterComplete(x,de.states,scaling;kwargs...), missed)

    for i = 1:length(missed)
        minds = get_indices(missed[i])
        newMinds = copy(minds)
        for ind1 in minds
            extras_=filterExtras(ind1,extras)
            for k = 1:length(extras_)
                if findall(x->isequal(x,ind1),extras_)[1] > k && extras_[k] ∉ newMinds #this might go somewhat easier, maybe delete ind2 out of extras after each replacement somehow
                    missed[i] = change_index(missed[i],ind1,extras_[k])
                    newMinds = get_indices(missed[i])
                    break
                elseif findall(x->isequal(x,ind1),extras_)[1] <= k
                    break
                end
            end
        end
    end
    missed = unique(missed) #no duplicates

    while !isempty(missed)
        ops_ = [SymbolicUtils.arguments(m)[1] for m in missed]
        me = indexed_meanfield(ops_,de.hamiltonian,de.jumps;
            Jdagger=de.jumps_dagger,
            rates=de.rates,
            simplify=simplify,
            multithread=multithread,
            order=order_,
            mix_choice=mix_choice,
            iv=de.iv,
            kwargs...)
        _append!(de, me)
        vhash_ = hash.(me.states)
        vs′hash_ = hash.(_conj.(me.states))
        append!(vhash, vhash_)
        for i=1:length(vhash_)
            vs′hash_[i] ∈ vhash_ || push!(vs′hash, vs′hash_[i])
        end
        missed = find_missing(me.equations, vhash, vs′hash; get_adjoints=false)
        missed = find_missing_sums(missed,de;extra_indices=extras,scaling=scaling,kwargs...)
        missed = findMissingSpecialTerms(missed,de;scaling=scaling,kwargs...)
        missed = sortByIndex.(missed)

        filter!(x -> (x isa Average),missed)
        isnothing(filter_func) || filter!(filter_func, missed) # User-defined filter
    
        filter!(x -> filterComplete(x,de.states,scaling;kwargs...), missed)
        for i = 1:length(missed)
            minds = get_indices(missed[i])
            newMinds = copy(minds)
            for ind1 in minds
                extras_=filterExtras(ind1,extras)
                for k = 1:length(extras_)
                    if ind1 ∉ extras_
                        missed[i] = change_index(missed[i],ind1,extras_[1])
                        newMinds = get_indices(missed[i])
                        break
                    end
                    if findall(x->isequal(x,ind1),extras_)[1] > k && extras_[k] ∉ newMinds #this might go somewhat easier, maybe delete ind2 out of extras after each replacement somehow
                        missed[i] = change_index(missed[i],ind1,extras_[k])
                        newMinds = get_indices(missed[i])
                        break
                    elseif findall(x->isequal(x,ind1),extras_)[1] <= k
                        break
                    end
                end
            end
        end
        missed = sortByIndex.(missed)
        filter!(x -> filterComplete(x,de.states,scaling;kwargs...), missed)
        missed = unique(missed)
        missed = elimRed(missed;scaling=scaling,kwargs...)
    end

    if !isnothing(filter_func)
        # Find missing values that are filtered by the custom filter function,
        # but still occur on the RHS; set those to 0
        missed = find_missing(de.equations, vhash, vs′hash; get_adjoints=false)
        
        missed = find_missing_sums(missed,de;extra_indices=extras,checking=false,kwargs...) #checkin dissabled, since there might be some terms, that are redundant, but also zero -> set them to zero aswell forsa fety
        missed = findMissingSpecialTerms(missed,de;kwargs...)
    
        missed = sortByIndex.(missed) #this one might not be right (?)
        filter!(x -> (x isa Average),missed)

        filter!(!filter_func, missed)
        missed_adj = map(_adjoint, missed)
        subs = Dict(vcat(missed, missed_adj) .=> 0)
        for i=1:length(de.equations)
            de.equations[i] = (SymbolicUtils.substitute(de.equations[i], subs))
            de.states[i] = de.equations[i].lhs
        end
    end

    return de
end
filterComplete(x,states,scaling;kwargs...) = (isNotIn(getOps(x;scaling=scaling,kwargs...),getOps.(states;scaling=scaling,kwargs...),scaling) && isNotIn(getOps(sortByIndex(_conj(x));scaling=scaling,kwargs...),getOps.(states;scaling=scaling,kwargs...),scaling))

#gets all indices that are used in the states of the meanfieldequations
function get_all_indices(vec::Vector)
    inds = []
    for i = 1:length(vec)
        indices_ = get_indices(vec[i])
        isempty(indices_) && continue
        for ind in indices_
            if ind ∉ inds
                push!(inds,ind)
            end
        end
    end
    return inds
end
function get_all_indices(me::AbstractMeanfieldEquations)
    eqs = me.equations
    inds = []
    for eq in eqs
        for ind in get_indices(eq.rhs)
            if ind ∉ inds
                push!(inds,ind)
            end
        end
    end
    for ind in get_all_indices(me.states)
        if ind ∉ inds
            push!(inds,ind)
        end
    end
    for ind in get_indices(me.hamiltonian)
        if ind ∉ inds 
            push!(inds,ind)
        end
    end
    for ind in get_indices(me.jumps)
        if ind ∉ inds
            push!(inds,ind)
        end
    end
    return inds
end
function get_indices_equations(me::AbstractMeanfieldEquations)
    eqs = me.equations
    inds = []
    for ind in get_all_indices(me.states)
        if ind ∉ inds
            push!(inds,ind)
        end
    end
    for eq in eqs
        for ind in get_indices(eq.rhs)
            if ind ∉ inds
                push!(inds,ind)
            end
        end
    end
    return inds
end
function containsMultiple(inds::Vector) #checks if in a list of indices, they act on different sub-hilbertSpaces
    isempty(inds) && return false
    length(inds) == 1 && return false
    for i=2:length(inds)
        if inds[1].specHilb != inds[i].specHilb
            return true
        end
    end
    return false
end

"""
    find_missing_sums(missed,de::MeanfieldEquations)

From a initial set of differential equation of averages, find all averages that are missing
and inside a Symbolic sum. If a missing average contains one of the summation indices used in
the equations, the [`Index`](@ref) will be exchanged according to the keyword argument
`extra_indices`. Uses [`find_missing`](@ref).

# Arguments
*`missed`: A initial Vector of averages, representing the missed averages before calling
    this method.
*`de`: The set of equations, in which the missing averages are searched in.

# Optional arguments
*`extra_indices`: A Vector of symbols, representing extra [`Index`](@ref) entities, which are
    needed and created in the process of finding missing terms. This argument is required, if the order
    of the Meanfield-Equations exceeds 1 and the number of given symbols must match the corresponding order.
*`checking`: A Bool defining if the algorithm checks for adjoint values and duplicates, before adding a found
    average into the `missed` vector.
*`scaling`: A Bool defining the way how averages are added to the `missed` vector. If true only averages, whose
    operators (without indices) are not already inside the `missed` vector will be added.

see also: [`find_missing`](@ref), [`indexed_meanfield`](@ref), [`meanfield`](@ref), [`find_missing_sums`](@ref)
"""
function find_missing_sums(missed,de::AbstractMeanfieldEquations;extra_indices::Vector=[],checking=true,scaling=false,kwargs...)
    missed_ = copy(missed)
    extras = copy(extra_indices)
    sums = Any[]
    for eq in de.equations
        sums = checkIfSum(eq.rhs)   #get all sums of the rhs
        for sum in sums
            extras_filtered = filterExtras(sum.metadata.sum_index,extras) #all the extraIndices, that have the same specific hilbertspace as the sum
            checkAndChange(missed_,sum,extras_filtered,de.states,checking,scaling;kwargs...)
        end
    end   
    return missed_
end
function checkAndChange(missed_,sum,extras,states,checking,scaling;kwargs...)
    avrgs = getAvrgs(sum) #get vector of all avrgs in the sum
    for avr in avrgs
        changed_ = nothing
        avrg_inds = get_indices(avr)
        for i = 1:length(extras)    # check if the extraindex is already in use for the term -> if so use the next in line
            if extras[i] ∉ avrg_inds #this is wrong for multiple indices!!!!
                changed_ = change_index(avr,sum.metadata.sum_index,extras[i])
                break
            end
        end
        if isequal(changed_,nothing) #if avrg consists of none indexed operators, or only of operators that already have the desired indices -> push along the whole average
            changed_ = avr
        end
        
        if changed_ isa Term{AvgSym, Nothing}
            changed = sortByIndex(changed_)    # this can be done, since terms inside the sum commute anyway
            if !(changed isa Average)
                continue
            end
            if checking     # checks if the missed term found by this function is already in the list of missed averages (checking = false is needed for substituting and filtering of adjoint summation terms since the same function is called when substituting inside sums)
                if filterComplete(changed,states,scaling;kwargs...) && filterComplete(changed,missed_,scaling;kwargs...)
                    push!(missed_,changed)
                end
            else
                push!(missed_,changed)
            end
        end
    end
end
# function that filters indices so that only indices get replaces that "live" in the same sub-hilbertspace
function filterExtras(wanted,extras)
    filtered = []
    for ind in extras
        if isequal(ind.specHilb,wanted.specHilb)
            push!(filtered,ind)
        end
    end
    return filtered
end

#This is basically just the find_missing function but for the extra defined SpecialTerms, which have extra index constraints and for that also some extra checks
function findMissingSpecialTerms(missed,me::AbstractMeanfieldEquations;scaling=false,kwargs...)
    missed_ = copy(missed)
    vs = me.states
    vhash = map(hash, vs)
    vs′=map(_conj, vs)
    vs′hash = map(hash, vs′)
    filter!(!in(vhash), vs′hash)
    missed_hashes = map(hash,missed_)
    for eq in me.equations
        if eq.rhs isa Number
            continue
        else
            for arg in arguments(eq.rhs)
                if typeof(arg) == Int64
                    continue
                end
                if typeof(arg) == SymbolicUtils.Sym{Parameter,SpecialIndexedAverage}
                    missed_ = find_missing!(missed_, missed_hashes, arg.metadata.term, vhash, vs′hash)
                end
                if typeof(arg) <: SymbolicUtils.Mul
                    for arg_ in arguments(arg)
                        if typeof(arg_) == Int64
                            continue
                        end
                        if typeof(arg_) == SymbolicUtils.Sym{Parameter,SpecialIndexedAverage}
                            missed_ = find_missing!(missed_, missed_hashes, arg_.metadata.term, vhash, vs′hash)
                        end
                    end
                end
            end
        end
        
    end
    return missed_
end

#Utility Functions
function sortByIndex(term::Term{AvgSym, Nothing})
    arg = arguments(term)[1]
    if typeof(arg) <: QMul
        args_nc = arg.args_nc
        sort!(args_nc,by=getIndName)
        if isempty(args_nc)
            return 0
        elseif length(args_nc) == 1
            return average(args_nc[1])
        end
        mult = *(args_nc...)
        return average(mult)
    end
    return term
end
#checks if there is a sum (or multiple) in the equation rhs term, if so it returns the sums (as symbol) as a vector
function checkIfSum(term)
    sums = Any[]
    if typeof(term) == SymbolicUtils.Sym{Parameter,IndexedAverageSum}
        return [term]
    elseif typeof(term) <: SymbolicUtils.Add || typeof(term) <: SymbolicUtils.Mul
        args = arguments(term)
        for arg in args
            sums = vcat(sums,checkIfSum(arg))
        end
    end
    return sums
end
function checkIfDSum(term)
    sums = Any[]
    if typeof(term) == SymbolicUtils.Sym{Parameter,IndexedAverageDoubleSum}
        return [term]
    elseif typeof(term) <: SymbolicUtils.Add || typeof(term) <: SymbolicUtils.Mul
        args = arguments(term)
        for arg in args
            sums = vcat(sums,checkIfDSum(arg))
        end
    end
    return sums
end
function getOps(term::SymbolicUtils.Mul;scaling::Bool=false,kwargs...)
    for arg in arguments(term)
        if arg isa Average
            return getOps(arg;scaling=scaling,kwargs...)
        end
    end
    return []
end
function getOps(ops::Term{AvgSym, Nothing};scaling::Bool=false,h=nothing,kwargs...)
    args = arguments(ops)[1]
    if typeof(args) <: QMul 
        return getOps(args;scaling=scaling,kwargs...)
    elseif scaling && (args isa NumberedOperator || args isa IndexedOperator)
        return Any[args.op]
    else #single op in term
        return Any[args]
    end
end
function getOps(ops::QMul;scaling::Bool=false,h=nothing,kwargs...)
    arr = Any[]
    for arg in ops.args_nc
        if scaling && arg isa NumberedOperator
            push!(arr,arg.op)
        elseif scaling && arg isa IndexedOperator
            if !=(h,nothing)
                if isequal(arg.ind.specHilb, h)
                    push!(arr,arg.op)
                else
                    push!(arr,arg)
                end
            else
                push!(arr,arg.op)
            end
        else
            push!(arr,arg)
        end
    end
    return arr
end
getOps(x;scaling::Bool=false,kwargs...) = Vector{Vector{Any}}()
function isNotIn(terms::Vector{Any},vects::Vector{Vector{Any}},scaling::Bool)
    for vect in vects
        if scaling
            if hasSameOps(vect,terms)
                return false
            end
        else
            if isequal(vect,terms)
                return false
            end
        end
    end
    return true
end
isNotIn(terms::Vector{Any},vects::Vector{Any},scaling) = isNotIn(terms,[vects],scaling)

function hasSameOps(vec1::Vector,vec2::Vector)
    length(vec1) != length(vec2) && return false
    for elem in vec1
        if elem ∉ vec2
            return false
        end
    end
    for elem in vec2
        if elem ∉ vec1
            return false
        end
    end
    return true
end

# this is most likely not the best method, works fine however
function elimRed(missed::Vector;scaling::Bool=false,kwargs...)
    all = getOps.(missed)
    for i=1:length(missed)
        if i > length(missed)
            break
        end
        for j=1:length(missed)
            if j > length(missed)
                break
            end
            if hasSameOps(getOps(_conj(missed[i]);scaling=scaling,kwargs...),getOps(missed[j];scaling=scaling,kwargs...)) && !hasSameOps(getOps(missed[i];scaling=scaling,kwargs...),getOps(missed[j];scaling=scaling,kwargs...))
                deleteat!(missed,j)
            end
        end
    end
    return missed
end

"""
    subst_reds(de::MeanfieldEquations)

Function that substitutes possible redundant conjugate averages inside the given Equations with their corresponding
average given as the conjugate of one of the left-hand-side (of the equations) averages.

# Optional Arguments
*`scaling`: A Bool defining the way how averages are added to the `missed` vector. If true only averages, whose
    operators (without indices) are not already inside the `missed` vector will be added.
"""
function subst_reds(me::AbstractMeanfieldEquations;scaling::Bool=false,kwargs...)
    eqs = Vector{Union{Missing,Symbolics.Equation}}(missing,length(me.equations))
    to_sub = find_missing(me)
    filter!(x->!(x in me.states)&&!(_conj(x) in me.states),to_sub)
    if scaling
        #brute force for now
        #this should not take too long
        #proportional to number of elems in to_sub
        states = me.states
        op_states = getOps.(states;scaling=scaling)
        to_insert = Vector{Any}(nothing,length(to_sub))
        counter = 1
        for elem in to_sub
            ind_ = findfirst(x->hasSameOps(x,getOps(elem;scaling=scaling)),op_states)
            if !=(ind_,nothing)
                to_insert[counter] = states[ind_]
                counter = counter + 1
                continue
            end
            ind_ = findfirst(x->hasSameOps(x,getOps(_conj(elem);scaling=scaling)),op_states)
            if !=(ind_,nothing)
                to_insert[counter] = _conj(states[ind_])
                counter = counter + 1
                continue
            end
        end
    else  
        to_insert = _conj.(to_sub)
    end
    subs = Dict(to_sub .=> to_insert)
    for i=1:length(me.equations)
        eqs[i] = SymbolicUtils.substitute(me.equations[i],subs)
    end
    return IndexedMeanfieldEquations(eqs,me.operator_equations,me.states,me.operators,me.hamiltonian,me.jumps,me.jumps_dagger,me.rates,me.iv,me.varmap,me.order)
end
#=
function subst_reds(term::Average,vect;kwargs...)
    term in vect && return term
    _conj(term) in vect && return conj(_conj(term))
    @warn "unknown term encountered $(term)"
    return term
end
function subst_reds(term,vect;kwargs...)
    if istree(term)
        op = operation(term)
        args = map(x->subst_reds(x,vect;kwargs...),arguments(term))
        return SymbolicUtils.similarterm(term, op, args)
    else
        return term
    end
end
=#
#=
function subst_reds(term,vect::Vector;scaling::Bool=false,kwargs...)
    if term isa Average
        for i = 1:length(vect)
            temp = getOps(_conj(vect[i]);scaling=scaling,kwargs...)
            if hasSameOps(getOps(term;scaling=scaling,kwargs...),temp) && !(hasSameOps(getOps(term;scaling=scaling,kwargs...),getOps(vect[i];scaling=scaling,kwargs...)))
                return conj(vect[i])
            end
            if hasSameOps(getOps(term;scaling=scaling,kwargs...),getOps(vect[i];scaling=scaling,kwargs...))
                return vect[i]
            end
        end
    elseif istree(term)
        op = operation(term)
        args = map(x->subst_reds(x,vect;scaling=scaling,kwargs...),arguments(term))
        return SymbolicUtils.similarterm(term, op, args)
    elseif typeof(term) == SymbolicUtils.Sym{Parameter,SpecialIndexedAverage}
        return SpecialIndexedAverage(subst_reds(term.metadata.term,vect;scaling=scaling,kwargs...),term.metadata.indexMapping)
    end
    return term
end
=#

complete(eqs::IndexedMeanfieldEquations;kwargs...) = indexed_complete(eqs;kwargs...)
complete!(eqs::IndexedMeanfieldEquations;kwargs...) = indexed_complete!(eqs;kwargs...) 

"""
    evaluate(eqs::IndexedMeanfieldEquations)

Function, that evaluates a given [`MeanfieldEquations`](@ref) entity and returns again equations,
where indices have been inserted and sums evaluated.

# Arguments
*`me::MeanfieldEquations`: A [`MeanfieldEquations`](@ref) entity, which shall be evaluated.

# Optional argumentes
*`mapping::Dict{SymbolicUtils.Sym,Int64}=Dict{Symbol,Int64}()`: A seperate dictionary, to
    specify any symbolic limits used when [`Index`](@ref) entities were defined. This needs
    to be specified, when the equations contain summations, for which the upper bound is given
    by a Symbolic.

see also: [`evalME`](@ref)
"""
function evaluate(eqs::IndexedMeanfieldEquations;mapping=nothing,kwargs...)
    if !=(mapping,nothing) && mapping isa Pair
        mapping_ = Dict{SymbolicUtils.Sym,Int64}(first(mapping) => last(mapping));
        mapping = mapping_
    end
    if mapping === nothing
        mapping =  Dict{SymbolicUtils.Sym,Int64}();
    end
    return subst_reds(evalME(eqs;mapping=mapping,kwargs...))
end