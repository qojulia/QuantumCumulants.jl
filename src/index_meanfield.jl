#File for adapting the meanfield, complete,... algorithms to indices. could be included in the already existing meanfield file.
#I did not want to change anything on the files in the already existing package, so I just copied and renamed a bunch of stuff here.


#function that takes indexed operators and double indexed varaibles to calculate the meanfield equations
#the jump operators have to have same indices as the indices specified by the double indexed variable
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
*`iv=ModelingToolkit.t`: The independent variable (time parameter) of the system.

"""
function indexed_meanfield(a::Vector,H,J;Jdagger::Vector=adjoint.(J),rates=ones(Int,length(J)),
    multithread=false,
    simplify::Bool=true,
    order=nothing,
    mix_choice=maximum,
    iv=MTK.t_nounits,
    kwargs...)

    ind_J = []
    for j in J
        ind = get_indices(j)
        if ind ∉ ind_J
            push!(ind_J,ind)
        end
    end

    for ind in get_indices(a)
        if ind in get_indices(H)
            error("Index $(ind.name) in operator-vector is already used in H!")
        end
        if ind in ind_J
            error("Index $(ind.name) in operator-vector is already used in J!")
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
    return me
end
indexed_meanfield(a::QNumber,args...;kwargs...) = indexed_meanfield([a],args...;kwargs...)
indexed_meanfield(a::Vector,H;kwargs...) = indexed_meanfield(a,H,[];Jdagger=[],kwargs...)

function indexed_master_lindblad(a_,J,Jdagger,rates)
    args = Any[]
    for k=1:length(J)
        if !isempty(get_indices(J[k])) && !(rates[k] isa BasicSymbolic{DoubleIndexedVariable})
            if J[k] isa SQA.Sums
                c1 = 0.5*rates[k]*Jdagger[k]*commutator(a_,J[k])
                c2 = 0.5*rates[k]*commutator(Jdagger[k],a_)*J[k]
                c = c1 + c2
                push!(args, c)
            else #J[k] indexed variable
                inds = get_indices(J[k])
                length(inds) != 1 && error("Jump operators with multiple different indices are not supported!")
                c1 = 0.5*rates[k]*Jdagger[k]*commutator(a_,J[k])
                c2 = 0.5*rates[k]*commutator(Jdagger[k],a_)*J[k]
                c = c1 + c2
                !SymbolicUtils._iszero(c) && push!(args, SingleSum(c,inds[1],Index[]))
            end
        else
            if rates[k] isa BasicSymbolic{DoubleIndexedVariable}
                meta = SymbolicUtils.metadata(rates[k])[DoubleIndexedVariable]
                if !(J[k] isa Vector) && !(isempty(get_indices(J[k])))
                    inds = get_indices(J[k])
                    length(inds) != 1 && error("Jump operators with multiple different indices are not supported!")
                    ind_ = inds[1]
                    if !isequal(ind_,meta.ind1) && !isequal(ind_,meta.ind2)
                        error("unequal index of jump operator and variable!")
                    elseif isequal(ind_,meta.ind1) && !isequal(ind_,meta.ind2)
                        vec = [J[k],change_index(J[k],ind_,meta.ind2)]
                    elseif !isequal(ind_,meta.ind1) && isequal(ind_,meta.ind2)
                        vec = [change_index(J[k],ind_,meta.ind1),J[k]]
                    end
                    vecdagger = adjoint.(vec)
                    c1 = vecdagger[1]*commutator(a_,vec[2])
                    c2 = commutator(vecdagger[1],a_)*vec[2]
                    c = 0.5*rates[k]*(c1+c2)
                    push!(args,∑(c,vec[1].ind,vec[2].ind))
                else
                    length(get_indices(J[k][1])) != 1 && error("Jump operators with multiple different indices are not supported!")
                    length(get_indices(J[k][2])) != 1 && error("Jump operators with multiple different indices are not supported!")
                    ind_1 = get_indices(J[k][1])[1]
                    ind_2 = get_indices(J[k][2])[1]
                    if !isequal(ind_1, meta.ind1)
                        error("unequal index of first jump operator and variable")
                    end
                    if !isequal(ind_2, meta.ind2)
                        error("unequal index of second jump operator and variable")
                    end
                    c1 = Jdagger[k][1]*commutator(a_,J[k][2])
                    c2 = commutator(Jdagger[k][1],a_)*J[k][2]
                    c = 0.5*rates[k]*(c1+c2)
                    push!(args,∑(c,ind_1,ind_2))
                end
            elseif isa(rates[k],SymbolicUtils.Symbolic) || isa(rates[k],Number) || isa(rates[k],Function)
                c1 = 0.5*rates[k]*Jdagger[k]*commutator(a_,J[k])
                c2 = 0.5*rates[k]*commutator(Jdagger[k],a_)*J[k]
                SQA.push_or_append_nz_args!(args, c1)
                SQA.push_or_append_nz_args!(args, c2)
            elseif isa(rates[k],Matrix)
                for i=1:length(J[k]), j=1:length(J[k])
                    c1 = 0.5*rates[k][i,j]*Jdagger[k][i]*commutator(a_,J[k][j])
                    c2 = 0.5*rates[k][i,j]*commutator(Jdagger[k][i],a_)*J[k][j]
                    SQA.push_or_append_nz_args!(args, c1)
                    SQA.push_or_append_nz_args!(args, c2)
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
    scaling=false,
    kwargs...)

    ops_indices = get_indices(de.operators)
    extra_indices = copy(extra_indices) # otherwise the global variable will be overwritten
    isempty(extra_indices) && error("can not complete equations with empty extra_indices!")

    for ind in extra_indices
        if typeof(ind) != typeof(extra_indices[1])
            error("Cannot use extra_indices of different types. Use either only Symbols or Index-Objects!")
        end
    end

    allInds = get_all_indices(de)

    if extra_indices[1] isa Symbol
        filter!(x -> x ∉ SQA.getIndName.(allInds),extra_indices)
    else
        filter!(x -> x ∉ allInds,extra_indices)
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

    if order isa Vector
        order_max = maximum(order)
    else
        if order === nothing
            order_max = order_
        else
            order_max = order
        end
    end

    num_inds =  length(unique([ops_indices;extra_indices]))
    num_ind_hilbert = length(unique([ind.aon for ind in get_indices([de.hamiltonian,de.jumps, de.jumps_dagger])]))

    if order_max*num_ind_hilbert > num_inds
        error("Too few extra_indices provided! Please make sure that for higher orders of cumulant expansion,
            you also use the extra_indices argument to provide additional indices for calculation. The Number of
            extra_indices provided should be at least $(order_max - num_inds).")
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

    if containsMultiple(allInds) && extra_indices[1] isa Symbol
        dic = split_inds(allInds)
        filter!(x->!isempty(last(x)),dic)
        dic2 = Dict{Int,Any}(i => Any[] for i in keys(dic))
        for k in keys(dic)
            dic2[k] = filter(x->isequal(x.aon,k),indices_lhs)
            ind = allInds[findfirst(x->isequal(x.aon,k),allInds)]
            while length(dic2[k]) < order_max
                push!(dic2[k],Index(ind.hilb,extra_indices[1],ind.range,k))
                deleteat!(extra_indices,1)
            end
        end
        extras = []
        for vec in values(dic2)
            extras = [extras;vec]
        end
    else
        #if there are no indices on the LHS -> use the first index of extra_indices as starting point
        #this first index is (if it is a symbol) created similar to the first occuring index in the Equations
        if isempty(indices_lhs) && extra_indices[1] isa Symbol
            for ind in allInds
                indices_lhs = [Index(ind.hilb,extra_indices[1],ind.range,ind.aon)]
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
                    push!(extras,Index(first.hilb,name,first.range,first.aon))
                    if length(extras) >= order_max
                        break
                    end
                end
            end
            if extra_indices[1] isa Index
                extras = [extras;extra_indices]
            end
        end
        if length(extras) < order_max
            error("More extra_indices are needed for order $(order_max)")
        end

    end

    #at this point extras is a list of extra_indices, sorted by their priority
    # (meaning that indices that were used in the ops of in indexed_meanfield come first)

    vhash = map(hash, vs)
    vs′ = map(_inconj, vs)
    vs′hash = map(hash, vs′)
    filter!(!in(vhash), vs′hash)
    missed = find_missing(de.equations, vhash, vs′hash; get_adjoints=false)
    missed = find_missing_sums(missed,de;extra_indices=extras,kwargs...)
    missed = find_missing_Dsums(missed,de;extra_indices=extras,kwargs...)

    missed = inorder!.(missed)
    filter!(x -> (x isa Average),missed)
    isnothing(filter_func) || filter!(filter_func, missed) # User-defined filter

    filter!(x -> filterComplete(x,de.states,scaling;kwargs...), missed) # filterComplete does for whatever reason interfere with the order in the averages...
    missed = inorder!.(missed)   # ...thats why we resort here the missed again

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

    missed = elimRed!(missed;kwargs...)
    missed = inorder!.(missed)

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
        vs′hash_ = hash.(_inconj.(me.states))
        append!(vhash, vhash_)
        for i=1:length(vhash_)
            vs′hash_[i] ∈ vhash_ || push!(vs′hash, vs′hash_[i])
        end
        missed = find_missing(me.equations, vhash, vs′hash; get_adjoints=false)
        missed = find_missing_sums(missed,de;extra_indices=extras,kwargs...)
        missed = find_missing_Dsums(missed,de;extra_indices=extras,kwargs...)

        missed = inorder!.(missed)

        filter!(x -> (x isa Average),missed)
        isnothing(filter_func) || filter!(filter_func, missed) # User-defined filter

        filter!(x -> filterComplete(x,de.states,scaling;kwargs...), missed)
        missed = inorder!.(missed)

        for i = 1:length(missed)
            minds = get_indices(missed[i])
            newMinds = copy(minds)
            for ind1 in minds
                extras_=filterExtras(ind1,extras)
                for k = 1:length(extras_)
                    ind1_inds = findall(x->isequal(x,ind1),extras_)
                    if (isempty(ind1_inds) || findall(x->isequal(x,ind1),extras_)[1] > k) && extras_[k] ∉ newMinds #this might go somewhat easier, maybe delete ind2 out of extras after each replacement somehow
                        missed[i] = change_index(missed[i],ind1,extras_[k])
                        newMinds = get_indices(missed[i])
                        break
                    elseif !isempty(ind1_inds) && findall(x->isequal(x,ind1),extras_)[1] <= k
                        break
                    end
                end
            end
        end

        missed = inorder!.(missed)
        missed = unique(missed)
        missed = filter!(x -> filterComplete(x,de.states,scaling;kwargs...), missed)
        missed = inorder!.(missed)
        missed = elimRed!(missed;scaling=scaling,kwargs...)

    end

    if !isnothing(filter_func)
        # Find missing values that are filtered by the custom filter function,
        # but still occur on the RHS; set those to 0
        missed = find_missing(de.equations, vhash, vs′hash; get_adjoints=false)
        missed = find_missing_sums(missed,de;extra_indices=extras,checking=false,kwargs...) #checkin dissabled, since there might be some terms, that are redundant, but also zero -> set them to zero aswell forsa fety
        missed = find_missing_Dsums(missed,de;extra_indices=extras,checking=false,kwargs...)
        missed = find_missing_sums(missed,de;checking=false,kwargs...)
        missed = find_missing_Dsums(missed,de;checking=false,kwargs...)

        missed = inorder!.(missed) #this one might not be right (?) -> not even needed i think
        filter!(x -> (x isa Average),missed)

        filter!(!filter_func, missed)
        missed_adj = map(_inconj, missed)
        subs = Dict(vcat(missed, missed_adj) .=> 0)
        for i=1:length(de.equations)
            de.equations[i] = (SymbolicUtils.substitute(de.equations[i], subs))
            de.states[i] = de.equations[i].lhs
        end
    end

    return de
end

#function that splits a vector of indices into a vector of vectors, in which all indices are acting on the same hilbertspace
function split_inds(inds::Vector)
    h = inds[1].hilb
    spaces = h.spaces
    dic = Dict{Int,Any}(i => Any[] for i = 1:length(spaces))
    for ind in inds
        push!(dic[ind.aon],ind)
    end
    return dic
end

filterComplete(x,states,scaling;kwargs...) = isNotIn(x,states,scaling;kwargs...) && isNotIn(_inconj(x),states,scaling;kwargs...)

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
    for ind in get_all_indices(me.rates)
        if ind ∉ inds
            push!(inds,ind)
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
        if inds[1].aon != inds[i].aon
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
    eqs = de.equations
    for eq in eqs
        sums = checkIfSum(eq.rhs)   #get all sums of the rhs
        for sum in sums
            meta = SymbolicUtils.metadata(sum)[IndexedAverageSum]
            extras_filtered = filterExtras(meta.sum_index,extras) #all the extraIndices, that have the same specific hilbertspace as the sum
            checkAndChange(missed_,sum,extras_filtered,de.states,checking,scaling;kwargs...)
        end
    end
    return inorder!.(missed_)
end
#this function might not be so fast
function find_missing_Dsums(missed,de::AbstractMeanfieldEquations;extra_indices::Vector=[],checking=true,scaling=false,kwargs...)
    missed_ = copy(missed)
    extras = copy(extra_indices)
    for eq in de.equations
        Dsums = getDSums(eq.rhs)
        for sum in Dsums
            checkAndChangeDsum(missed_,sum,extras,de.states,checking,scaling;kwargs...)
        end
    end
    return inorder!.(missed_)
end
function checkAndChange(missed_,sum,extras,states,checking,scaling;kwargs...)
    meta = SymbolicUtils.metadata(sum)[IndexedAverageSum]
    avrgs = getAvrgs(sum) #get vector of all avrgs in the sum
    if avrgs isa Average
        avrgs = [avrgs]
    end
    for avr in avrgs
        changed_ = nothing
        avrg_inds = get_indices(avr)
        for i = 1:length(extras)    # check if the extraindex is already in use for the term -> if so use the next in line
            if extras[i] ∉ avrg_inds #this is wrong for multiple indices!!!!
                changed_ = change_index(avr,meta.sum_index,extras[i])
                break
            end
        end
        if changed_ === nothing #if avrg consists of none indexed operators, or only of operators that already have the desired indices -> push along the whole average
            changed_ = avr
        end

        if changed_ isa Average
            changed = inorder!(changed_)    # this can be done, since terms inside the sum commute anyway
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
function checkAndChangeDsum(missed_,sum,extras,states,checking,scaling;kwargs...)
    avrgs = getAvrgs(sum) #get vector of all avrgs in the sum
    meta = SymbolicUtils.metadata(sum)[IndexedAverageDoubleSum]
    inner = meta.innerSum
    outerInd = meta.sum_index
    meta_inner = SymbolicUtils.metadata(sum)[IndexedAverageSum]
    innerInd = meta_inner.sum_index
    for avr in avrgs
        changed_ = nothing
        avrg_inds = get_indices(avr)
        if outerInd in avrg_inds
            filtered = filterExtras(outerInd,extras)
            for i = 1:length(filtered)    # check if the extraindex is already in use for the term -> if so use the next in line
                if filtered[i] ∉ avrg_inds
                    changed_ = change_index(avr,outerInd,filtered[i])
                    break
                end
            end
        end
        if innerInd in avrg_inds
            filtered = filterExtras(innerInd,extras)
            for i = 1:length(filtered)
                if filtered[i] ∉ avrg_inds
                    changed_ = change_index(avr,innerInd,filtered[i])
                    break
                end
            end
        end
        if changed_ === nothing #if avrg consists of none indexed operators, or only of operators that already have the desired indices -> push along the whole average
            changed_ = avr
        end
        if changed_ isa Average
            changed = inorder!(changed_)    # this can be done, since terms inside the sum commute anyway
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
    return filter(x -> isequal(x.aon,wanted.aon),extras)
end

function find_missing!(missed,missed_hashes,term::BasicSymbolic{SpecialIndexedAverage},vhash,vshash;kwargs...)
    meta = SymbolicUtils.metadata(term)[SpecialIndexedAverage]
    return find_missing!(missed,missed_hashes,meta.term,vhash,vshash;kwargs...)
end

#Utility Functions
#checks if there is a sum (or multiple) in the equation rhs term, if so it returns the sums as a vector
function checkIfSum(term)
    sums = Any[]
    if term isa BasicSymbolic{IndexedAverageSum}
        return copy([term])
    elseif iscall(term)
        args = copy(arguments(term))
        for arg in args
            sums = vcat(sums,checkIfSum(arg))
        end
    end
    return sums
end
function getDSums(term)
    Dsums = Any[]
    if iscall(term)
        for arg in arguments(term)
            Dsums = vcat(Dsums,getDSums(arg))
        end
    elseif term isa BasicSymbolic{IndexedAverageDoubleSum}
        return [term]
    end
    return Dsums
end

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
function elimRed!(missed::Vector;scaling::Bool=false,kwargs...)
    counter = 2
    while counter <= length(missed)
        deleted = false
        for j=1:(counter-1)
            if scaling
                if isscaleequal(missed[counter],missed[j];kwargs...) || isscaleequal(inorder!(_inconj(missed[counter])),missed[j];kwargs...)
                    deleteat!(missed,counter)
                    deleted = true
                    break
                end
            else
                if isequal(missed[counter],missed[j]) || isequal(inorder!(_inconj(missed[counter])),missed[j])
                    deleteat!(missed,counter)
                    deleted = true
                    break
                end
            end
        end
        if !(deleted)
            counter = counter + 1
        end
    end
    return missed
end

complete(eqs::IndexedMeanfieldEquations;kwargs...) = indexed_complete(eqs;kwargs...)
complete!(eqs::IndexedMeanfieldEquations;kwargs...) = indexed_complete!(eqs;kwargs...)

"""
    evaluate(eqs::IndexedMeanfieldEquations;limits)
    evaluate(corr::CorrelationFunction;limits)
    evaluate(x;limits)

Function, that evaluates a given [`MeanfieldEquations`](@ref) entity and returns again equations,
where indices have been inserted and sums evaluated. Can also be called on individual terms and a [`CorrelationFunction`](@ref) entity, to
evaluate any summations inside these terms.

# Arguments
*`me::MeanfieldEquations`: A [`MeanfieldEquations`](@ref) entity, which shall be evaluated.

# Optional argumentes
*`limits::Dict{BasicSymbolic,Int64}=Dict{Symbol,Int64}()`: A seperate dictionary, to
    specify any symbolic limits used when [`Index`](@ref) entities were defined. This needs
    to be specified, when the equations contain summations, for which the upper bound is given
    by a Symbolic.
*`h`: A HilbertSpace, Vector of Hilbertspaces or Numbers, specifying the specific Hilbertspaces,
    that shall be evaluated. Does not evaluate any other Hilbertspace, other than the given ones.
"""
function evaluate(eqs::IndexedMeanfieldEquations;limits=nothing,h=nothing,kwargs...)
    hilb = hilbert(arguments(eqs[1].lhs)[1]) #hilbertspace of the whole system
    if !=(h,nothing)
        if !(h isa Vector)
            h=[h]
        end
        h_ = Vector{Any}(nothing,length(h))
        for i = 1:length(h)
            if h[i] isa HilbertSpace
                h_[i] = findfirst(x->isequal(x,h[i]),hilb.spaces)
            else
                h_[i] = h[i]
            end
        end
        h = h_
    end
    if !=(limits,nothing) && limits isa Pair
        limits_ = Dict{BasicSymbolic,Int64}(first(limits) => last(limits));
        limits = limits_
    end
    if limits === nothing
        limits =  Dict{BasicSymbolic,Int64}();
    end
    return subst_reds_eval(evalME(eqs;limits=limits,h=h,kwargs...);limits=limits,h=h,kwargs...)
end
function evaluate(term;limits=nothing,kwargs...)
    if !=(limits,nothing) && limits isa Pair
        limits_ = Dict{BasicSymbolic,Int64}(first(limits) => last(limits));
        limits = limits_
    end
    if limits === nothing
        limits =  Dict{BasicSymbolic,Int64}();
    end
    return eval_term(term;limits=limits)
end
