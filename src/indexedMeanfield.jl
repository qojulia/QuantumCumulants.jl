#File for adapting the meanfield, complete,... algorithms to indices. could be included in the already existing meanfield file.
#I did not want to change anything on the files in the already existing package, so I just copied and renamed a bunch of stuff here.


#function that takes indexed operators and double indexed varaibles to calculate the meanfield equations
#the jump operators have to have same indices as the indices specified by the double indexed variable
"""
    indexedMeanfield(ops::Vector,H::QNumber,J::Vector;
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
function indexedMeanfield(a::Vector,H,J;Jdagger::Vector=adjoint.(J),rates=ones(Int,length(J)),
    multithread=false,
    simplify::Bool=true,
    order=nothing,
    mix_choice=maximum,
    iv=SymbolicUtils.Sym{Real}(:t))

    # Derive operator equations
    rhs = Vector{Any}(undef, length(a))
    imH = im*H
    for i=1:length(a)
        rhs_ = commutator(imH,a[i])
        rhs_diss = indexed_master_lindblad(a[i],J,Jdagger,rates)
        indices = getIndices(a[i])
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

    me = MeanfieldEquations(eqs_avg,eqs,vs,a,H,J,Jdagger,rates,iv,varmap,order)
    if has_cluster(H)
        return scale(me;simplify=simplify,order=order,mix_choice=mix_choice)
    else
        return me
    end
end

function indexed_master_lindblad(a_,J,Jdagger,rates)
    args = Any[]
    for k=1:length(J)
        if typeof(J[k]) == IndexedOperator
            c1 = orderByIndex(0.5*rates[k]*Jdagger[k]*commutator(a_,J[k]))
            c2 = orderByIndex(0.5*rates[k]*commutator(Jdagger[k],a_)*J[k])
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
                push!(args, IndexedSingleSum(c,J[k].ind,Index[]))
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
                push!(args,IndexedDoubleSum(IndexedSingleSum((c),J[k][1].ind,Index[]),J[k][2].ind,Index[]))
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
    indexedComplete(de::MeanfieldEquations)

From a set of differential equation of averages, find all averages that are missing
and derive the corresponding equations of motion. Uses [`find_missing`](@ref)
and [`indexedMeanfield`](@ref) to do so.

Optional arguments
==================

*`order=de.order`: The order at which the [`cumulant_expansion`](@ref) is performed
    on the newly derived equations. If `nothing`, the order is inferred from the
    existing equations.
*`filter_func=nothing`: Custom function that specifies whether some averages should
    be ignored when completing a system. This works by calling `filter!(filter_func, missed)`
    where `missed` is the vector resulting from [`find_missing`](@ref). Occurrences
    of averages for which `filter_func` returns `false` are substituted to 0.
*`extraIndices`: A Vector of symbols, representing extra [`Index`](@ref) entities, which are
    needed and created in the process of finding missing terms.
*`kwargs...`: Further keyword arguments are passed on to [`indexedMeanfield`](@ref) and
    simplification.

see also: [`find_missing`](@ref), [`indexedMeanfield`](@ref), [`meanfield`](@ref), [`findMissingSumTerms`](@ref)
"""
function indexedComplete(de::AbstractMeanfieldEquations;kwargs...)
    de_ = deepcopy(de)
    indexedComplete!(de_;kwargs...)
    return de_
end

"""
    indexedComplete!(de::MeanfieldEquations)

In-place version of [`indexedComplete`](@ref)
"""
function indexedComplete!(de::AbstractMeanfieldEquations;
    order=de.order,
    multithread=false,
    filter_func=nothing,
    mix_choice=maximum,
    simplify=true,
    extraIndices::Vector=[],
    scaling::Bool=false,
    kwargs...)
    sort!(extraIndices)
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

    indices_ = nothing
    for i = 1:length(de.states)
        indices_ = getIndices(de.states[i])
        isempty(indices_) || break
    end
    #08.08.22
    if isempty(indices_)
        for op in de.jumps
            if op isa IndexedOperator
                indices_ = [Index(op.ind.hilb,extraIndices[1],op.ind.rangeN,op.ind.specHilb)]
                deleteat!(extraIndices,1)
                break
            end
        end
    end

    vhash = map(hash, vs)
    vs′ = map(_conj, vs)
    vs′hash = map(hash, vs′)
    filter!(!in(vhash), vs′hash)
    missed = find_missing(de.equations, vhash, vs′hash; get_adjoints=false)

    
    missed = findMissingSumTerms(missed,de;extraIndices=extraIndices,scaling=scaling,indices=indices_)
    missed = findMissingSpecialTerms(missed,de;scaling=scaling)
    
    missed = sortByIndex.(missed)
    filter!(x -> (x isa Average),missed)

    isnothing(filter_func) || filter!(filter_func, missed) # User-defined filter
    
    filter!(x -> filterComplete(x,de.states,scaling), missed)
    

    sort!(indices_,by=getIndName)

    for i = 1:length(missed)
        mInd_ = getIndices(missed[i])
        isempty(mInd_) && continue
        if indices_[1] ∉ mInd_ #term on lhs does not have the initial index -> change first occuring index into that one
            missed[i] = changeIndex(missed[i],mInd_[1],indices_[1]) #replace missed ops with changed indexed ones
        end
    end
    missed = unique(missed) #no duplicates

    while !isempty(missed)
        ops_ = [SymbolicUtils.arguments(m)[1] for m in missed]
        me = indexedMeanfield(ops_,de.hamiltonian,de.jumps;
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

        missed = findMissingSumTerms(missed,de;extraIndices=extraIndices,scaling=scaling,indices=indices_)
        missed = findMissingSpecialTerms(missed,de;scaling=scaling)

        missed = sortByIndex.(missed)
        filter!(x -> (x isa Average),missed)
        isnothing(filter_func) || filter!(filter_func, missed) # User-defined filter
    
        filter!(x -> filterComplete(x,de.states,scaling), missed)
        for i = 1:length(missed)
            mInd_ = getIndices(missed[i])
            isempty(mInd_) && continue
            if indices_[1] ∉ mInd_ #term on lhs does not have the initial index -> change first occuring index into that one
                missed[i] = changeIndex(missed[i],mInd_[1],indices_[1]) #replace missed ops with changed indexed ones
            end
        end
        filter!(x -> filterComplete(x,de.states,scaling), missed)
        missed = unique(missed)
        missed = elimRed(missed;scaling=scaling)
    end

    if !isnothing(filter_func)
        # Find missing values that are filtered by the custom filter function,
        # but still occur on the RHS; set those to 0
        missed = find_missing(de.equations, vhash, vs′hash; get_adjoints=false)
        #if order != 1
            missed = findMissingSumTerms(missed,de;extraIndices=extraIndices,checking=false,indices=indices_) #checkin dissabled, since there might be some terms, that are redundant, but also zero -> set them to zero aswell forsa fety
            missed = findMissingSpecialTerms(missed,de)
        #end
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
filterComplete(x,states,scaling) = (isNotIn(getOps(x;scaling=scaling),getOps.(states;scaling=scaling),scaling) && isNotIn(getOps(sortByIndex(_conj(x));scaling=scaling),getOps.(states;scaling=scaling),scaling))

"""
    findMissingSumTerms(missed,de::MeanfieldEquations)

From a initial set of differential equation of averages, find all averages that are missing
and inside a Symbolic sum. If a missing average contains one of the summation indices used in
the equations, the [`Index`](@ref) will be exchanged according to the keyword argument
`extraIndices`. Uses [`find_missing`](@ref).

# Arguments
*`missed`: A initial Vector of averages, representing the missed averages before calling
    this method.
*`de`: The set of equations, in which the missing averages are searched in.

# Optional arguments
*`extraIndices`: A Vector of symbols, representing extra [`Index`](@ref) entities, which are
    needed and created in the process of finding missing terms. This argument is required, if the order
    of the Meanfield-Equations exceeds 1 and the number of given symbols must match the corresponding order.
*`checking`: A Bool defining if the algorithm checks for adjoint values and duplicates, before adding a found
    average into the `missed` vector.
*`scaling`: A Bool defining the way how averages are added to the `missed` vector. If true only averages, whose
    operators (without indices) are not already inside the `missed` vector will be added.

see also: [`find_missing`](@ref), [`indexedMeanfield`](@ref), [`meanfield`](@ref), [`findMissingSumTerms`](@ref)
"""
function findMissingSumTerms(missed,de::MeanfieldEquations;extraIndices::Vector=[],checking=true,scaling=false,indices::Vector=[])

    # This might not work if there are only double sums inside the meanfield equations
    # also: this might not work if double sums combine transition and bosonic operators -> NEEDS SOME CHECKING

    missed_ = copy(missed)

    if de.order > 1 &&  de.order - 1 > length(extraIndices)
        error("Wrong number of extraIndices! Please make sure that for higher orders of cumulant expansion, 
            you also use the extraIndices argument to provide additional indices for calculation. The Number of
            extraIndices provided should be at least (expansion order - 1).
        ")
        #TODO: change above notifcation and checks
    end

    extras = []

    for name in extraIndices    # create actual indices out of the extras provided
        push!(extras,Index(indices[1].hilb,name,indices[1].rangeN,indices[1].specHilb))
    end
    sums = Any[]
    for eq in de.equations
        sums = checkIfSum(eq.rhs)   #get all sums of the rhs
        for sum in sums
            avrgs = getAvrgs(sum) #get vector of all avrgs in the sum
            for avr in avrgs
                changed_ = nothing
                for ind in indices  #indices is a vector of indices already given freely by the system (e.g. by the operators on the lhs or the operators, used for creating the meanfieldEq)
                    if ind ∉ getIndices(avr)    # if one of these indices is free -> change the summation index into that one
                        changed_ = changeIndex(avr,sum.metadata.sumIndex,ind)
                        break
                    end
                end
                if isequal(changed_,nothing)
                    for i = 1:length(extras)    # check if the extraindex is already in use for the term -> if so use the next in line
                        if extras[i] ∉ getIndices(avr)
                            changed_ = changeIndex(avr,sum.metadata.sumIndex,extras[i])
                            break
                        end
                    end
                    if isequal(changed_,nothing) #if avrg consists of none indexed operators, or only of operators that already have the desired indices -> push along the whole average
                        changed_ = avr
                    end
                end
                if typeof(changed_) == Term{AvgSym, Nothing}
                    changed = sortByIndex(changed_)    # this can be done, since terms inside the sum commute anyway
                    if !(changed isa Average)
                        continue
                    end
                    if checking     # checks if the missed term found by this function is already in the list of missed averages (checking = false is needed for substituting and filtering of adjoint summation terms since the same function is called when substituting inside sums)
                        if filterComplete(changed,de.states,scaling) && filterComplete(changed,missed_,scaling)
                            push!(missed_,changed)
                        end
                    else
                        push!(missed_,changed)
                    end
                end
            end
        end
    end   
    return missed_
end
#TODO: write also function for double sums and findMissingDoubleSumTerms


#This is basically just the find_missing function but for the extra defined SpecialTerms, which have extra index constraints and for that also some extra checks
function findMissingSpecialTerms(missed,me::MeanfieldEquations;scaling=false)
    missed_ = copy(missed)
    vs = me.states
    vhash = map(hash, vs)
    vs′=map(_conj, vs)
    vs′hash = map(hash, vs′)
    filter!(!in(vhash), vs′hash)
    missed_hashes = map(hash,missed_)
    for eq in me.equations
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
function getOps(term::SymbolicUtils.Mul;scaling::Bool=false)
    for arg in arguments(term)
        if arg isa Average
            return getOps(arg;scaling=scaling)
        end
    end
    return []
end
function getOps(ops::Term{AvgSym, Nothing};scaling::Bool=false)
    args = arguments(ops)[1]
    if typeof(args) <: QMul 
        return getOps(args;scaling=scaling)
    elseif scaling && (args isa NumberedOperator || args isa IndexedOperator)
        return Any[args.op]
    else #single op in term
        return Any[args]
    end
end
function getOps(ops::QMul;scaling::Bool=false)
    arr = Any[]
    for arg in ops.args_nc
        if scaling && (arg isa NumberedOperator || arg isa IndexedOperator)
            push!(arr,arg.op)
        else
            push!(arr,arg)
        end
    end
    return arr
end
getOps(x;scaling::Bool=false) = Vector{Vector{Any}}()
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
function elimRed(missed::Vector;scaling::Bool=false)
    all = getOps.(missed)
    for i=1:length(missed)
        if i > length(missed)
            break
        end
        for j=1:length(missed)
            if j > length(missed)
                break
            end
            if hasSameOps(getOps(_conj(missed[i]);scaling=scaling),getOps(missed[j];scaling=scaling)) && !hasSameOps(getOps(missed[i];scaling=scaling),getOps(missed[j];scaling=scaling))
                deleteat!(missed,j)
            end
        end
    end
    return missed
end

"""
    substReds(de::MeanfieldEquations)

Function that substitutes possible redundant conjugate averages inside the given Equations with their corresponding
average given as the conjugate of one of the left-hand-side (of the equations) averages.

# Optional Arguments
*`scaling`: A Bool defining the way how averages are added to the `missed` vector. If true only averages, whose
    operators (without indices) are not already inside the `missed` vector will be added.
"""
function substReds(me::MeanfieldEquations;scaling::Bool=false)
    eqs = []
    for eq in me.equations
        push!(eqs,Symbolics.Equation(eq.lhs,substReds(eq.rhs,me.states;scaling=scaling)))
    end
    return MeanfieldEquations(eqs,me.operator_equations,me.states,me.operators,me.hamiltonian,me.jumps,me.jumps_dagger,me.rates,me.iv,me.varmap,me.order)
end
function substReds(term,vect::Vector;scaling::Bool=false)
    if term isa Average
        for i = 1:length(vect)
            temp = getOps(_conj(vect[i]);scaling=scaling)
            if hasSameOps(getOps(term;scaling=scaling),temp) && !(hasSameOps(getOps(term;scaling=scaling),getOps(vect[i];scaling=scaling)))
                return conj(vect[i])
            end
            if hasSameOps(getOps(term;scaling=scaling),getOps(vect[i];scaling=scaling))
                return vect[i]
            end
        end
    elseif istree(term)
        op = operation(term)
        args = map(x->substReds(x,vect;scaling=scaling),arguments(term))
        return SymbolicUtils.similarterm(term, op, args)
    elseif typeof(term) == SymbolicUtils.Sym{Parameter,IndexedAverageSum}
        return IndexedAverageSum(substReds(term.metadata.term,vect;scaling=scaling),term.metadata.sumIndex,term.metadata.nonEqualIndices)
    elseif typeof(term) == SymbolicUtils.Sym{Parameter,IndexedAverageDoubleSum}
        return IndexedAverageDoubleSum(substReds(term.metadata.innerSum,vect;scaling=scaling),term.metadata.sumIndex,term.metadata.NEI)
    elseif typeof(term) == SymbolicUtils.Sym{Parameter,SpecialIndexedAverage}
        return SpecialIndexedAverage(substReds(term.metadata.term,vect;scaling=scaling),term.metadata.indexMapping)
    end
    return term
end