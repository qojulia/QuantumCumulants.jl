#File for adapting the meanfield, complete,... algorithms to indices. could be included in the already existing meanfield file.
#I did not want to change anything on the files in the already existing package, so I just copied and renamed a bunch of stuff here.

#include("doubleSums.jl")

#function that takes indexed operators and double indexed varaibles to calculate the meanfield equations
#the jump operators have to have same indices as the indices specified by the double indexed variable
function indexedMeanfield(a::Vector,H,J;Jdagger::Vector=adjoint.(J),rates=ones(Int,length(J)),
    multithread=false,
    simplify=true,
    order=nothing,
    mix_choice=maximum,
    iv=SymbolicUtils.Sym{Real}(:t))
    #print(a)

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

function indexedComplete(de::AbstractMeanfieldEquations;kwargs...)
    de_ = deepcopy(de)
    indexedComplete!(de_;kwargs...)
    return de_
end

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

    vhash = map(hash, vs)
    vs′ = map(_conj, vs)
    vs′hash = map(hash, vs′)
    filter!(!in(vhash), vs′hash)
    missed = find_missing(de.equations, vhash, vs′hash; get_adjoints=false)

    if order != 1
        missed = findMissingSumTerms(missed,de;extraIndices=extraIndices,scaling=scaling)
        missed = findMissingSpecialTerms(missed,de;scaling=scaling)
    end
    missed = sortByIndex.(missed)
    filter!(x -> (x isa Average),missed)

    isnothing(filter_func) || filter!(filter_func, missed) # User-defined filter
    
    filter!(x -> (isNotIn(getOps(x),getOps.(de.states),scaling) && isNotIn(getOps(sortByIndex(_conj(x))),getOps.(de.states),scaling)), missed)
    indices_ = nothing
    for i = 1:length(de.states)
        indices_ = getIndices(de.states[i])
        isempty(indices_) || break
    end
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

        if order != 1
            missed = findMissingSumTerms(missed,de;extraIndices=extraIndices,scaling=scaling)
            missed = findMissingSpecialTerms(missed,de;scaling=scaling)
        end
        missed = sortByIndex.(missed)
        filter!(x -> (x isa Average),missed)
        isnothing(filter_func) || filter!(filter_func, missed) # User-defined filter
    
        filter!(x -> (isNotIn(getOps(x),getOps.(de.states),scaling) && isNotIn(getOps(sortByIndex(_conj(x))),getOps.(de.states),scaling)), missed)
        for i = 1:length(missed)
            mInd_ = getIndices(missed[i])
            isempty(mInd_) && continue
            if indices_[1] ∉ mInd_ #term on lhs does not have the initial index -> change first occuring index into that one
                missed[i] = changeIndex(missed[i],mInd_[1],indices_[1]) #replace missed ops with changed indexed ones
            end
        end
        filter!(x -> (isNotIn(getOps(x),getOps.(de.states),scaling) && isNotIn(getOps(sortByIndex(_conj(x))),getOps.(de.states),scaling)), missed)
        missed = unique(missed)
        missed = elimRed(missed)
    end

    if !isnothing(filter_func)
        # Find missing values that are filtered by the custom filter function,
        # but still occur on the RHS; set those to 0
        missed = find_missing(de.equations, vhash, vs′hash; get_adjoints=false)
        if order != 1
            missed = findMissingSumTerms(missed,de;extraIndices=extraIndices,checking=false) #checkin dissabled, since there might be some terms, that are redundant, but also zero -> set them to zero aswell forsa fety
            missed = findMissingSpecialTerms(missed,de)
        end
        missed = sortByIndex.(missed) #this one might not be right (?)
        filter!(x -> (x isa Average),missed)

        filter!(!filter_func, missed)
    
        missed_adj = map(_adjoint, missed)
        subs = Dict(vcat(missed, missed_adj) .=> 0)
        for i=1:length(de.equations)
            de.equations[i] = substitute(de.equations[i], subs)
            de.states[i] = de.equations[i].lhs
        end
    end

    return de
end

#TODO: Clean up code, the functions below look and are very (!!!!) messy
# Function for extending find_missing function onto summation terms
function findMissingSumTerms(missed,de::MeanfieldEquations;extraIndices::Vector=[],checking=true,scaling=false)
    #TODO: simplify 3rd order to a loop over 2nd order algorithm (somehow)
    missed_ = copy(missed)
    indices = nothing #gets initial indices, that are on the lhs
    for i = 1:length(de.states)
        indices = getIndices(de.states[i])
        isempty(indices) || break
    end
    if de.order > 1 &&  de.order - 1 != length(extraIndices)
        error("Please make sure that for higher orders of cumulant expansion, you also use the indices argument to provide additional indices for calculation. You can do this by calling: indexedComplete(ME;extraIndices=[Vector of Symbols with length being (order-1)])")
    end
    extraIndex = nothing
    extras = []
    if isempty(indices)
        for op in de.jumps
            if typeof(op) == IndexedOperator
                extraIndex = op.ind
                push!(extras,extraIndex)
                break
            end
        end
    else
        for name in extraIndices
            push!(extras,Index(indices[1].hilb,name,indices[1].rangeN,indices[1].transition))
        end
    end
    #=
    if de.order == 1 #2nd order
        if !isempty(extras)
            extraIndex = extras[1]
        end
        sums = Any[]
        for eq in de.equations
            sums = checkIfSum(eq.rhs)
            for sum in sums
                avrgs = getAvrgs(sum) #get vector of all avrgs in the sum
                for avr in avrgs
                    changed_ = changeIndex(avr,sum.metadata.sumIndex,extraIndex)
                    if typeof(changed_) == Term{AvgSym, Nothing}
                        changed = sortByIndex(changed_)    # this can be done, since terms inside the sum commute anyway
                        if checking 
                            if ((isNotIn(getOps(changed),getOps.(de.states)) && isNotIn(getOps(sortByIndex(_conj(changed))),getOps.(de.states))) && 
                                    isNotIn(getOps(changed),getOps.(missed_)) && isNotIn(getOps(sortByIndex(_conj(changed))),getOps.(missed_)))
                                push!(missed_,changed)
                            end
                        else
                            push!(missed_,changed)
                        end
                    end
                end
            end
        end
    elseif de.order == 2 #third order
        =#
        
    sums = Any[]
    for eq in de.equations
        sums = checkIfSum(eq.rhs)
        for sum in sums
            avrgs = getAvrgs(sum) #get vector of all avrgs in the sum
            for avr in avrgs
                changed_ = nothing
                for ind in indices
                    if ind ∉ getIndices(avr)
                        changed_ = changeIndex(avr,sum.metadata.sumIndex,ind)
                        break
                    end
                end
                if isequal(changed_,nothing)
                    for i = 1:length(extras)
                        if extras[i] ∉ getIndices(avr)
                            changed_ = changeIndex(avr,sum.metadata.sumIndex,extras[i])
                            break
                        end
                    end
                end
                if typeof(changed_) == Term{AvgSym, Nothing}
                    changed = sortByIndex(changed_)    # this can be done, since terms inside the sum commute anyway
                    if !(changed isa Average)
                        continue
                    end
                    if checking
                        if ((isNotIn(getOps(changed),getOps.(de.states),scaling) && isNotIn(getOps(sortByIndex(_conj(changed))),getOps.(de.states),scaling)) && 
                                isNotIn(getOps(changed),getOps.(missed_),scaling) && isNotIn(getOps(sortByIndex(_conj(changed))),getOps.(missed_),scaling))
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
function findMissingSpecialTerms(missed,me::MeanfieldEquations;scaling=false)
    
    #does not work (?)
    #=

    missed_ = copy(missed)
    for eq in me.equations
        for arg in arguments(eq.rhs)
            if typeof(arg) == Int64
                continue
            end
            if typeof(arg) == SymbolicUtils.Sym{Parameter,SpecialIndexedAverage}
                for avr in getAvrgs(arg)
                    if ((isNotIn(getOps(avr),getOps.(me.states)) && isNotIn(getOps(sortByIndex(_conj(avr))),getOps.(me.states))) && 
                        isNotIn(getOps(avr),getOps.(missed_)) && isNotIn(getOps(sortByIndex(_conj(avr))),getOps.(missed_)))
                        push!(missed_,avr)
                    end
                end
            end
            if typeof(arg) <: SymbolicUtils.Mul
                for arg_ in arguments(arg)
                    if typeof(arg_) == Int64
                        continue
                    end
                    if typeof(arg_) == SymbolicUtils.Sym{Parameter,SpecialIndexedAverage}
                        for avr in getAvrgs(arg_)
                            if ((isNotIn(getOps(avr),getOps.(me.states)) && isNotIn(getOps(sortByIndex(_conj(avr))),getOps.(me.states))) && 
                                isNotIn(getOps(avr),getOps.(missed_)) && isNotIn(getOps(sortByIndex(_conj(avr))),getOps.(missed_)))
                                push!(missed_,avr)
                            end
                        end
                    end
                end
            end
        end
    end
    return missed_
    =#
    


    
    #old code:

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
function getOps(term::SymbolicUtils.Mul)
    for arg in arguments(term)
        if arg isa Average
            return getOps(arg)
        end
    end
    return []
end
function getOps(ops::Term{AvgSym, Nothing})
    args = arguments(ops)[1]
    if typeof(args) <: QMul 
        return getOps(args)
    else #single op in term
        return typeof(args) == IndexedOperator ? Any[args] : Any[args]
    end
end
function getOps(ops::QMul)
    arr = Any[]
    for arg in ops.args_nc
        if typeof(arg) == IndexedOperator
            push!(arr,arg)
        elseif typeof(arg) == NumberedOperator
            push!(arr,arg)
        else
            push!(arr,arg)
        end
    end
    return arr
end
getOps(x) = Vector{Vector{Any}}()
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
    for elem in vec1
        if elem ∉ vec2
            return false
        end
    end
    return true
end

# this is most likely not the best method, works fine however
function elimRed(missed::Vector)
    all = getOps.(missed)
    for i=1:length(missed)
        if i > length(missed)
            break
        end
        for j=1:length(missed)
            if j > length(missed)
                break
            end
            if hasSameOps(getOps(_conj(missed[i])),getOps(missed[j])) && !hasSameOps(getOps(missed[i]),getOps(missed[j]))
                deleteat!(missed,j)
            end
        end
    end
    return missed
end

function substReds(me::MeanfieldEquations)
    eqs = []
    for eq in me.equations
        push!(eqs,Symbolics.Equation(eq.lhs,substReds(eq.rhs,me.states)))
    end
    return MeanfieldEquations(eqs,me.operator_equations,me.states,me.operators,me.hamiltonian,me.jumps,me.jumps_dagger,me.rates,me.iv,me.varmap,me.order)
end

function substReds(term,vect::Vector)
    if term isa Average
        for i = 1:length(vect)
            temp = getOps(_conj(vect[i]))
            if hasSameOps(getOps(term),temp) && !(hasSameOps(getOps(term),getOps(vect[i])))
                return conj(vect[i])
            end
            if hasSameOps(getOps(term),getOps(vect[i]))
                return vect[i]
            end
        end
    elseif istree(term)
        op = operation(term)
        args = map(x->substReds(x,vect),arguments(term))
        return SymbolicUtils.similarterm(term, op, args)
    elseif typeof(term) == SymbolicUtils.Sym{Parameter,IndexedAverageSum}
        return IndexedAverageSum(substReds(term.metadata.term,vect),term.metadata.sumIndex,term.metadata.nonEqualIndices)
    elseif typeof(term) == SymbolicUtils.Sym{Parameter,IndexedAverageDoubleSum}
        return IndexedAverageDoubleSum(substReds(term.metadata.innerSum,vect),term.metadata.sumIndex,term.metadata.NEI)
    elseif typeof(term) == SymbolicUtils.Sym{Parameter,SpecialIndexedAverage}
        return SpecialIndexedAverage(substReds(term.metadata.term,vect),term.metadata.indexMapping)
    end
    return term
end

#=
if term isa Average
        for i = 1:length(vect)
            temp = getOps(_conj(vect[i]))
            if hasSameOps(getOps(term),temp)
                return vect[i]
            end
            if hasSameOps(getOps(term),getOps(vect[i]))
                return vect[i]
            end
        end
    elseif typeof(term) <: SymbolicUtils.Mul
        mults = []
        avrg = 1
        conj_ = false
        for arg in arguments(term)
            if arg isa Average
                for i = 1:length(vect)
                    temp = getOps(_conj(vect[i]))
                    if hasSameOps(getOps(term),temp) #is complex-conjugated of original
                        avrg = vect[i]
                        conj_ = true
                    end
                    if hasSameOps(getOps(term),getOps(vect[i])) #is same as original
                        avrg = vect[i]
                    end
                end
            else
                push!(mults,arg)
            end
        end
        if conj_
            for i=1:length(mults)
                if mults[i] isa Complex
                    mults[i] = conj(mults[i])
                end
            end
            return *(mults...)*avrg
        end
        return *(mults...)*avrg
        =#