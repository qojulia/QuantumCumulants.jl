# indexed_arithmetic(a_,J,Jdagger,rate, arithmetic)

# This function tries to perform a summation over some arithmetic expression
# given three operators a_, J and Jdagger and numbers/symbols/a matrix rate.

# #Example
# The terms that occur due to a single jump operator J in the master equation can
# be generated with (for the equation of motion for operator a_)

# master_arithmetic = (a, J, Jdagger,rate) -> 0.5*rate*Jdagger*commutator(a,J)+0.5*rate*commutator(Jdagger,a)*J
# indexed_arithmetic(a_, Jloc, Jdaggerloc,rate, master_arithmetic)
function indexed_arithmetic(a_, J, Jdagger, rate, arithmetic)
    inds1 = get_indices(J)
    inds2 = get_indices(Jdagger)
    rinds = get_indices(rate)
    all_indices = collect(Set(vcat(inds1, inds2, rinds)))

    if isnothing(all_indices) || length(all_indices) == 0#No indices

        if isa(rate, Matrix)
            args = Any[]
            for i = 1:length(J[k]), j = 1:length(J[k])
                push_or_append_nz_args!(args, arithmetic(a_, J, Jdagger, rate))
            end
            return QAdd(args)
        elseif isa(rate, SymbolicUtils.Symbolic) || isa(rate, Number) || isa(rate, Function)
            return arithmetic(a_, J, Jdagger, rate)
        else
            error("Unknown rates type!")
        end

    elseif length(all_indices) == 1#Single index
        return ∑(arithmetic(a_, J, Jdagger, rate), all_indices[1])
    elseif length(all_indices) == 2# Double index
        return ∑(arithmetic(a_, J, Jdagger, rate), all_indices[1], all_indices[2])
    else
        error("Too many indices occurring in multiplication")
    end
end

# Indexed master equation calculation with the refactored functionality
function _simplified_indexed_master(a_, J, Jdagger, rates)
    out = nothing
    master_arithmetic =
        (a, J, Jdagger, rate) ->
            0.5*rate*Jdagger*commutator(a, J)+0.5*rate*commutator(Jdagger, a)*J

    for k = 1:length(J)

        # Construct the proper decay operators
        rate = rates[k]
        Jloc = J[k]
        Jdaggerloc = Jdagger[k]

        # For a vector user input just assign the proper decay operators
        if J[k] isa Vector
            Jloc = Jloc[2]
            Jdaggerloc = Jdaggerloc[1]
        elseif !(isempty(get_indices(J[k]))) # Define decay operator but with different index if no vector was given
            inds = get_indices(Jloc)
            length(inds) != 1 &&
                error("Jump operators with multiple different indices are not supported!")
            ind_ = inds[1]
            rind = get_indices(rate)
            if !isequal(ind_, rind[1]) && !isequal(ind_, rind[2])
                error("unequal index of jump operator and variable!")
            elseif isequal(ind_, rind[1]) && !isequal(ind_, rind[2])
                vec = [J[k], change_index(Jloc, ind_, rind[2])]
            elseif !isequal(ind_, rind[1]) && isequal(ind_, rind[2])
                vec = [change_index(Jloc, ind_, rind[2]), J[k]]
            end
            Jloc = vec[2]
            Jdaggerloc = adjoint(vec[1])
        end

        res = indexed_arithmetic(a_, Jloc, Jdaggerloc, rate, master_arithmetic)
        if out===nothing
            out = res
        else
            out+=res
        end
    end
    return out
end

function indexed_noise(a_, J, Jdagger, rates, efficiencies)
    out = nothing
    noise_arithmetic(a, J, Jdagger, rate) =
        sqrt(rate)*(Jdagger*a-average(Jdagger)*a) + sqrt(rate)*(a*J-a*average(J))

    for k = 1:length(J)
        if isequal(efficiencies[k], 0)
            continue
        end
        rate = rates[k] * efficiencies[k]
        all_indices = collect(
            Set(vcat(get_indices(J[k]), get_indices(Jdagger[k]), get_indices(rate))),
        )

        if !isnothing(all_indices) &&
           !isnothing(get_indices(rate)) &&
           length(get_indices(rate)) > 1
            error("Rates with multiple indices not supported for measurement backaction")
        end

        res = indexed_arithmetic(a_, J[k], Jdagger[k], rate, noise_arithmetic)
        if out===nothing
            out = res
        else
            out+=res
        end
    end

    return out
end

# This function generates a self mapping of indices onto itself without the diagonal, i.e.
# (i,j,k) would yield the map [(i,j), (i,k),(j,i),(j,k),(k,i),(k,j)]
function generate_index_self_mapping(indices)
    mapping = Tuple{Index,Index}[]
    for j = 1:length(indices)
        for k = 1:j
            if k != j
                push!(mapping, (indices[k], indices[j]))
            end
        end
    end
    return mapping
end

# Creates an array of equations from two arrays of left and righthand sides for equations
create_equation_array(lhs, rhs) = [Symbolics.Equation(l, r) for (l, r) in zip(lhs, rhs)]

function check_index_collision(a::Vector, H, J)
    for ind in get_indices(a)
        if ind in get_indices(H)
            error("Index $(ind.name) in operator-vector is already used in H!")
        end
        if ind in get_indices(J)
            error("Index $(ind.name) in operator-vector is already used in J!")
        end
    end
end

"""
    IndexedMeanfieldNoiseEquations

Like a [`MeanfieldNoiseEquations`](@ref), but with symbolic indices.
"""
struct IndexedMeanfieldNoiseEquations <: AbstractMeanfieldEquations
    equations::Vector{Symbolics.Equation}
    operator_equations::Vector{Symbolics.Equation}
    noise_equations::Vector{Symbolics.Equation}
    operator_noise_equations::Vector{Symbolics.Equation}
    states::Vector
    operators::Vector{QNumber}
    hamiltonian::QNumber
    jumps::Vector
    jumps_dagger::Any
    rates::Vector
    efficiencies::Vector
    iv::MTK.Num
    varmap::Vector{Pair}
    order::Union{Int,Vector{<:Int},Nothing}
end

function _append!(lhs::IndexedMeanfieldNoiseEquations, rhs::IndexedMeanfieldNoiseEquations)
    append!(lhs.noise_equations, rhs.noise_equations)
    append!(lhs.operator_noise_equations, rhs.operator_noise_equations)
    append!(lhs.equations, rhs.equations)
    append!(lhs.operator_equations, rhs.operator_equations)
    append!(lhs.states, rhs.states)
    append!(lhs.operators, rhs.operators)
    append!(lhs.varmap, rhs.varmap)
end

"""
indexed_meanfield_backaction(ops::Vector,H::QNumber,J::Vector;
    Jdagger::Vector=adjoint.(J),rates=ones(length(J)),efficiencies=zeros(Int,length(J)))

Compute the set of equations for the indexed-operators [`IndexedOperator`](@ref) in `ops` under the Hamiltonian
`H` and with loss operators contained in `J` and measurement backaction due to detectors with efficiencies 'efficiencies'.
The resulting equation is equivalent to the Quantum-Langevin equation where noise is neglected.
This is a modified version of the [`indexed_meanfield`](@ref) function, that also generates the measurement backaction
terms if nonzero efficiencies are given.
See also: [`indexed_meanfield`](@ref).

# Arguments
*`ops::Vector`: The operators of which the equations are to be computed.
*`H::QNumber`: The Hamiltonian describing the reversible dynamics of the
    system.
*`J::Vector{<:QNumber}`: A vector containing the collapse operators of
    the system. A term of the form
    ``\\sum_i J_i^\\dagger O J_i - \\frac{1}{2}\\left(J_i^\\dagger J_i O + OJ_i^\\dagger J_i\\right)``
    is added to the Heisenberg equation.

# Optional arguments
*`Jdagger::Vector=adjoint.(J)`: Vector containing the hermitian conjugates of
    the collapse operators.
*`rates=ones(length(J))`: Decay rates corresponding to the collapse operators in `J`.
*`efficiencies=zeros(length(J))`: Efficiencies for the operators measuring the quantum jump associated with 'J'.
*`multithread=false`: Specify whether the derivation of equations for all operators in `ops`
    should be multithreaded using `Threads.@threads`.
*`simplify=true`: Specify whether the derived equations should be simplified.
*`order=nothing`: Specify to which `order` a [`cumulant_expansion`](@ref) is performed.
    If `nothing`, this step is skipped.
*`mix_choice=maximum`: If the provided `order` is a `Vector`, `mix_choice` determines
    which `order` to prefer on terms that act on multiple Hilbert spaces.
*`iv=ModelingToolkit.t`: The independent variable (time parameter) of the system.
"""
function indexed_meanfield_backaction(
    a::Vector,
    H,
    J;
    Jdagger::Vector = adjoint.(J),
    rates = ones(Int, length(J)),
    efficiencies = zeros(Int, length(J)),
    multithread = false,
    simplify::Bool = true,
    order = nothing,
    mix_choice = maximum,
    iv = MTK.t_nounits,
    kwargs...,
)

    check_index_collision(a, H, J)

    # Derive operator equations
    rhs = Vector{Any}(undef, length(a))
    rhs_noise = Vector{Any}(undef, length(a))
    imH = im*H

    function calculate_term(i)
        try
            rhs[i] =
                commutator(imH, a[i]) + indexed_master_lindblad(a[i], J, Jdagger, rates)
            rhs_noise[i] = indexed_noise(a[i], J, Jdagger, rates, efficiencies)
            indices = get_indices(a[i])
            if length(indices) > 1 #everything on lhs commutes -> reorder corresponding terms on rhs
                rhs[i] = reorder(rhs[i], generate_index_self_mapping(indices))
                rhs_noise[i] = reorder(rhs_noise[i], generate_index_self_mapping(indices))
            end
        catch err
            println("could not calculate meanfield-equations for operator $(a[i])")
            rethrow(err)
        end
    end

    if multithread
        Threads.@threads for i = 1:length(a)
            calculate_term(i)
        end
    else
        for i = 1:length(a)
            calculate_term(i)
        end
    end

    # Average
    vs = map(average, a)
    rhs_avg, rhs = take_function_averages(rhs, simplify)
    rhs_noise_avg, rhs_noise = take_function_averages(rhs_noise, simplify)
    eqs_avg = create_equation_array(vs, rhs_avg)
    eqs = create_equation_array(a, rhs)
    eqs_noise_avg = create_equation_array(vs, rhs_noise_avg)
    eqs_noise = create_equation_array(a, rhs_noise)

    cumulant_expand_equations!(eqs_avg, order; mix_choice = mix_choice, simplify = simplify)
    cumulant_expand_equations!(
        eqs_noise_avg,
        order;
        mix_choice = mix_choice,
        simplify = simplify,
    )

    varmap = make_varmap(vs, iv)

    me = IndexedMeanfieldNoiseEquations(
        eqs_avg,
        eqs,
        eqs_noise_avg,
        eqs_noise,
        vs,
        a,
        H,
        J,
        Jdagger,
        rates,
        efficiencies,
        iv,
        varmap,
        order,
    )
    return me
end

# This function makes some sanity checks on the input extra indices
function check_extra_indices(extra_indices)
    if isempty(extra_indices)
        error("can not complete equations with empty extra_indices!")
    end

    for ind in extra_indices
        if typeof(ind) != typeof(extra_indices[1])
            error(
                "Cannot use extra_indices of different types. Use either only Symbols or Index-Objects!",
            )
        end
    end
end

# filters indices from extra_indices which are contained in allInds
function filer_indices!(extra_indices, allInds)
    if extra_indices[1] isa Symbol
        filter!(x -> x ∉ SQA.getIndName.(allInds), extra_indices)
    else
        filter!(x -> x ∉ allInds, extra_indices)
    end
end

# determines the current and maximum order that is to be calculated using the current system
function determine_orders(vs, eqs; order = de.order)
    order_lhs = maximum(get_order.(vs))
    order_rhs = 0

    # Determine RHS order
    for i = 1:length(eqs)
        k = get_order(eqs[i].rhs)
        k > order_rhs && (order_rhs = k)
    end

    # Set order to max of LHS and RHS if not provided
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

    maximum(order_) >= order_lhs || error(
        "Cannot form cumulant expansion of derivative; you may want to use a higher order!",
    )

    return order, order_max
end

# $his generates a list of extra indices ordered by priority
function generate_extras(vs, extra_indices, allInds, order_max)

    indices_lhs = get_all_indices(vs) #indices that are used in the beginning -> have priority

    if containsMultiple(allInds) && extra_indices[1] isa Symbol
        dic = split_inds(allInds)
        filter!(x->!isempty(last(x)), dic)
        dic2 = Dict{Int,Any}(i => Any[] for i in keys(dic))
        for k in keys(dic)
            dic2[k] = filter(x->isequal(x.aon, k), indices_lhs)
            ind = allInds[findfirst(x->isequal(x.aon, k), allInds)]
            while length(dic2[k]) < order_max
                push!(dic2[k], Index(ind.hilb, extra_indices[1], ind.range, k))
                deleteat!(extra_indices, 1)
            end
        end
        extras = []
        for vec in values(dic2)
            extras = [extras; vec]
        end
    else
        #if there are no indices on the LHS -> use the first index of extra_indices as starting point
        #this first index is (if it is a symbol) created similar to the first occurring index in the Equations
        if isempty(indices_lhs) && extra_indices[1] isa Symbol
            for ind in allInds
                indices_lhs = [Index(ind.hilb, extra_indices[1], ind.range, ind.aon)]
                deleteat!(extra_indices, 1)
                break
            end
        elseif isempty(indices_lhs) && extra_indices[1] isa Index
            indices_lhs = [extra_indices[1]]
            deleteat!(extra_indices, 1)
        end

        #create (or append to the extras vector) all other extra_indices using the indices on the LHS
        extras = indices_lhs
        if !isempty(extra_indices)
            if extra_indices[1] isa Symbol
                first = extras[1]
                for name in extra_indices
                    push!(extras, Index(first.hilb, name, first.range, first.aon))
                    if length(extras) >= order_max
                        break
                    end
                end
            end
            if extra_indices[1] isa Index
                extras = [extras; extra_indices]
            end
        end
        if length(extras) < order_max
            error("More extra_indices are needed for order $(order_max)")
        end

    end

    return extras
end


function find_missing_sums(
    missed,
    eqs,
    vs;
    extra_indices::Vector = [],
    checking = true,
    scaling = false,
    kwargs...,
)
    missed_ = copy(missed)
    extras = copy(extra_indices)
    for eq in eqs
        sums = checkIfSum(eq.rhs)   #get all sums of the rhs
        for sum in sums
            meta = TermInterface.metadata(sum)[IndexedAverageSum]
            extras_filtered = filterExtras(meta.sum_index, extras) #all the extraIndices, that have the same specific hilbertspace as the sum
            checkAndChange(missed_, sum, extras_filtered, vs, checking, scaling; kwargs...)
        end
    end
    return inorder!.(missed_)
end

#this function might not be so fast
function find_missing_Dsums(
    missed,
    eqs;
    extra_indices::Vector = [],
    checking = true,
    scaling = false,
    kwargs...,
)
    missed_ = copy(missed)
    extras = copy(extra_indices)
    for eq in eqs
        Dsums = getDSums(eq.rhs)
        for sum in Dsums
            checkAndChangeDsum(
                missed_,
                sum,
                extras,
                de.states,
                checking,
                scaling;
                kwargs...,
            )
        end
    end
    return inorder!.(missed_)
end

# Determines which operators are missing from the system and the updated system and filters all redundant/unneeded operators
function find_index_missed(vs, eqs_de, eqs_me, extras, filter_func, scaling; kwargs...)
    vhash = map(hash, vs)
    vs′ = map(_inconj, vs)
    vs′hash = map(hash, vs′)
    filter!(!in(vhash), vs′hash)


    missed = find_missing(eqs_me, vhash, vs′hash; get_adjoints = false)
    missed = find_missing_sums(missed, eqs_de, vs; extra_indices = extras, kwargs...)
    missed = find_missing_Dsums(missed, eqs_de; extra_indices = extras, kwargs...)

    missed = inorder!.(missed)

    filter!(x -> (x isa Average), missed)
    isnothing(filter_func) || filter!(filter_func, missed) # User-defined filter

    filter!(x -> filterComplete(x, vs, scaling; kwargs...), missed)
    missed = inorder!.(missed)

    for i = 1:length(missed)
        minds = get_indices(missed[i])
        newMinds = copy(minds)
        for ind1 in minds
            extras_=filterExtras(ind1, extras)
            for k = 1:length(extras_)
                if ind1 ∉ extras_ #this check might be redundant ?
                    println(typeof(missed[i]))
                    println(typeof(change_index(missed[i], ind1, extras_[1])))

                    missed[i] = change_index(missed[i], ind1, extras_[1])
                    newMinds = get_indices(missed[i])
                    break
                end
                if findall(x->isequal(x, ind1), extras_)[1] > k && extras_[k] ∉ newMinds #this might go somewhat easier, maybe delete ind2 out of extras after each replacement somehow
                    missed[i] = change_index(missed[i], ind1, extras_[k])
                    newMinds = get_indices(missed[i])
                    break
                elseif findall(x->isequal(x, ind1), extras_)[1] <= k
                    break
                end
            end
        end
    end

    missed = inorder!.(missed)
    missed = unique(missed)
    missed = filter!(x -> filterComplete(x, vs, scaling; kwargs...), missed)
    missed = inorder!.(missed)
    missed = elimRed!(missed; scaling = scaling, kwargs...)

    return missed
end

# Generates the substitutions such that the operators included in filter_func that occur in eqs or vs are set to zero
function generate_substitutions(vs, eqs, filter_func, extras)
    vhash = map(hash, vs)
    vs′ = map(_inconj, vs)
    vs′hash = map(hash, vs′)
    # Find missing values that are filtered by the custom filter function,
    # but still occur on the RHS; set those to 0
    missed = find_missing(eqs, vhash, vs′hash; get_adjoints = false)
    missed =
        find_missing_sums(missed, de; extra_indices = extras, checking = false, kwargs...) #checkin disabled, since there might be some terms, that are redundant, but also zero -> set them to zero aswell forsa fety
    missed =
        find_missing_Dsums(missed, de; extra_indices = extras, checking = false, kwargs...)
    missed = find_missing_sums(missed, de; checking = false, kwargs...)
    missed = find_missing_Dsums(missed, de; checking = false, kwargs...)

    missed = inorder!.(missed) #this one might not be right (?) -> not even needed i think
    filter!(x -> (x isa Average), missed)

    filter!(!filter_func, missed)
    missed_adj = map(_inconj, missed)
    subs = Dict(vcat(missed, missed_adj) .=> 0)
    return subs
end

# Cumulant expand all equations in eqs
function cumulant_expand_equations!(eqs, order; mix_choice = maximum, simplify = true)
    for i = 1:length(eqs)
        lhs = eqs[i].lhs
        rhs = cumulant_expansion(
            eqs[i].rhs,
            order;
            mix_choice = mix_choice,
            simplify = simplify,
        )
        eqs[i] = Symbolics.Equation(lhs, rhs)
    end
end

function simplified_indexed_complete!(
    de::AbstractMeanfieldEquations;
    order = de.order,
    multithread = false,
    filter_func = nothing,
    mix_choice = maximum,
    simplify = true,
    extra_indices::Vector = [:i, :j, :k, :l, :m, :n, :p, :q, :r, :s, :t],
    scaling = false,
    kwargs...,
)

    check_extra_indices(extra_indices)
    maxNumb = maximum(length.(get_indices.(de.operators)))
    allInds = get_all_indices(de)
    filer_indices!(extra_indices, allInds)
    vs = de.states
    order_, order_max = determine_orders(vs, de.equations; order)

    if order_max > maxNumb && order_max - maxNumb > length(extra_indices)
        error(
            "Too few extra_indices provided! Please make sure that for higher orders of cumulant expansion,
          you also use the extra_indices argument to provide additional indices for calculation. The Number of
          extra_indices provided should be at least $(order_max - maxNumb).
      ",
        )
    end

    if order_ != de.order
        cumulant_expand_equations!(
            de.equations,
            order_;
            mix_choice = mix_choice,
            simplify = simplify,
        )
    end

    extras = generate_extras(vs, extra_indices, allInds, order_max)
    missed = find_index_missed(
        vs,
        de.equations,
        de.equations,
        extras,
        filter_func,
        scaling;
        kwargs...,
    )

    while !isempty(missed)
        ops_ = [SymbolicUtils.arguments(m)[1] for m in missed]
        me = indexed_meanfield(
            ops_,
            de.hamiltonian,
            de.jumps;
            Jdagger = de.jumps_dagger,
            rates = de.rates,
            simplify = simplify,
            multithread = multithread,
            order = order_,
            mix_choice = mix_choice,
            iv = de.iv,
            kwargs...,
        )
        _append!(de, me)
        missed = find_index_missed(
            vs,
            de.equations,
            me.equations,
            extras,
            filter_func,
            scaling;
            kwargs...,
        )
    end

    if !isnothing(filter_func)
        subs = generate_substitutions(vs, de.equations, filter_func, extras)
        for i = 1:length(de.equations)
            de.equations[i] = (SymbolicUtils.substitute(de.equations[i], subs))
            de.states[i] = de.equations[i].lhs
        end
    end

    return de
end


"""
    indexed_complete(de::IndexedMeanfieldNoiseEquations)

From a set of differential equation of averages, find all averages that are missing
and derive the corresponding equations of motion. Uses [`find_missing`](@ref)
and [`indexed_meanfield`](@ref) to do so. Implementation for IndexedMeanfieldNoiseEquations.

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
function indexed_complete(de::IndexedMeanfieldNoiseEquations; kwargs...)
    de_ = deepcopy(de)
    indexed_complete!(de_; kwargs...)
    return de_
end


"""
    indexed_complete!(de::MeanfieldEquations)

In-place version of [`indexed_complete`](@ref)
"""
function indexed_complete!(
    de::IndexedMeanfieldNoiseEquations;
    order = de.order,
    multithread = false,
    filter_func = nothing,
    mix_choice = maximum,
    simplify = true,
    extra_indices::Vector = [:i, :j, :k, :l, :m, :n, :p, :q, :r, :s, :t],
    scaling = false,
    kwargs...,
)

    check_extra_indices(extra_indices)
    maxNumb = maximum(length.(get_indices.(de.operators)))
    allInds = get_all_indices(de)
    filer_indices!(extra_indices, allInds)
    vs = de.states
    order_, order_max = determine_orders(vs, vcat(de.equations, de.noise_equations); order)

    if order_max > maxNumb && order_max - maxNumb > length(extra_indices)
        error(
            "Too few extra_indices provided! Please make sure that for higher orders of cumulant expansion,
          you also use the extra_indices argument to provide additional indices for calculation. The Number of
          extra_indices provided should be at least $(order_max - maxNumb).
      ",
        )
    end

    if order_ != de.order
        cummulant_expand_equations!(
            de.equations,
            order_;
            mix_choice = mix_choice,
            simplify = simplify,
        )
        cummulant_expand_equations!(
            de.noise_equations,
            order_;
            mix_choice = mix_choice,
            simplify = simplify,
        )
    end

    extras = generate_extras(vs, extra_indices, allInds, order_max)
    missed = find_index_missed(
        vs,
        vcat(de.equations, de.noise_equations),
        vcat(de.equations, de.noise_equations),
        extras,
        filter_func,
        scaling;
        kwargs...,
    )
    #missed = find_index_missed(vs, de.equations, de.equations, extras, filter_func, scaling;kwargs...)

    while !isempty(missed)

        ops_ = [SymbolicUtils.arguments(m)[1] for m in missed]
        me = indexed_meanfield_backaction(
            ops_,
            de.hamiltonian,
            de.jumps;
            Jdagger = de.jumps_dagger,
            rates = de.rates,
            efficiencies = de.efficiencies,
            simplify = simplify,
            multithread = multithread,
            order = order_,
            mix_choice = mix_choice,
            iv = de.iv,
            kwargs...,
        )
        _append!(de, me)
        missed = find_index_missed(
            vs,
            vcat(de.equations, de.noise_equations),
            vcat(me.equations, me.noise_equations),
            extras,
            filter_func,
            scaling;
            kwargs...,
        )
    end

    if !isnothing(filter_func)
        subs = generate_substitutions(vs, de.equations, filter_func, extras)
        subs_noise = generate_substitutions(vs, de.noise_equations, filter_func, extras)
        for i = 1:length(de.equations)
            de.equations[i] = (SymbolicUtils.substitute(de.equations[i], subs))
            de.noise_equations[i] =
                (SymbolicUtils.substitute(de.noise_equations[i], subs_noise))
            de.states[i] = de.equations[i].lhs
        end
    end

    return de
end

MTK.complete(eqs::IndexedMeanfieldNoiseEquations; kwargs...) =
    indexed_complete(eqs; kwargs...)
MTK.complete!(eqs::IndexedMeanfieldNoiseEquations; kwargs...) =
    indexed_complete!(eqs; kwargs...)


function scale_equations(eqs; kwargs...)
    newEqs = []
    for eq in eqs
        tempEq = scaleEq(eq; kwargs...)
        if tempEq.lhs in getLHS.(newEqs)
            continue
        elseif isNotIn(tempEq.lhs, getLHS.(newEqs), true; kwargs...) &&
               isNotIn(_inconj(tempEq.lhs), getLHS.(newEqs), true; kwargs...)
            push!(newEqs, tempEq)
        end
    end

    return newEqs
end

"""
    scaleME(me::IndexedMeanfieldNoiseEquations)

Function, that evaluates a given [`IndexedMeanfieldNoiseEquations`](@ref) entity and returns again equations,
where indices have been inserted and sums evaluated, regarding the same relations, as done when calculating
with oparators using a [`ClusterSpace`](@ref).

# Arguments
*`me::IndexedMeanfieldNoiseEquations`: A [`IndexedMeanfieldNoiseEquations`](@ref) entity, which shall be evaluated.

"""
function scaleME(me::IndexedMeanfieldNoiseEquations; kwargs...)
    newEqs = scale_equations(me.equations; kwargs...)
    newNoiseEqs = scale_equations(me.noise_equations; kwargs...)
    vs = getLHS.(newEqs)
    varmap = make_varmap(vs, me.iv)
    ops = undo_average.(vs)
    return IndexedMeanfieldNoiseEquations(
        newEqs,
        me.operator_equations,
        newNoiseEqs,
        me.operator_noise_equations,
        vs,
        ops,
        me.hamiltonian,
        me.jumps,
        me.jumps_dagger,
        me.rates,
        me.efficiencies,
        me.iv,
        varmap,
        me.order,
    )
end

# substitute redundant operators in scaled equations
function subst_reds_scale_equation(states, eqs; kwargs...)

    #states = me.states
    vhash = map(hash, states)
    states′ = map(_inconj, states)
    vs′hash = map(hash, states′)
    missed = find_missing(eqs, vhash, vs′hash; get_adjoints = false)
    inorder!.(missed)

    to_sub = missed
    filter!(x->x ∉ states, to_sub)

    to_insert = Vector{Any}(nothing, length(to_sub))

    counter = 1
    while counter <= length(to_sub)
        elem = to_sub[counter]
        ind_ = findfirst(x -> isscaleequal(elem, x; kwargs...), states)
        if !=(ind_, nothing)
            to_insert[counter] = states[ind_]
            counter = counter + 1
        else
            ind_ =
                findfirst(x -> isscaleequal(inorder!(_inconj(elem)), x; kwargs...), states)
            if !=(ind_, nothing)
                to_insert[counter] = conj(states[ind_])
                counter = counter + 1
            else
                deleteat!(to_insert, counter) # these deletes are for consistency only -> it is possible that not all terms are fully evaluated
                deleteat!(to_sub, counter)   # yet in the system -> leftovers in the find_missing
            end
        end
    end

    subs = Dict(to_sub .=> to_insert)
    return [substitute(eq, subs) for eq in eqs]
end

function scale(eqs::IndexedMeanfieldNoiseEquations; h = nothing, kwargs...)
    hilb = hilbert(arguments(eqs[1].lhs)[1]) #hilbertspace of the whole system
    if !=(h, nothing)
        if !(h isa Vector)
            h=[h]
        end
        h_ = Vector{Any}(nothing, length(h))
        for i = 1:length(h)
            if h[i] isa HilbertSpace
                h_[i] = findfirst(x->isequal(x, h[i]), hilb.spaces)
            else
                h_[i] = h[i]
            end
        end
        h = h_
    end

    me = scaleME(eqs; h = h, kwargs...)
    newEqs = subst_reds_scale_equation(me.states, me.equations; h = h, kwargs...)
    newNoiseEqs = subst_reds_scale_equation(me.states, me.noise_equations; h = h, kwargs...)
    return IndexedMeanfieldNoiseEquations(
        newEqs,
        me.operator_equations,
        newNoiseEqs,
        me.operator_noise_equations,
        me.states,
        me.operators,
        me.hamiltonian,
        me.jumps,
        me.jumps_dagger,
        me.rates,
        me.efficiencies,
        me.iv,
        me.varmap,
        me.order,
    )
end

function split_equations(
    eqin::IndexedMeanfieldNoiseEquations,
)::Tuple{IndexedMeanfieldEquations,IndexedMeanfieldEquations}
    determ = IndexedMeanfieldEquations(
        eqin.equations,
        eqin.operator_equations,
        eqin.states,
        eqin.operators,
        eqin.hamiltonian,
        eqin.jumps,
        eqin.jumps_dagger,
        eqin.rates,
        eqin.iv,
        eqin.varmap,
        eqin.order,
    )
    noise = IndexedMeanfieldEquations(
        eqin.noise_equations,
        eqin.operator_noise_equations,
        eqin.states,
        eqin.operators,
        eqin.hamiltonian,
        eqin.jumps,
        eqin.jumps_dagger,
        eqin.efficiencies,
        eqin.iv,
        eqin.varmap,
        eqin.order,
    )
    return determ, noise
end

function merge_equations(
    determ::IndexedMeanfieldEquations,
    noise::IndexedMeanfieldEquations,
)::IndexedMeanfieldNoiseEquations
    return IndexedMeanfieldNoiseEquations(
        determ.equations,
        determ.operator_equations,
        noise.equations,
        noise.operator_equations,
        determ.states,
        determ.operators,
        determ.hamiltonian,
        determ.jumps,
        determ.jumps_dagger,
        determ.rates,
        noise.rates,
        determ.iv,
        determ.varmap,
        determ.order,
    )
end


function MTK.SDESystem(
    me::Union{MeanfieldNoiseEquations,IndexedMeanfieldNoiseEquations},
    iv = me.iv,
    vars = map(last, me.varmap),
    pars = nothing;
    complete_sys = true,
    kwargs...,
)
    determ, noise = split_equations(me)
    eqs = MTK.equations(determ)
    neqs = MTK.equations(noise)
    neqs_rhs = map(x->x.rhs, neqs)
    pars = isnothing(pars) ? extract_parameters(vcat(eqs..., neqs...), iv) : pars

    sys =  MTK.SDESystem(eqs, neqs_rhs, iv, vars, pars; kwargs...,)
    return complete_sys ? complete(sys) : sys
end

function MTK.System(de::Union{MeanfieldNoiseEquations,IndexedMeanfieldNoiseEquations};
    kwargs...,
)
    determ, noise = split_equations(de)
    return MTK.System(determ; kwargs...)
end

function Base.show(
    io::IO,
    de::Union{MeanfieldNoiseEquations,IndexedMeanfieldNoiseEquations},
)
    for i = 1:length(de.equations)
        write(io, "∂ₜ(")
        show(io, de.equations[i].lhs)
        write(io, ") = ")
        show(io, de.equations[i].rhs)
        write(io, " + dW(t)/dt[")
        show(io, de.noise_equations[i].rhs)
        write(io, "]")
        write(io, "\n")
    end
end
