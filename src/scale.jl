
### scaling ###
get_avgs(x) = []
get_avgs(avg::Average) = [avg]
function get_avgs(t::NumberTerm)
    avgs = []
    for arg in t.arguments
         append!(avgs, get_avgs(arg))
    end
    return unique(avgs)
end

"""
    scale(de::DifferentialEquation, identical_ops::Vector{AbstractOperator}, interaction_op::AbstractOperator, N::Number)

For a list of operators `identical_ops`, which represent idendical particles in the system, the set of equations `de` is adapted such
that it corresponds to a system with `N` identical particles. This means a factor of `N` (or a modification of it) is multiplied at
proper positions. The rules for the additional factors in the equations depend on the operator interacting with the identical operators `interaction_op`.
For a n-th order cumulant expansion at least n identical operators (acting on different Hilbert) spaces are needed.
"""
function scale(de::DifferentialEquation, identical_ops::Vector{<:AbstractOperator}, interaction_op::AbstractOperator, N::Number)
    identical_aons = unique(acts_on.(identical_ops))
    interaction_aon = acts_on(interaction_op)
    scale(de::DifferentialEquation, identical_aons, interaction_aon, N)
end
function scale(de::Union{ScaleDifferentialEquation, DifferentialEquation}, identical_aons::Vector, interaction_aon::Vector, N::Vector{<:Number})
    @assert length(interaction_aon) == length(N) == length(identical_aons)
    de_scale = de
    for it=1:length(N)
        de_scale = scale(de_scale, identical_aons[it], interaction_aon[it], N[it])
    end
    return de_scale
end
function scale(de::ScaleDifferentialEquation, identical_aons::Vector, interaction_aon, N::Number)
    he = scale(get_DiffEq_from_scale(de), identical_aons, interaction_aon, N)
    he_scale = ScaleDifferentialEquation(he.lhs,he.rhs,he.hamiltonian,he.jumps,he.rates,[de.factors;he.factors],[de.identicals;he.identicals],[de.interactions;he.interactions],[de.dictionaries; he.dictionaries])
end
function scale(de::DifferentialEquation, identical_aons::Vector{Int}, interaction_aon::Int, N::Number)
    names = get_names(de)
    sort!(identical_aons)
    de_lhs = Number[]
    de_rhs = Number[]
    for it=1:length(de)
        all_identical_filter(de.lhs[it], identical_aons) && (push!(de_lhs, de.lhs[it]); push!(de_rhs, de.rhs[it]))
    end
    for it=1:length(de_lhs)
        avg_aon_lhs = [acts_on(de_lhs[it])...]
        avg_aon_lhs_all_id = filter(x->(x∈identical_aons), avg_aon_lhs)
        if (interaction_aon ∉ avg_aon_lhs) || (length(identical_aons) == length(avg_aon_lhs_all_id))
            N_sc = 1
        else
            N_sc = (N-length(avg_aon_lhs_all_id))/(length(identical_aons) - length(avg_aon_lhs_all_id)) #double sum -> extension to RHS??
        end
        tmp_dict = Dict()
        avgs_rhs = get_avgs(de_rhs[it])
        for avg_rhs in avgs_rhs
            avg_aon_rhs = [acts_on(avg_rhs)...]
            avg_aon_rhs_all_id = filter(x->x∈identical_aons, avg_aon_rhs)
            if !isempty(avg_aon_rhs_all_id)
                it_max = findfirst(x->isequal(x,maximum(avg_aon_rhs_all_id)), identical_aons)
                new_avg = copy(avg_rhs)
                for it_aon=1:length(avg_aon_rhs_all_id)
                    new_avg = swap_aon(new_avg, avg_aon_rhs_all_id[it_aon], identical_aons[it_aon], names[identical_aons[it_aon]])
                end
                if length(avg_aon_rhs_all_id) > length(avg_aon_lhs_all_id)
                    tmp_dict[avg_rhs] = N_sc*new_avg
                else
                    tmp_dict[avg_rhs] = new_avg
                end
            end
        end
        de_rhs[it] = substitute(de_rhs[it], tmp_dict)
    end
    he = DifferentialEquation(de_lhs, de_rhs, de.hamiltonian, de.jumps, de.rates)
    he, scale_dict =  filter_redundant_all_id_aon(he, identical_aons, names)
    he_scale = ScaleDifferentialEquation(he.lhs,he.rhs,he.hamiltonian,he.jumps,he.rates,[N],[identical_aons],[interaction_aon],[scale_dict])
end

function all_identical_filter(avg::Average, identical_aons::Vector{Int})
    avg_aon = [acts_on(avg)...]
    avg_identical_aons = filter(x->x∈identical_aons, avg_aon)
    if isempty(avg_identical_aons)
        return true
    else
        it_max = findfirst(x->isequal(x,maximum(avg_identical_aons)), identical_aons)
        return isequal(avg_identical_aons, identical_aons[1:it_max])
    end
end
function swap_aon(op::Union{Destroy, Create, Transition}, aon1, aon2, op_name_aon2)
    if acts_on(op) != aon1
        return op
    elseif isa(op, Destroy)
        return Destroy(op.hilbert, op_name_aon2, aon2)
    elseif isa(op, Create)
        return Create(op.hilbert, op_name_aon2, aon2)
    elseif isa(op, Transition)
        return Transition(op.hilbert, op_name_aon2, op.i, op.j, aon2)
    end
end
function swap_aon(op::Union{Destroy, Create, Transition}, aon1::Vector, aon2::Vector, op_name_aon2::Vector)
    @assert length(aon1) == length(aon2) == length(op_name_aon2)
    if acts_on(op) ∉ aon1
        return op
    else
        it_aon = findfirst(x->x==acts_on(op), aon1)
        if isa(op, Destroy)
            return Destroy(op.hilbert, op_name_aon2[it_aon], aon2[it_aon])
        elseif isa(op, Create)
            return Create(op.hilbert, op_name_aon2[it_aon], aon2[it_aon])
        elseif isa(op, Transition)
            return Transition(op.hilbert, op_name_aon2[it_aon], op.i, op.j, aon2[it_aon])
        end
    end
end
swap_aon(x::Average, aon1, aon2, op_name_aon2) = average(swap_aon(x.operator, aon1, aon2, op_name_aon2))
function swap_aon(op::OperatorTerm, aon1, aon2, op_name_aon2)
    args = op.arguments
    args_ = [swap_aon(arg, aon1, aon2, op_name_aon2) for arg in args]
    return simplify_operators((op.f)(args_...))
end
function swap_aon(op::OperatorTerm, aon1::Vector, aon2::Vector, op_name_aon2::Vector)
    @assert length(aon1) == length(aon2) == length(op_name_aon2)
    args = op.arguments
    args_ = [swap_aon(arg, aon1, aon2, op_name_aon2) for arg in args]
    return simplify_operators((op.f)(args_...))
end

### all identical replace reduandant averages
function aon_in_all_id_aon(x, all_id_aon)
    aon_x = [acts_on(x)...]
    filter!(x->x∈all_id_aon, aon_x)
end
function get_all_id_permuted_avgs(x, all_id_aon, names)
    aon_ls = aon_in_all_id_aon(x, all_id_aon)
    aon_ls_permus = permutations(aon_ls)
    [swap_aon(x, aon_ls, aon_p, names[aon_p]) for aon_p in aon_ls_permus]
end
function get_redundant_all_id_aon(lhs, all_id_aon, names) #only for averages with more then one all_id_aon
    lhsf = filter(x->length(aon_in_all_id_aon(x, all_id_aon))>1, lhs)
    lhs_redundants = []
    dict_lhs_redundant = Dict()
    len_lhsf = length(lhsf)
    for itl=1:length(lhsf)
        lhsf1 = lhsf[itl]
        if lhsf1 ∉ lhs_redundants #if it is in, it will be deleted later
            id_avgs = get_all_id_permuted_avgs(lhsf1, all_id_aon, names)[2:end] #TODO
            id_avgs_adj = adjoint.(id_avgs)
            filter!(x->x≠lhsf1, id_avgs)
            filter!(x->x≠lhsf1, id_avgs_adj)
            push!(lhs_redundants, [id_avgs; id_avgs_adj]...)
            [dict_lhs_redundant[id_avg] = lhsf1 for id_avg in id_avgs]
            [dict_lhs_redundant[id_avg'] = lhsf1' for id_avg in id_avgs] #maybe not needed TODO
        end
    end
    return lhs_redundants, dict_lhs_redundant
end
function filter_redundant_all_id_aon(de::DifferentialEquation, all_id_aon, names)
    lhs_ = de.lhs; rhs_ = de.rhs
    lhsf = eltype(lhs_)[]; rhsf = eltype(rhs_)[]
    reds, dict_reds = get_redundant_all_id_aon(lhs_, all_id_aon, names)
    for it=1:length(lhs_)
        lhs_[it] ∉ reds && (push!(lhsf, lhs_[it]); push!(rhsf, substitute(rhs_[it], dict_reds)))
    end
    return DifferentialEquation(lhsf, rhsf, de.hamiltonian, de.jumps, de.rates), dict_reds
end

function get_names(de::Union{DifferentialEquation, ScaleDifferentialEquation})
    lhs_ops = unique(collect(Iterators.flatten(get_operators.(de.lhs))))
    H_ops = unique(get_operators(de.hamiltonian))
    J_ops = unique(collect(Iterators.flatten(get_operators.(de.jumps))))
    ops = unique(collect(Iterators.flatten(get_operators.([lhs_ops; H_ops; J_ops]))))
    sort!(ops; by=acts_on)
    indices = [findfirst(x->acts_on(x)==aon_, ops) for aon_=acts_on(ops[1]):acts_on(ops[end])]
    ops = ops[indices]
    names = [op.name for op in ops]
end


### scale_complete() ###

function get_ref_avg_and_all_identical_avgs_adj(avg, all_id_aon, names)
    aon_ls = aon_in_all_id_aon(avg, all_id_aon)
    if isempty(aon_ls)
        return avg, []
    end
    len = length(aon_ls)
    aon_ls_permus = unique(sort.(permutations(all_id_aon, len)))
    ref_avg = swap_aon(avg, aon_ls, all_id_aon[1:len], names[all_id_aon[1:len]])
    aon_ls_ref = aon_in_all_id_aon(ref_avg, all_id_aon)
    all_ids =  [swap_aon(ref_avg, aon_ls_ref, aon_p, names[aon_p]) for aon_p in aon_ls_permus]
    if len > 1
        avg_permus = get_all_id_permuted_avgs(ref_avg, all_id_aon, names)
        filter!(x->x≠ref_avg, avg_permus)
        for it=1:length(avg_permus)
            all_ids_ = [swap_aon(avg_permus[it], aon_ls_ref, aon_p, names[aon_p]) for aon_p in aon_ls_permus]
            push!(all_ids, all_ids_...)
        end
    end
    all_ids = [all_ids; adjoint.(all_ids)]
    unique!(all_ids)
    filter!(x->x≠ref_avg, all_ids)
    return ref_avg, all_ids
end

function complete(de::ScaleDifferentialEquation{<:Number,<:Number};kwargs...)
    names = get_names(de)
    return scale_complete(de.rhs, de.lhs, de.hamiltonian, de.jumps, de.rates, de.identicals, de.interactions, de.factors, names; kwargs...)
end

function scale_complete(rhs::Vector{<:Number}, vs::Vector{<:Number}, H, J, rates, identical_aons::Vector, interaction_aons::Vector, N::Vector, names; order=nothing, filter_func=nothing, mix_choice=maximum, kwargs...)
    order_lhs = maximum(get_order.(vs))
    order_rhs = maximum(get_order.(rhs))
    if order isa Nothing
        order_ = max(order_lhs, order_rhs)
    else
        order_ = order
    end
    maximum(order_) >= order_lhs || error("Cannot form cumulant expansion of derivative; you may want to use a higher order!")

    lhs_init_ops = getfield.(vs, :operator)
    de_ops_init = average(heisenberg(lhs_init_ops, H, J; rates=rates, kwargs...),order)
    vs_ = copy(de_ops_init.lhs)
    redundants = [] #create identical specific redundants later
    ref_avgs = []
    feed_redundants(redundants,identical_aons,vs_,names)
    rhs_ = [cumulant_expansion(r, order_) for r in de_ops_init.rhs]
    missed = unique_ops(find_missing(rhs_, vs_))
    filter!(x->isa(x,Average),missed)
    isnothing(filter_func) || filter!(filter_func, missed) # User-defined filter
    filter!(x->x∉redundants, missed)
    ### redundant
    feed_redundants(redundants,identical_aons,missed,names)
    filter!(x->x∉redundants, missed)

    while !isempty(missed)
        ops = getfield.(missed, :operator)
        he = isempty(J) ? heisenberg(ops,H; kwargs...) : heisenberg(ops,H,J;rates=rates, kwargs...)
        he_avg = average(he,order_;mix_choice=mix_choice, kwargs...)
        rhs_ = [rhs_;he_avg.rhs]
        vs_ = [vs_;he_avg.lhs]
        missed = unique_ops(find_missing(rhs_,vs_))
        filter!(x->isa(x,Average),missed)
        isnothing(filter_func) || filter!(filter_func, missed) # User-defined filter
        filter!(x->x∉redundants, missed)
        feed_redundants(redundants,identical_aons,missed,names)
        filter!(x->x∉redundants, missed)
    end
    if !isnothing(filter_func)
        # Find missing values that are filtered by the custom filter function,
        # but still occur on the RHS; set those to 0
        missed = unique_ops(find_missing(rhs_, vs_))
        filter!(x->isa(x,Average),missed)
        filter!(!filter_func, missed)
        subs = Dict(missed .=> 0)
        rhs_ = [substitute(r, subs) for r in rhs_]
    end
    he = DifferentialEquation(vs_, rhs_, H, J, rates)
    he_scale = scale(he, identical_aons, interaction_aons, N)
    return he_scale
end

function feed_redundants(redundants::Vector, identical_aons, avg_ls, names)
    for it=1:length(identical_aons)
        for itm=1:length(avg_ls)
            ref_avg, id_avgs = get_ref_avg_and_all_identical_avgs_adj(avg_ls[itm], identical_aons[it], names)
            if ref_avg ∉ redundants
                id_avgs_adj = unique([id_avgs; adjoint.(id_avgs)])
                filter!(x->x≠ref_avg, id_avgs_adj)
                push!(redundants, id_avgs...)
                unique!(redundants)
            end
        end
    end
end
