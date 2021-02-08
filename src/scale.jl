function scale(he::DifferentialEquation{<:AbstractOperator, <:AbstractOperator})
    lhs_new = []
    rhs_new = []
    names = get_names(he)
    h = hilbert(he.lhs[1])
    cluster_Ns, cluster_orders, cluster_aons = get_cluster_stuff(h)
    he_avg = average(he)
    redundants = Average[]
    ref_avgs = Average[]
    dict_red = Dict()
    for it=1:length(he_avg) #dictionary and redundants from LHS
        if (he_avg.lhs[it] ∉ ref_avgs) && (he_avg.lhs[it] ∉ redundants)
            ref_avg, reds = get_ref_avg(he_avg.lhs[it], cluster_aons, names; no_adj=true)
            for red in reds
                dict_red[red] = ref_avg
            end
            all_reds = unique([reds; adjoint.(reds)])
            filter!(!isequal(ref_avg), all_reds)
            push!(lhs_new, ref_avg)
            push!(rhs_new, he_avg.rhs[it]) #apply dictionary later
        end
    end
    he_avg_new = DifferentialEquation(lhs_new, rhs_new, he.hamiltonian, he.jumps, he.rates)
    for it=1:length(he_avg_new) #dictionary from RHS
        avgs_rhs = unique(get_avgs(he_avg_new.rhs[it]))
        for avg_rhs in avgs_rhs
            if (avg_rhs ∉ ref_avgs) && (avg_rhs ∉ redundants)
                ref_avg, reds = get_ref_avg(avg_rhs, cluster_aons, names; no_adj=true)
                for red in reds
                    dict_red[red] = ref_avg
                end
            end
        end
    end
    he_avg_new = substitute(he_avg_new, dict_red)
    for it=1:length(cluster_aons)
        he_avg_new.rhs = set_scale_factors_rhs(he_avg_new, cluster_orders[it], cluster_aons[it], cluster_Ns[it])
    end
    he_avg_new = simplify_constants(he_avg_new)
    he_new = simplify_operators(unaverage(he_avg_new))
    he_new_scale = ScaleDifferentialEquation(he_new.lhs,he_new.rhs,he_new.hamiltonian,he_new.jumps,he_new.rates,Dict())
end

function set_scale_factors_rhs(he::DifferentialEquation, order, cluster_aon, N)
    rhs_new = eltype(he.rhs)[]
    for it=1:length(he)
        lhs = he.lhs[it]
        rhs = expand(copy(he.rhs[it]))
        aon_lhs_clus = intersect(acts_on(lhs), cluster_aon)
        len_alc = length(aon_lhs_clus)
        if order <= len_alc
            N_ = 1.0 #no scaling possible if all ops from cluster are on the LHS
            # TODO think about cluster*cluster
        else
            N_ = (N-len_alc) / (order-len_alc)
            # TODO # double sum -> look at RHS σ's: include r_c somehow (notes 19.01)
        end
        rhs_args = eltype(rhs.arguments)[]
        rhs_f = rhs.f # could be * (normal it is +)
        for rhs_arg in rhs.arguments
            aon_rhs_clus = intersect(acts_on(rhs_arg), cluster_aon)
            filter!(x->x∉aon_lhs_clus, aon_rhs_clus)
            len_arc = length(aon_rhs_clus)
            push!(rhs_args, rhs_arg*N_^len_arc)
        end
        push!(rhs_new, simplify_operators(rhs_f(rhs_args...)))
    end
    return rhs_new
end

function get_cluster_stuff(h::HilbertSpace)
    cluster_idx = findall(x->x isa ClusterSpace, h.spaces)
    cluster_Ns = getfield.(h.spaces[cluster_idx], :N)
    cluster_orders = getfield.(h.spaces[cluster_idx], :order)
    cluster_aons = []
    for i = 1:length(cluster_idx)
        aon_i = []
        for j=1:cluster_orders[i]
            push!(aon_i, ClusterAon(cluster_idx[i],j))
        end
        push!(cluster_aons, aon_i)
    end
    return cluster_Ns, cluster_orders, cluster_aons
end

"""
    scale(de::DifferentialEquation, identical_ops::Vector{AbstractOperator}, interaction_op::AbstractOperator, N::Number)

For a list of operators `identical_ops`, which represent idendical particles in the system, the set of equations `de` is adapted such
that it corresponds to a system with `N` identical particles. This means a factor of `N` (or a modification of it) is multiplied at
proper positions. The rules for the additional factors in the equations depend on the operator interacting with the identical operators `interaction_op`.
For a n-th order cumulant expansion at least n identical operators (acting on different Hilbert) spaces are needed.
"""

# Complete system skipping over unnecessary averages when scaling
function complete(de::ScaleDifferentialEquation{<:AbstractOperator,<:AbstractOperator};
            order=nothing, filter_func=nothing, mix_choice=maximum, kwargs...)

    names = get_names(de)
    J = de.jumps; H = de.hamiltonian; rates = de.rates
    h = hilbert(de.lhs[1])
    cluster_Ns, cluster_orders, cluster_aons = get_cluster_stuff(h)
    (order isa Nothing) ? (order_ = maximum(cluster_orders)) : order_ = order
    (maximum(order_) < maximum(get_order.(de.lhs))) && error("Cannot form cumulant expansion of derivative; you may want to use a higher order!")

    he_avg0 = average(de,order_)
    redundants = Average[]
    feed_redundants!(redundants,cluster_aons,he_avg0.lhs,names)
    missed = unique_ops(find_missing(he_avg0.rhs, he_avg0.lhs))
    filter!(x->isa(x,Average),missed)
    isnothing(filter_func) || filter!(filter_func, missed) # User-defined filter
    filter!(!in(redundants), missed)
    feed_redundants!(redundants,cluster_aons,missed,names)
    filter!(!in(redundants), missed)
    missed = unique_ops(missed)
    lhs_ = he_avg0.lhs
    rhs_ = he_avg0.rhs
    while !isempty(missed)
        ops = getfield.(missed, :operator)
        he = isempty(J) ? heisenberg(ops,H; kwargs...) : heisenberg(ops,H,J;rates=rates, kwargs...)
        he_avg = average(he,order_;mix_choice=mix_choice, kwargs...)
        rhs_ = [rhs_;he_avg.rhs]
        lhs_ = [lhs_;he_avg.lhs]
        missed = unique_ops(find_missing(rhs_,lhs_))
        filter!(x->isa(x,Average),missed)
        isnothing(filter_func) || filter!(filter_func, missed) # User-defined filter
        filter!(!in(redundants), missed)
        missed = unique_ops(missed)
        feed_redundants!(redundants,cluster_aons,missed,names)
        filter!(!in(lhs_),missed)
        filter!(!in(redundants), missed)
        missed = unique_ops(missed)
    end
    if !isnothing(filter_func)
        # Find missing values that are filtered by the custom filter function,
        # but still occur on the RHS; set those to 0
        missed = unique_ops(find_missing(rhs_, lhs_))
        filter!(x->isa(x,Average),missed)
        filter!(!filter_func, missed)
        subs = Dict(missed .=> 0)
        rhs_ = [substitute(r, subs) for r in rhs_]
    end
    dict = Dict()
    for it=1:length(lhs_)
        ref_avg, all_ids = get_ref_avg(lhs_[it], cluster_aons, names; no_adj=true)
        for id in all_ids
            dict[id] = ref_avg
        end
    end
    he_avg_scale = ScaleDifferentialEquation(lhs_, rhs_, H, J, rates, dict)
    return substitute(he_avg_scale, dict)
end

function complete(de::ScaleDifferentialEquation{<:Number,<:Number};order=nothing, filter_func=nothing, kwargs...)
    J = de.jumps; H = de.hamiltonian; rates = de.rates
    ops = getfield.(de.lhs, :operator)
    de0 = isempty(J) ? heisenberg(ops,H; kwargs...) : heisenberg(ops,H,J;rates=rates, kwargs...)
    return complete(de0;order=order,filter_func=filter_func,kwargs...)
end

### Auxiliary function
get_avgs(x) = Average[]
get_avgs(avg::Average) = [avg]
function get_avgs(t::NumberTerm)
    avgs = Average[]
    for arg in t.arguments
         append!(avgs, get_avgs(arg))
    end
    return unique(avgs)
end

function swap_aon(op::T, aon1, aon2, op_name_aon2::Symbol) where T<:BasicOperator
    acts_on(op) == aon1 || return op
    return T(op.hilbert, op_name_aon2, aon2)
end
function swap_aon(op::T, aon1, aon2, op_name_aon2::Symbol) where T<:Transition
    acts_on(op) == aon1 || return op
    return T(op.hilbert, op_name_aon2, op.i, op.j, aon2)
end
function swap_aon(op::BasicOperator, aon1::Vector, aon2::Vector, op_name_aon2::Vector)
    idx = findfirst(isequal(acts_on(op)), aon1)
    isnothing(idx) && return op
    return swap_aon(op, aon1[idx], aon2[idx], op_name_aon2[idx])
end
swap_aon(x::Average, aon1, aon2, op_name_aon2) = average(swap_aon(x.operator, aon1, aon2, op_name_aon2))
function swap_aon(op::OperatorTerm, aon1, aon2, op_name_aon2)
    args = op.arguments
    args_ = [swap_aon(arg, aon1, aon2, op_name_aon2) for arg in args]
    return simplify_operators((op.f)(args_...))
end
function swap_aon(op::OperatorTerm, aon1::Vector, aon2::Vector, op_name_aon2::Vector)
    args = [swap_aon(arg, aon1, aon2, op_name_aon2) for arg in op.arguments]
    return simplify_operators((op.f)(args...))
end

### all identical replace redundant averages
function intersect_aon(x, all_id_aon)
    aon_x = [acts_on(x)...]
    intersect!(aon_x, all_id_aon)
    return aon_x
end
function get_permuted_avgs(x, all_id_aon, names)
    aon_ls = intersect_aon(x, all_id_aon)
    avgs = Average[]
    for p in permutations(aon_ls)
        isequal(p, aon_ls) && continue
        push!(avgs, swap_aon(x, aon_ls, p, extract_names(names,p)))
    end
    return avgs
end

function get_names(de::Union{DifferentialEquation, ScaleDifferentialEquation})
    lhs_ops = unique(Iterators.flatten(get_operators.(de.lhs)))
    H_ops = unique(get_operators(de.hamiltonian))
    J_ops = unique(Iterators.flatten(get_operators.(de.jumps)))
    ops = vcat(lhs_ops, H_ops, J_ops)
    unique!(ops)
    ops = collect(Iterators.flatten(get_operators.(expand.(ops))))
    unique!(ops)
    filter!(x -> x isa BasicOperator, ops)
    sort!(ops; by=acts_on)
    aon_ls = acts_on.(ops)
    unique!(aon_ls)
    indices = [findfirst(x->acts_on(x)==aon_, ops) for aon_=aon_ls]
    ops = ops[indices]
    ops_names = getfield.(ops, :name)
    aons = acts_on.(ops)
    i_aons = [(isa(aon,Int) ? aon : aon.i) for aon in aons]
    ops_names_new = []
    for it=1:length(aons)
        aon = aons[it]
        i_aon = i_aons[it]
        if isa(aon, Int)
            idx = findfirst(x->x==i_aon, i_aons)
            push!(ops_names_new, ops_names[idx])
        elseif aon.j == 1 #ClusterAon
            idx = findall(x->x==i_aon, i_aons)
            push!(ops_names_new, ops_names[idx])
        end
    end
    return ops_names_new
end

### scale_complete() ###

function get_ref_avg(avg, all_id_aons, names; no_adj=false) #all_id_aons is a list of lists
    aon_ls = [intersect_aon(avg, all_id_aon) for all_id_aon in all_id_aons]
    if all(isempty.(aon_ls))
        return avg, []
    end
    len = length.(aon_ls)
    aon_ls_permus = [unique(sort.(permutations(all_id_aons[it], len[it]))) for it=1:length(all_id_aons)]
    ### get ref_avg, considering all "clusters"
    ref_avg = copy(avg)
    for it=1:length(all_id_aons)
        id_aon_ = all_id_aons[it][1:len[it]]
        name_idx = all_id_aons[it][1:len[it]]
        ref_avg = swap_aon(ref_avg, aon_ls[it], id_aon_, extract_names(names, name_idx))
    end
    aon_ls_ref = [intersect_aon(ref_avg, all_id_aon) for all_id_aon in all_id_aons]
    ### get all identical avgs (identical to ref_avg)
    all_ids = [ref_avg]
    for it=1:length(all_id_aons)
        if len[it] > 1
            avg_permus = []
            for it_a = 1:length(all_ids)
                permu_avgs = get_permuted_avgs(all_ids[it_a], all_id_aons[it], names)
                push!(avg_permus, permu_avgs...)
            end
            unique!(avg_permus)
            push!(all_ids, avg_permus...)
        end
    end
    for it=1:length(all_id_aons)
        all_ids_ = [swap_aon(all_ids[it_], aon_ls_ref[it], aon_p, extract_names(names, aon_p)) for it_ = 1:length(all_ids) for aon_p in aon_ls_permus[it]]
        push!(all_ids, all_ids_...)
    end
    filter!(!isequal(ref_avg), all_ids)
    unique!(all_ids)
    if no_adj # true -> to create dictionary
        return ref_avg, all_ids
    end
    all_ids = [all_ids; adjoint.(all_ids)]
    unique!(all_ids)
    filter!(!isequal(ref_avg), all_ids)
    return ref_avg, all_ids
end

function feed_redundants!(redundants::Vector, cluster_aons, avg_ls, names)
    for itm=1:length(avg_ls)
        ref_avg, id_avgs = get_ref_avg(avg_ls[itm], cluster_aons, names)
        if ref_avg ∉ redundants
            avg_ls[itm] = ref_avg
            id_avgs_adj = unique([id_avgs; adjoint.(id_avgs)])
            filter!(!isequal(ref_avg), id_avgs_adj)
            push!(redundants, id_avgs...)
            unique!(redundants)
        end
    end
    return redundants
end

unaverage(x) = x
unaverage(x::Average) = getfield(x, :operator)
unaverage(x::NumberTerm) = (x.f)(unaverage.(x.arguments)...)
function unaverage(de::DifferentialEquation{<:Number,<:Number})
    lhs_new = []; rhs_new = []
    for it=1:length(de)
        push!(lhs_new, unaverage(de.lhs[it]))
        push!(rhs_new, unaverage(de.rhs[it]))
    end
    return DifferentialEquation(lhs_new, rhs_new, de.hamiltonian, de.jumps, de.rates)
end
