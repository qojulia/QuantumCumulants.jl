function scale(he::DifferentialEquation{<:AbstractOperator, <:AbstractOperator})
    lhs_new = eltype(he.lhs)[]
    rhs_new = eltype(he.rhs)[]
    names = get_names(he)

    h = hilbert(he.lhs[1])
    cluster_Ns, cluster_orders, cluster_aons = get_cluster_stuff(h)
    he_avg = average(he)
    redundants = Average[]
    ref_avgs = Average[]
    dict_red = Dict()
    for it=1:length(he_avg)
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
    he_avg_new = substitute(he_avg_new, dict_red)
    for it=1:length(cluster_aons)
        he_avg_new.rhs = set_scale_factors_rhs(he_avg_new, cluster_orders[it], cluster_aons[it], cluster_Ns[it])
    end
    he_new = unaverage(he_avg_new)
    # TODO create ScaleDifferentialEquation from he_new
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
            N_ = (N-l_c) / (order-l_c)
            # TODO # double sum -> look at RHS σ's: include r_c somehow (notes 19.01)
        end
        rhs_args = eltype(rhs.arguments)[]
        rhs_f = rhs.f # could be * (normal it is +)
        for rhs_arg in rhs.arguments
            aon_rhs_clus = intersect(acts_on(rhs_arg), cluster_aon)
            filter!(x->aon_lhs_clus, aon_rhs_clus)
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
        for j=1:c_order[i]
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
function scale(de::AbstractEquation{<:Number,<:Number}, cluster_aons, interaction_aon, N)
    names = get_names(de)
    # TODO: clean up dispatch
    # TODO: don't use variable names identical to function arguments!!
    !isa(N, Vector) && (cluster_aons = [cluster_aons]; interaction_aon = [interaction_aon]; N = [N]) #convert to vector
    if isa(interaction_aon, Vector{<:AbstractOperator}) #convert operator to acts_on-number
        cluster_aons = [unique(acts_on.(id_aons)) for id_aons in cluster_aons]
        interaction_aon = acts_on.(interaction_aon)
    end
    cluster_aons = sort.(cluster_aons)
    de_s = scale_nosub(de, cluster_aons, interaction_aon, N, names)
    sub_ref_avg(de_s, cluster_aons, interaction_aon, N, names)
end

function scale(de::AbstractEquation{<:AbstractOperator,<:AbstractOperator})
    h = hilbert(de.lhs[1])
    # (h isa ClusterSpace) && error("Not yet implemented!") # TODO: implement single cluster interacting with itself
    cluster_idx = findall(x->x isa ClusterSpace, h.spaces)
    c_N = getfield.(h.spaces[cluster_idx], :N)
    c_order = getfield.(h.spaces[cluster_idx], :order)
    cluster_aons = []
    for i = 1:length(cluster_idx)
        aon_i = []
        for j=1:c_order[i]
            push!(aon_i, ClusterAon(cluster_idx[i],j))
        end
        push!(cluster_aons, aon_i)
    end
    H = simplify_operators(de.hamiltonian)

    cluster_aons_new = []
    N_new = Number[]
    for i=1:length(cluster_aons)
        c_aon_ = cluster_aons[i][1]
        for inter_ ∈ interacting_aons_
            if c_aon_ ∈ inter_
                idx = findall(x->isa(x, Int), inter_)
                # length(idx)==1 || error("Also not implemented!") # TODO
                push!(cluster_aons_new, cluster_aons[i])
                push!(N_new, c_N[i])
            end
        end
    end
    interaction_aon = collect(Iterators.flatten(interacting_aons))
    return scale(de, cluster_aons_new, interaction_aon, N_new)
end

function _interacting_aons(t::OperatorTerm{<:typeof(+)})
    aon = []
    for arg in t.arguments
        push!(aon, _interacting_aons(arg))
    end
    unique!(aon)
    sort!.(aon)
    function _f(x)
        length(x)==2 || return false
        if x[1] isa Int
            return x[2] isa ClusterAon
        elseif x[2] isa Int
            return x[1] isa ClusterAon
        else # TODO: two interacting clusters
            return false
        end
    end
    filter!(_f, aon)
    return aon
end
_interacting_aons(t::OperatorTerm{<:typeof(*)}) = acts_on(t)
_interacting_aons(x) = []

function scale_nosub(de::AbstractEquation{<:Number,<:Number}, cluster_aons::Vector, interaction_aon::Vector, N::Vector{<:Number}, names)
    @assert length(interaction_aon) == length(N) == length(cluster_aons)
    de_scale = de
    for it=1:length(N)
        de_scale = scale_nosub(get_DiffEq_from_scale(de_scale), cluster_aons[it], interaction_aon[it], N[it], names)
    end
    return de_scale
end
function scale_nosub(de::DifferentialEquation, cluster_aons_::Vector, interaction_aon::Int, N::Number, names)
    cluster_aons = sort(cluster_aons_)
    idx = [filter_identical(x, cluster_aons) for x in de.lhs]
    de_lhs = de.lhs[idx]
    de_rhs = de.rhs[idx]
    for it=1:length(de_lhs)
        aon_lhs = [acts_on(de_lhs[it])...]
        aon_lhs_id = intersect(cluster_aons, aon_lhs)
        if (interaction_aon ∉ aon_lhs) || (length(cluster_aons) == length(aon_lhs_id))
            N_sc = 1
        else
            N_sc = (N-length(aon_lhs_id))/(length(cluster_aons) - length(aon_lhs_id)) #double sum -> extension to RHS??
        end
        tmp_dict = Dict{Average,Number}()
        avgs_rhs = get_avgs(de_rhs[it])
        for avg_rhs in avgs_rhs
            avg_aon_rhs = [acts_on(avg_rhs)...]
            intersect!(avg_aon_rhs, cluster_aons)
            if !isempty(avg_aon_rhs)
                new_avg = copy(avg_rhs)
                for it_aon=1:length(avg_aon_rhs)
                    name_idx = cluster_aons[it_aon].i + cluster_aons[it_aon].j - 1
                    new_avg = swap_aon(new_avg, avg_aon_rhs[it_aon], cluster_aons[it_aon], names[name_idx])
                end
                if length(avg_aon_rhs) > length(aon_lhs_id)
                    tmp_dict[avg_rhs] = N_sc*new_avg
                else
                    tmp_dict[avg_rhs] = new_avg
                end
            end
        end
        de_rhs[it] = substitute(de_rhs[it], tmp_dict)
    end
    he = DifferentialEquation(de_lhs, de_rhs, de.hamiltonian, de.jumps, de.rates)
    return he
end

function sub_ref_avg(de_s::AbstractEquation{<:Number,<:Number}, cluster_aons, interaction_aon, N, names)
    dict = Dict()
    for it=1:length(de_s.lhs)
        ref_avg, all_ids = get_ref_avg(de_s.lhs[it], cluster_aons, names; no_adj=true)
        for id in all_ids
            dict[id] = ref_avg
        end
    end
    he_s = substitute(de_s, dict)
    return ScaleDifferentialEquation(he_s.lhs,he_s.rhs,he_s.hamiltonian,he_s.jumps,he_s.rates,N,cluster_aons,interaction_aon,dict)
end

# Complete system skipping over unnecessary averages when scaling
function scale_complete(rhs::Vector{<:Number}, vs::Vector{<:Number}, H, J, rates,
            cluster_aons::Vector, interaction_aons::Vector, N::Vector, names;
            order=nothing, filter_func=nothing, mix_choice=maximum,
            kwargs...)

    order_lhs = maximum(get_order.(vs))
    order_rhs = maximum(get_order.(rhs))
    if order isa Nothing
        order_ = max(order_lhs, order_rhs)
    else
        order_ = order
    end
    maximum(order_) >= order_lhs || error("Cannot form cumulant expansion of derivative; you may want to use a higher order!")

    lhs_init_ops = getfield.(vs, :operator)
    de_ops_init = average(heisenberg(lhs_init_ops, H, J; rates=rates, kwargs...),order_)
    vs_ = copy(de_ops_init.lhs)
    redundants = Average[] #create identical specific redundants later
    feed_redundants!(redundants,cluster_aons,vs_,names)
    rhs_ = [cumulant_expansion(r, order_) for r in de_ops_init.rhs]
    missed = unique_ops(find_missing(rhs_, vs_))
    filter!(x->isa(x,Average),missed)
    isnothing(filter_func) || filter!(filter_func, missed) # User-defined filter
    filter!(!in(redundants), missed)
    ### redundant
    feed_redundants!(redundants,cluster_aons,missed,names)
    filter!(!in(redundants), missed)
    missed = unique_ops(missed)

    while !isempty(missed)
        ops = getfield.(missed, :operator)
        he = isempty(J) ? heisenberg(ops,H; kwargs...) : heisenberg(ops,H,J;rates=rates, kwargs...)
        he_avg = average(he,order_;mix_choice=mix_choice, kwargs...)
        rhs_ = [rhs_;he_avg.rhs]
        vs_ = [vs_;he_avg.lhs]
        missed = unique_ops(find_missing(rhs_,vs_))
        filter!(x->isa(x,Average),missed)
        isnothing(filter_func) || filter!(filter_func, missed) # User-defined filter
        filter!(!in(redundants), missed)
        missed = unique_ops(missed)
        feed_redundants!(redundants,cluster_aons,missed,names)
        filter!(!in(vs_),missed)
        filter!(!in(redundants), missed)
        missed = unique_ops(missed)
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

    dict = Dict()
    for it=1:length(vs_)
        ref_avg, all_ids = get_ref_avg(vs_[it], cluster_aons, names; no_adj=true)
        for id in all_ids
            dict[id] = ref_avg
        end
    end
    return ScaleDifferentialEquation(vs_, rhs_, H, J, rates, N, cluster_aons, interaction_aons, dict)
end
function complete(de::ScaleDifferentialEquation{<:Number,<:Number};kwargs...)
    names = get_names(de)
    return scale_complete(de.rhs, de.lhs, de.hamiltonian, de.jumps, de.rates, de.identicals, de.interactions, de.factors, names; kwargs...)
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

function filter_identical(avg::Average, cluster_aons::Vector)
    avg_aon = acts_on(avg)
    avg_cluster_aons = intersect(cluster_aons, avg_aon)
    if isempty(avg_cluster_aons)
        return true
    else
        it_max = findfirst(isequal(maximum(avg_cluster_aons)), cluster_aons)
        return isequal(avg_cluster_aons, cluster_aons[1:it_max])
    end
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
        name_idx = getfield.(p, :i) .+ getfield.(p, :j) .- 1
        push!(avgs, swap_aon(x, aon_ls, p, names[name_idx]))
    end
    return avgs
end
function get_redundant(lhs, all_id_aon, names) #only for averages with more than one all_id_aon
    lhsf = filter(x->length(intersect_aon(x, all_id_aon))>1, lhs)
    lhs_redundants = []
    dict_lhs_redundant = Dict{Average,Average}()
    for itl=1:length(lhsf)
        lhsf1 = lhsf[itl]
        if lhsf1 ∉ lhs_redundants #if it is in, it will be deleted later
            id_avgs = get_permuted_avgs(lhsf1, all_id_aon, names)
            id_avgs_adj = adjoint.(id_avgs)
            filter!(!isequal(lhsf1), id_avgs)
            filter!(!isequal(lhsf1), id_avgs_adj)
            append!(lhs_redundants, id_avgs)
            append!(lhs_redundants, id_avgs_adj)
            for id_avg in id_avgs
                dict_lhs_redundant[id_avg] = lhsf1
            end
        end
    end
    return lhs_redundants, dict_lhs_redundant
end
function filter_redundant(de::DifferentialEquation, cluster_aons, names)
    lhs_ = de.lhs; rhs_ = de.rhs
    lhsf = eltype(lhs_)[]; rhsf = eltype(rhs_)[]
    reds, dict_reds = get_redundant(lhs_, cluster_aons, names)
    for it=1:length(lhs_)
        lhs_[it] ∈ reds || (push!(lhsf, lhs_[it]); push!(rhsf, substitute(rhs_[it], dict_reds)))
    end
    return DifferentialEquation(lhsf, rhsf, de.hamiltonian, de.jumps, de.rates), dict_reds
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
    i_aons = acts_on.(ops); i_aons = [(isa(i_aon,Int) ? i_aon : i_aon.i) for i_aon in i_aons]
    ops_names_new = []
    for it=1:i_aons[end]
        idx = findall(x->x==it, i_aons); (length(idx)==1 && (idx = idx[1]))
        push!(ops_names_new, ops_names[idx])
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
        name_idx = getfield.(id_aon_, :i) .+ getfield.(id_aon_, :j) .- 1
        ref_avg = swap_aon(ref_avg, aon_ls[it], id_aon_, names[name_idx])
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
        all_ids_ = [swap_aon(all_ids[it_], aon_ls_ref[it], aon_p, names[getfield.(aon_p, :i) .+ getfield.(aon_p, :j) .- 1]) for it_ = 1:length(all_ids) for aon_p in aon_ls_permus[it]]
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
function unaverage(x::NumberTerm)
    return simplify_operators((x.f)(unaverage.(x.arguments)...))
end
