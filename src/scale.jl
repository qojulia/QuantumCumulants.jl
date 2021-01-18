"""
    scale(de::DifferentialEquation, identical_ops::Vector{AbstractOperator}, interaction_op::AbstractOperator, N::Number)

For a list of operators `identical_ops`, which represent idendical particles in the system, the set of equations `de` is adapted such
that it corresponds to a system with `N` identical particles. This means a factor of `N` (or a modification of it) is multiplied at
proper positions. The rules for the additional factors in the equations depend on the operator interacting with the identical operators `interaction_op`.
For a n-th order cumulant expansion at least n identical operators (acting on different Hilbert) spaces are needed.
"""
function scale(de::AbstractEquation{<:Number,<:Number}, identical_aons, interaction_aon, N)
    names = get_names(de)
    # TODO: clean up dispatch
    # TODO: don't use variable names identical to function arguments!!
    !isa(N, Vector) && (identical_aons = [identical_aons]; interaction_aon = [interaction_aon]; N = [N]) #convert to vector
    if isa(interaction_aon, Vector{<:AbstractOperator}) #convert operator to acts_on-number
        identical_aons = [unique(acts_on.(id_aons)) for id_aons in identical_aons]
        interaction_aon = acts_on.(interaction_aon)
    end
    identical_aons = sort.(identical_aons)
    de_s = scale_nosub(de, identical_aons, interaction_aon, N, names)
    sub_ref_avg(de_s, identical_aons, interaction_aon, N, names)
end

function scale(de::AbstractEquation{<:Number,<:Number})
    h = hilbert(de.lhs[1])
    (h isa ClusterSpace) && error("Not yet implemented!") # TODO: implement single cluster interacting with itself
    c_idx = findall(x->x isa ClusterSpace, h.spaces)
    N = getfield.(h.spaces[c_idx], :N)
    c_order = getfield.(h.spaces[c_idx], :order)
    identical_aons = []
    for i = 1:length(c_idx)
        aon_i = []
        for j=1:c_order[i]
            push!(aon_i, ClusterAon(c_idx[i],j))
        end
        push!(identical_aons, aon_i)
    end

    H = simplify_operators(de.hamiltonian)
    interacting_aons_ = _interacting_aons(H)

    interacting_aons = []
    identical_aons_new = []
    N_new = Number[]
    for i=1:length(identical_aons)
        c_aon_ = identical_aons[i][1]
        for inter_ ∈ interacting_aons_
            if c_aon_ ∈ inter_
                idx = findall(x->isa(x, Int), inter_)
                length(idx)==1 || error("Also not implemented!") # TODO
                push!(interacting_aons, inter_[idx])
                push!(identical_aons_new, identical_aons[i])
                push!(N_new, N[i])
            end
        end
    end
    interaction_aon = collect(Iterators.flatten(interacting_aons))
    return scale(de, identical_aons_new, interaction_aon, N_new)
end

function _interacting_aons(t::OperatorTerm{<:typeof(+)})
    aon = []
    for arg in t.arguments
        push!(aon, _interacting_aons(arg))
    end
    unique!(aon)
    sort!.(aon)
    return aon
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

function scale_nosub(de::AbstractEquation{<:Number,<:Number}, identical_aons::Vector, interaction_aon::Vector, N::Vector{<:Number}, names)
    @assert length(interaction_aon) == length(N) == length(identical_aons)
    de_scale = de
    for it=1:length(N)
        de_scale = scale_nosub(get_DiffEq_from_scale(de_scale), identical_aons[it], interaction_aon[it], N[it], names)
    end
    return de_scale
end
function scale_nosub(de::DifferentialEquation, identical_aons_::Vector, interaction_aon::Int, N::Number, names)
    identical_aons = sort(identical_aons_)
    idx = [filter_identical(x, identical_aons) for x in de.lhs]
    de_lhs = de.lhs[idx]
    de_rhs = de.rhs[idx]
    for it=1:length(de_lhs)
        aon_lhs = [acts_on(de_lhs[it])...]
        aon_lhs_id = intersect(identical_aons, aon_lhs)
        if (interaction_aon ∉ aon_lhs) || (length(identical_aons) == length(aon_lhs_id))
            N_sc = 1
        else
            N_sc = (N-length(aon_lhs_id))/(length(identical_aons) - length(aon_lhs_id)) #double sum -> extension to RHS??
        end
        tmp_dict = Dict{Average,Number}()
        avgs_rhs = get_avgs(de_rhs[it])
        for avg_rhs in avgs_rhs
            avg_aon_rhs = [acts_on(avg_rhs)...]
            intersect!(avg_aon_rhs, identical_aons)
            if !isempty(avg_aon_rhs)
                new_avg = copy(avg_rhs)
                for it_aon=1:length(avg_aon_rhs)
                    name_idx = identical_aons[it_aon].i + identical_aons[it_aon].j - 1
                    new_avg = swap_aon(new_avg, avg_aon_rhs[it_aon], identical_aons[it_aon], names[name_idx])
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

function sub_ref_avg(de_s::AbstractEquation{<:Number,<:Number}, identical_aons, interaction_aon, N, names)
    dict = Dict()
    for it=1:length(de_s.lhs)
        ref_avg, all_ids = get_ref_avg(de_s.lhs[it], identical_aons, names; no_adj=true)
        for id in all_ids
            dict[id] = ref_avg
        end
    end
    he_s = substitute(de_s, dict)
    return ScaleDifferentialEquation(he_s.lhs,he_s.rhs,he_s.hamiltonian,he_s.jumps,he_s.rates,N,identical_aons,interaction_aon,dict)
end

# Complete system skipping over unnecessary averages when scaling
function scale_complete(rhs::Vector{<:Number}, vs::Vector{<:Number}, H, J, rates,
            identical_aons::Vector, interaction_aons::Vector, N::Vector, names;
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
    feed_redundants!(redundants,identical_aons,vs_,names)
    rhs_ = [cumulant_expansion(r, order_) for r in de_ops_init.rhs]
    missed = unique_ops(find_missing(rhs_, vs_))
    filter!(x->isa(x,Average),missed)
    isnothing(filter_func) || filter!(filter_func, missed) # User-defined filter
    filter!(!in(redundants), missed)
    ### redundant
    feed_redundants!(redundants,identical_aons,missed,names)
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
        feed_redundants!(redundants,identical_aons,missed,names)
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
        ref_avg, all_ids = get_ref_avg(vs_[it], identical_aons, names; no_adj=true)
        for id in all_ids
            dict[id] = ref_avg
        end
    end
    return ScaleDifferentialEquation(vs_, rhs_, H, J, rates, N, identical_aons, interaction_aons, dict)
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

function filter_identical(avg::Average, identical_aons::Vector)
    avg_aon = acts_on(avg)
    avg_identical_aons = intersect(identical_aons, avg_aon)
    if isempty(avg_identical_aons)
        return true
    else
        it_max = findfirst(isequal(maximum(avg_identical_aons)), identical_aons)
        return isequal(avg_identical_aons, identical_aons[1:it_max])
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
function filter_redundant(de::DifferentialEquation, identical_aons, names)
    lhs_ = de.lhs; rhs_ = de.rhs
    lhsf = eltype(lhs_)[]; rhsf = eltype(rhs_)[]
    reds, dict_reds = get_redundant(lhs_, identical_aons, names)
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
    filter!(x -> x isa BasicOperator, ops)
    sort!(ops; by=acts_on)
    aon_ls = map(acts_on, ops)
    unique!(aon_ls)
    indices = [findfirst(x->acts_on(x)==aon_, ops) for aon_=aon_ls]
    ops = ops[indices]
    return getfield.(ops, :name)
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

function feed_redundants!(redundants::Vector, identical_aons, avg_ls, names)
    for itm=1:length(avg_ls)
        ref_avg, id_avgs = get_ref_avg(avg_ls[itm], identical_aons, names)
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
