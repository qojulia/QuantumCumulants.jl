function scale(he::HeisenbergEquation, scale_aons, N; simplify=true, kwargs...)
    M = length(scale_aons)
    rhs = QNumber[]
    names = get_names(he)
    he_avg = average(he) # Ensure correct form of expressions required for scaling

    reference_terms = []
    append!(reference_terms, filter(x->!isempty(intersect(acts_on(x), scale_aons)), he_avg.lhs))
    filter!(x->SymbolicUtils.istree(x.arguments[1]), reference_terms)

    rw_redundants = _redundant_rewriter(scale_aons,names,reference_terms;kwargs...)
    f = SymbolicUtils.Fixpoint(rw_redundants)

    for i=1:length(he)
        rhs_ = _scale(he_avg.lhs[i], he_avg.rhs[i], scale_aons, N, M, names, reference_terms)
        rhs_ = f(rhs_)
        push!(rhs, _undo_average(rhs_))
    end

    he_scaled = ScaledHeisenbergEquation(he.lhs,rhs,he.hamiltonian,he.jumps,he.rates,
                scale_aons, names, ones(Bool, length(he)), reference_terms
    )

    if simplify
        return qsimplify(he_scaled; kwargs...)
    else
        return he_scaled
    end
end

function _scale(lhs, rhs, scale_aons, N, M, names, reference_terms)
    aon_lhs = acts_on(lhs)
    scale_aon_ = intersect(aon_lhs, scale_aons)
    M_ = length(scale_aon_)
    if length(aon_lhs)==1 && aon_lhs[1]==scale_aons[1]
        N_ = 1
    elseif iszero(M_)
        N_ = N / M
    elseif M_ == M # all clusters on lhs
        N_ = 1
    else
        N_ = (N - (M-M_)) / (M - M_)
    end
    rhs_ = _scaling_rewrite(rhs, N_, scale_aons, aon_lhs, names)
    return rhs_
end

function scale(he::ScaledHeisenbergEquation, N; kwargs...)
    # TODO remove this function?
    M = length(scale_aons)
    rhs = []
    reference_terms = copy(he.reference_terms)
    he_avg = average(he)

    rw_redundants = _redundant_rewriter(he.scale_aons,names,reference_terms;kwargs...)
    f = SymbolicUtils.Fixpoint(rw_redundants)

    for i=1:length(he)
        if he.was_scaled[i]
            push!(rhs, he.rhs[i])
        else
            rhs_ = _scale(he_avg.lhs[i], he_avg.rhs[i], scale_aons, N, M, he.names; kwargs...)
            rhs_ = f(rhs_)
            push!(rhs, _undo_average(rhs_))
        end
    end
    return ScaledHeisenbergEquation(he.lhs, rhs, he.hamiltonian, he.jumps, he.rates,
            he.scale_aons, he.names, ones(Bool, length(he)), reference_terms
    )
end

_undo_average(x::Number) = x
function _undo_average(x::SymbolicUtils.Symbolic)
    if SymbolicUtils.istree(x)
        f = SymbolicUtils.operation(x)
        if f === average
            return SymbolicUtils.arguments(x)[1]
        else
            args = map(_undo_average, SymbolicUtils.arguments(x))
            T = SymbolicUtils._promote_symtype(f, args)
            if T <: QNumber
                return QTerm(f, args)
            else
                return SymbolicUtils.similarterm(x, f, args)
            end
        end
    else
        return x
    end
end

function get_names(he)
    H = he.hamiltonian
    J = he.jumps
    ops = get_operators(H)
    for j ∈ J
        append!(ops, get_operators(j))
    end
    for l ∈ he.lhs
        append!(ops, get_operators(l))
    end

    # The following can fail if an operator for one specific acts_on is missing
    names = Symbol[]
    aon = 1
    idx = findfirst(x->isequal(acts_on(x),aon), ops)
    while !isnothing(idx)
        push!(names, ops[idx].name)
        aon += 1
        idx = findfirst(x->isequal(acts_on(x),aon), ops)
    end
    return names
end

function _scaling_rewrite(rhs, N, scale_aons, aon_lhs, names)
    rw = let was_scaled=UInt[]

        SCALE_RULES = [
            SymbolicUtils.@rule(~x::should_scale(scale_aons,aon_lhs,was_scaled) => N*~x)
        ]

        rule_tree = [
            SymbolicUtils.If(
                SymbolicUtils.is_operation(average), SymbolicUtils.Chain(SCALE_RULES)
                ),
        ] |> SymbolicUtils.Chain
        SymbolicUtils.Postwalk(rule_tree)
    end

    f = SymbolicUtils.Fixpoint(rw)
    return f(rhs)
end

function substitute_redundants(x,scale_aons,names,reference_terms;kwargs...)
    # TODO: do we need this function?
    rw = _redundant_rewriter(scale_aons,names,reference_terms;kwargs...)
    return qsimplify(x;rewriter=rw,kwargs...)
end

function _redundant_rewriter(scale_aons,names,reference_terms;kwargs...)
    rw = let adj=false,idx=nothing

        SUBS_RULES = [
            SymbolicUtils.@rule(~x::is_redundant(scale_aons,reference_terms;idx=idx,adj=adj) => _replace_redundant(~x,scale_aons,names,reference_terms,idx,adj))
        ]

        rule_tree = [
            SymbolicUtils.If(
                SymbolicUtils.is_operation(average), SymbolicUtils.Chain(SUBS_RULES)
            ),
        ] |> SymbolicUtils.Chain

        SymbolicUtils.Postwalk(rule_tree)
    end
    return rw
end

should_scale(scale_aons,aon_lhs,was_scaled) = (x->should_scale(x, scale_aons, aon_lhs, was_scaled))
function should_scale(x, scale_aons, aon_lhs, was_scaled)
    h = hash(x)
    h ∈ was_scaled && return false # x was already scaled
    aon = acts_on(x)
    all(a ∈ aon_lhs for a ∈ aon) && return false # x acts only on things contained in lhs
    !any(a ∈ scale_aons for a ∈ aon) && return false # no scaling should occur if x is not part of a cluster
    should_scale_ = any(a ∈ scale_aons for a ∈ aon)
    should_scale_ && push!(was_scaled, h)
    return should_scale_
end

is_redundant(scale_aons,reference_terms;kwargs...) = x->is_redundant(x,scale_aons,reference_terms;kwargs...)
function is_redundant(x,scale_aons,reference_terms;idx=nothing,adj=false)
    # Judge whether a term is redundant from its acts_on
    aon = acts_on(x)
    idx_aon = findall(in(scale_aons), aon)
    isempty(idx_aon) && return false

    if length(idx_aon)>1
        for i=1:length(idx_aon)-1
            (idx_aon[i]+1 == idx_aon[i+1]) || return true
        end
    end

    aon[idx_aon[1]] == scale_aons[1] || return true

    # For more involved expressions, we need to check against other terms explicitly
    length(idx_aon)==1 && return false
    _in(x, reference_terms) && return false
    _in(_adjoint(x), reference_terms) && return false

    # TODO add adjoint(x) to reference_terms so we can lose the adj Bool?
    adj, idx = _findfirst_redundant(x,scale_aons,reference_terms)
    if isnothing(idx)
        push!(reference_terms, x)
        return false
    else
        return true
    end
end

function _findfirst_redundant(x,scale_aons,reference_terms)
    idx = nothing
    args = SymbolicUtils.arguments(SymbolicUtils.arguments(x)[1])
    args_std = filter(x->acts_on(x)∉scale_aons, args)
    args_cluster = filter(x->acts_on(x)∈scale_aons, args)
    adj = false
    for i=1:length(reference_terms)
        adj, _is_equal = _cluster_compare(args_std, args_cluster, reference_terms[i], scale_aons)
        if _is_equal
            idx = i
            break
        end
    end
    return adj, idx
end

function _cluster_compare(args_std, args_cluster, ref_term, scale_aons)
    # Split arguments in non-cluster (std) and clusters
    ref_args = SymbolicUtils.arguments(SymbolicUtils.arguments(ref_term)[1])
    ref_args_std = filter(x->acts_on(x)∉scale_aons, ref_args)
    length(args_std)==length(ref_args_std) || return (false,false)
    ref_args_cluster = filter(x->acts_on(x)∈scale_aons, ref_args)
    length(ref_args_cluster) == length(args_cluster) || return (false,false)

    # Compare non-cluster arguments first
    check_std = true
    for i=1:length(args_std)
        check_std = isequal(args_std[i], ref_args_std[i])
        check_std || break
    end

    # If it fails, see if the adjoints match
    if check_std
        check_std_adj = false
    else
        check_std_adj = true
        for i=1:length(args_std)
            check_std_adj = isequal(_adjoint(args_std[i]), ref_args_std[i])
            check_std_adj || break
        end
    end
    (check_std || check_std_adj) || return (false,false)

    # Only if the non-cluster arguments match, check the cluster arguments
    if check_std
        return false, _check_redundant_cluster_args(args_cluster, ref_args_cluster)
    else
        args_cluster_adj = map(_adjoint, args_cluster)
        return true, _check_redundant_cluster_args(args_cluster_adj, ref_args_cluster)
    end
end

function _check_redundant_cluster_args(cluster::Vector{<:Transition}, refs)
    check_cluster = true
    p = permutations(cluster)
    for arg ∈ p
        for i=1:length(cluster)
            check_cluster = _compare_without_name_or_aon(arg[i], refs[i])
            check_cluster || break
        end
        check_cluster && break
    end
    return check_cluster
end
_compare_without_name_or_aon(a::Transition,b::Transition) = (a.i==b.i && a.j==b.j)

function _check_redundant_cluster_args(cluster, refs)
    @assert cluster[1] isa Destroy || cluster[1] isa Create # TODO handle by dispatch
    # To avoid false positives of e.g. a_1'*a_1*a_2 compared to a_1*a_1*a_2'
    # we need to count the number of Create/Destroy per acts_on

    # TODO: clean up this mess here

    destroys_c, creates_c = _count_destroy_create_per_aon(cluster)
    destroys_r, creates_r = _count_destroy_create_per_aon(refs)

    length(destroys_c) == length(destroys_r) || return false
    length(creates_c) == length(creates_r) || return false

    # Check destroys
    check = true
    p = permutations(destroys_c)
    for arg ∈ p
        for i=1:length(destroys_c)
            check = isequal(arg[i], destroys_r[i])
            check || break
        end
        check && break
    end
    check || return false

    p = permutations(creates_c)
    for arg ∈ p
        for i=1:length(creates_c)
            check = isequal(arg[i], creates_r[i])
            check || break
        end
        check && break
    end

    return check
end

function _count_destroy_create_per_aon(args)
    num_destroy = Int[]
    num_create = Int[]

    aon_i = 0
    idx_count = 0
    idx_args = 1
    while idx_args <= length(args)
        if aon_i == acts_on(args[idx_args])
            num_destroy[idx_count] += args[idx_args] isa Destroy
            num_create[idx_count] += args[idx_args] isa Create
            idx_args += 1
        else
            push!(num_destroy, 0)
            push!(num_create, 0)
            aon_i = copy(acts_on(args[idx_args]))
            idx_count += 1
        end
    end

    return num_destroy, num_create
end

function _replace_redundant(x,scale_aons,names,reference_terms,idx,adj)
    if isnothing(idx) # Replace by changing acts_on
        aon = acts_on(x)
        idx_aon = findall(in(scale_aons), aon)
        if length(aon)==1
            aon_subs = scale_aons[idx_aon[1]]
        else
            aon_subs = copy(aon)
            aon_subs[idx_aon] .= scale_aons[1:length(idx_aon)]
        end
        return _swap_aon_and_name(x, aon, aon_subs, names)
    else
        if adj
            return _adjoint(reference_terms[idx])
        else
            return reference_terms[idx]
        end
    end
end

_swap_aon_and_name(x::Average, aon1, aon2, names) = _average(_swap_aon_and_name(x.arguments[1], aon1, aon2, names))
function _swap_aon_and_name(t::QTerm{<:typeof(*)}, aon1, aon2, names)
    args = []
    for arg in SymbolicUtils.arguments(t)
        idx = findfirst(isequal(acts_on(arg)), aon1)
        push!(args, _swap_aon_and_name(arg, aon1[idx], aon2[idx], names[idx]))
    end
    sort!(args, by=acts_on)
    return QTerm(*, args)
end
_swap_aon_and_name(op::T, aon1, aon2, name::Symbol) where T<:QSym = T(op.hilbert, name, aon2)
function _swap_aon_and_name(op::Transition, aon1, aon2, name::Symbol)
    Transition(op.hilbert, name, op.i, op.j, aon2)
end
_swap_aon_and_name(op::QSym, aon1, aon2, names) = _swap_aon_and_name(op, aon1, aon2, names[aon2])
