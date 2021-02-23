function scale(he::HeisenbergEquation, scale_aons, N; simplify=true, kwargs...)
    M = length(scale_aons)
    rhs = QNumber[]
    names = get_names(he)
    he_avg_ = average(he) # Ensure correct form of expressions required for scaling
    he_avg = substitute_redundants(he_avg_, scale_aons, names)

    for i=1:length(he)
        rhs_ = _scale(he_avg.lhs[i], he_avg.rhs[i], scale_aons, N, M, names)
        rhs_ = substitute_redundants(rhs_, scale_aons, names)
        push!(rhs, _undo_average(rhs_))
    end

    he_scaled = ScaledHeisenbergEquation(he.lhs,rhs,he.hamiltonian,he.jumps,he.rates,
                scale_aons, names, ones(Bool, length(he))
    )

    if simplify
        return qsimplify(he_scaled; kwargs...)
    else
        return he_scaled
    end
end

function _scale(lhs, rhs, scale_aons, N, M, names)
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


## Scaling terms by the correct factor

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

should_scale(scale_aons,aon_lhs,was_scaled) = (x->should_scale(x, scale_aons, aon_lhs, was_scaled))
function should_scale(x, scale_aons, aon_lhs, was_scaled)
    h = hash(x)
    h ∈ was_scaled && return false # x was already scaled
    aon = acts_on(x)
    all(a ∈ aon_lhs for a ∈ aon) && return false # x acts only on things contained in lhs
    should_scale_ = any(a ∈ scale_aons for a ∈ aon) # no scaling should occur if x is not part of a cluster
    should_scale_ && push!(was_scaled, h)
    return should_scale_
end


## Dealing with redundant averages

function substitute_redundants(he::HeisenbergEquation, args...)
    lhs = []
    rhs = []
    for i=1:length(he)
        lhs_ = substitute_redundants(he.lhs[i], args...)
        rhs_ = substitute_redundants(he.rhs[i], args...)
        push!(lhs, lhs_)
        push!(rhs, rhs_)
    end
    # TODO substitute jumps and hamiltonian?
    return HeisenbergEquation(lhs, rhs, he.hamiltonian, he.jumps, he.rates)
end

function substitute_redundants(t::SymbolicUtils.Symbolic,scale_aons,names)
    if SymbolicUtils.istree(t)
        f = SymbolicUtils.operation(t)
        if f === average
            op = SymbolicUtils.arguments(t)[1]
            op_ = substitute_redundants(op,scale_aons,names)
            return average(op_)
        else
            args = []
            for arg in SymbolicUtils.arguments(t)
                push!(args, substitute_redundants(arg,scale_aons,names))
            end
            return SymbolicUtils.similarterm(t, f, args)
        end
    else
        return t
    end
end
substitute_redundants(x::Number,args...) = x

function substitute_redundants(t::QTerm{<:typeof(*)},scale_aons,names)
    aon = acts_on(t)
    idx_aon = findall(in(scale_aons), aon)
    isempty(idx_aon) && return t
    if is_redundant_aon(t,scale_aons)
        if length(aon)==1
            aon_subs = scale_aons[idx_aon[1]]
        else
            aon_subs = copy(aon)
            aon_subs[idx_aon] .= scale_aons[1:length(idx_aon)]
        end
        return substitute_redundants(
                            _swap_aon_and_name(t, aon, aon_subs, names),
                            scale_aons, names
                            )
    else
        args = SymbolicUtils.arguments(t)
        args_cluster = filter(x->acts_on(x)∈scale_aons, args)

        # Get the proper ordering
        p = sortperm_ref_order(args_cluster)

        if issorted(p) # arguments are already in correct order
            return t
        else
            # acts_on and relevant names
            aon_subs = copy(aon)
            names_ = names[aon]

            # Permute cluster part to proper reference order
            aon_subs[idx_aon] = scale_aons[p]
            names_[idx_aon] = names_[idx_aon][p]

            # Swap and return
            return _swap_aon_and_name(t, aon, aon_subs, names_)
        end
    end
end
function substitute_redundants(x::QSym,scale_aons,names)
    if is_redundant_aon(x,scale_aons)
        aon = acts_on(x)
        idx_aon = findfirst(in(scale_aons), aon)
        aon_sub = scale_aons[idx_aon]
        return _swap_aon_and_name(x,aon,aon_sub,names)
    else
        return x
    end
end

function is_redundant_aon(x,scale_aons)
    # Judge whether a term is redundant from its acts_on
    aon = acts_on(x)
    idx_aon = findall(in(scale_aons), aon)
    isempty(idx_aon) && return false

    for i=1:length(idx_aon)-1
        (idx_aon[i]+1 == idx_aon[i+1]) || return true
    end

    return aon[idx_aon[1]] != scale_aons[1]
end

function sortperm_ref_order(args_cluster)
    if args_cluster[1] isa Transition
        return sortperm(args_cluster, lt=lt_reference_order)
    else
        # non-unique acts_on
        arg_aon = map(acts_on, args_cluster)

        # Separate arguments into blocks with equal acts_on
        i = 1
        blocks = []
        while i <= length(args_cluster)
            idx = findall(isequal(arg_aon[i]),arg_aon[i:end])
            n = length(idx)
            block = args_cluster[i:i+n-1]
            push!(blocks, block)
            i += n
        end

        # Sort blocks by length and count of Destroy/Create
        return sortperm(blocks, lt=_lt_num_destroy_create)
    end
end

function _lt_num_destroy_create(args1, args2)
    if length(args1) != length(args2) # Sort by length first
        return length(args1) > length(args2)
    else # Equal length is sorted by number of Destroy (lower ones to the right)
        f = x->isa(x,Destroy)
        idx_destroys1 = findall(f,args1)
        ndest1 = length(idx_destroys1)
        idx_destroys2 = findall(f,args2)
        ndest2 = length(idx_destroys2)
        return ndest1 <= ndest2
    end
end

function lt_reference_order(t1::Transition, t2::Transition)
    isequal(acts_on(t1),acts_on(t2)) && return false
    if t1.i < t2.i
        return true
    elseif t1.i == t2.i
        return t1.j < t2.j
    else
        return t1.j > t2.j
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
