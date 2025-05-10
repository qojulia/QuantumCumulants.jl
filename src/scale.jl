function get_cluster_stuff(h::ClusterSpace, aon=1)
    M = h.order
    aons = [ClusterAon(aon, i) for i in 1:M]
    return [aons], [h.N]
end
function get_cluster_stuff(h::ProductSpace)
    idx = findall(x->isa(x, ClusterSpace), h.spaces)
    aons = Vector{<:ClusterAon}[]
    N = []
    for i in idx
        aon_, N_ = get_cluster_stuff(h.spaces[i], i)
        append!(aons, aon_)
        append!(N, N_)
    end
    return aons, N
end

function _scale(lhs, rhs, scale_aons, N, M, names)
    aon_lhs = acts_on(lhs)
    scale_aon_ = intersect(aon_lhs, scale_aons)
    M_ = length(scale_aon_)
    if iszero(M_)
        N_ = N / M
    elseif M_ == M # all clusters on lhs
        N_ = 1
    else
        N_ = (N - (M-M_)) / (M - M_)
    end
    rhs = _scaling_rewrite(rhs, N_, scale_aons, aon_lhs, names)
    return rhs
end

function _get_names_ops(ops)
    hs_ = hilbert(ops[1])
    if isa(hs_, ProductSpace)
        hs = hilbert(ops[1]).spaces
    else
        hs = [hs_]
    end
    names = Vector{Any}(undef, length(hs))
    for op in ops
        aon = acts_on(op)
        if isa(aon, Int)
            names[aon] = op.name
        else #ClusterAon
            h = hs[aon.i]
            order = h.order
            op_name = h.op_name[]
            cluster_names = [Symbol(op_name, :_, i) for i in 1:order]
            names[aon.i] = cluster_names
        end
    end
    return names
end

function get_names(q::Union{QSym,QMul})
    ops = get_operators(q)
    unique_ops!(ops)
    _get_names_ops(ops)
end
function get_names(he)
    H = he.hamiltonian
    J = he.jumps
    ops = get_operators(H)
    for j in J
        append!(ops, get_operators(j))
    end
    for l in he.operators
        append!(ops, get_operators(l))
    end
    unique_ops!(ops)
    _get_names_ops(ops)
end

function unique_i_aons(aon)
    seen = eltype(aon)[]
    for a in aon
        if !(a∈seen)
            if a isa Integer
                push!(seen, a)
            else # ClusterAon
                isone(a.j) && push!(seen, a)
            end
        end
    end
    return seen
end

get_operators(q::QSym) = [q]
function get_operators(q::QMul)
    ops = QSym[]
    seen_hashes = UInt[]
    for arg in q.args_nc
        h = hash(arg)
        if !(h ∈ seen_hashes)
            push!(seen_hashes, h)
            push!(ops, arg)
        end
    end
    return ops
end
function get_operators(q::QAdd)
    ops = QSym[]
    for arg in q.arguments
        append!(ops, get_operators(arg))
    end
    unique_ops!(ops)
    return ops
end

## Scaling terms by the correct factor

function _scaling_rewrite(rhs, N, scale_aons, aon_lhs, names)
    rw = let was_scaled=UInt[]
        SCALE_RULES = [
            SymbolicUtils.@rule(~x::should_scale(scale_aons, aon_lhs, was_scaled) => N*~x)
        ]

        rule_tree = SymbolicUtils.Chain([
            SymbolicUtils.If(
                SymbolicUtils.is_operation(sym_average),
                SymbolicUtils.Chain(SCALE_RULES),
            ),
        ])
        SymbolicUtils.Postwalk(rule_tree)
    end

    f = SymbolicUtils.Fixpoint(rw)
    return f(rhs)
end

function should_scale(scale_aons, aon_lhs, was_scaled)
    (x->should_scale(x, scale_aons, aon_lhs, was_scaled))
end
function should_scale(x, scale_aons, aon_lhs, was_scaled)
    h = hash(x)
    h ∈ was_scaled && return false # x was already scaled
    aon = [acts_on(x)...]
    intersect!(aon, scale_aons)
    all(a in aon_lhs for a in aon) && return false # x acts only on things contained in lhs
    should_scale_ = any(a in scale_aons for a in aon) # no scaling should occur if x is not part of a cluster
    should_scale_ && push!(was_scaled, h)
    return should_scale_
end

## Dealing with redundant averages

function substitute_redundants(
    t::T, scale_aons::Vector{<:Vector}, names
) where {T<:SymbolicUtils.Symbolic}
    if SymbolicUtils.iscall(t)
        f = SymbolicUtils.operation(t)
        if f === sym_average
            op = deepcopy(SymbolicUtils.arguments(t)[1])
            for j in 1:length(scale_aons)
                op = substitute_redundants(op, scale_aons[j], names)
            end
            return average(op)
        else
            args = []
            for arg in SymbolicUtils.arguments(t)
                push!(args, substitute_redundants(arg, scale_aons, names))
            end
            return SymbolicUtils.maketerm(T, f, args, TermInterface.metadata(t))
        end
    else
        return t
    end
end
substitute_redundants(x::Number, args...) = x
substitute_redundants(x, scale_aons, names) = substitute_redundants(x, [scale_aons], names)

function substitute_redundants(t::QMul, scale_aons, names)
    aon = acts_on(t)
    idx_aon = findall(in(scale_aons), aon)
    isempty(idx_aon) && return t
    if is_redundant_aon(t, scale_aons)
        if length(aon)==1
            aon_subs = scale_aons[idx_aon[1]]
        else
            aon_subs = copy(aon)
            aon_subs[idx_aon] .= scale_aons[1:length(idx_aon)]
        end
        name_idx = map(get_i, aon_subs)
        return substitute_redundants(
            _swap_aon_and_name(t, aon, aon_subs, names[name_idx]), scale_aons, names
        )
    else
        args = t.args_nc
        args_cluster = filter(x->acts_on(x)∈scale_aons, args)

        # Get the proper ordering
        p = sortperm_ref_order(args_cluster)

        if issorted(p) # arguments are already in correct order
            return t
        else
            # acts_on and relevant names
            aon_subs = copy(aon)
            names_ = names[map(get_i, aon)]

            # Permute cluster part to proper reference order
            aon_subs[idx_aon[p]] = scale_aons[1:length(idx_aon)]
            names_[idx_aon] = names_[idx_aon][p]

            # Swap and return
            return _swap_aon_and_name(t, aon, aon_subs, names_)
        end
    end
end
function substitute_redundants(x::QSym, scale_aons, names)
    if is_redundant_aon(x, scale_aons)
        aon = acts_on(x)
        if aon isa ClusterAon
            i = aon.i
            idx_aon = findfirst(x->x.i==i, scale_aons)
        else
            idx_aon = findfirst(in(scale_aons), aon)
        end
        aon_sub = scale_aons[idx_aon]
        return _swap_aon_and_name(x, aon, aon_sub, names)
    else
        return x
    end
end

function is_redundant_aon(x, scale_aons)
    # Judge whether a term is redundant from its acts_on
    aon = acts_on(x)
    if aon isa ClusterAon
        aon ∈ scale_aons || return false
        return aon != scale_aons[1]
    else
        idx_aon = findall(in(aon), scale_aons)
        isempty(idx_aon) && return false
        for i in 1:(length(idx_aon) - 1)
            (idx_aon[i]+1 == idx_aon[i + 1]) || return true
        end
        idx = findfirst(in(scale_aons), aon)
        return aon[idx] != scale_aons[1]
    end
end

function sortperm_ref_order(args_cluster)
    if args_cluster[1] isa Transition
        return sortperm(args_cluster; lt=lt_reference_order)
    else
        # non-unique acts_on
        arg_aon = map(acts_on, args_cluster)

        # Separate arguments into blocks with equal acts_on
        i = 1
        blocks = []
        while i <= length(args_cluster)
            idx = findall(isequal(arg_aon[i]), arg_aon[i:end])
            n = length(idx)
            block = args_cluster[i:(i + n - 1)]
            push!(blocks, block)
            i += n
        end

        # Sort blocks by length and count of Destroy/Create
        return sortperm(blocks; lt=_lt_num_destroy_create)
    end
end

function _lt_num_destroy_create(args1, args2)
    if length(args1) != length(args2) # Sort by length first
        return length(args1) > length(args2)
    else # Equal length is sorted by number of Destroy (lower ones to the right)
        f = x->isa(x, Destroy)
        idx_destroys1 = findall(f, args1)
        ndest1 = length(idx_destroys1)
        idx_destroys2 = findall(f, args2)
        ndest2 = length(idx_destroys2)
        return ndest1 <= ndest2
    end
end

function lt_reference_order(t1::Transition, t2::Transition)
    aon = acts_on(t1)
    isequal(aon, acts_on(t2)) && return false
    lvls = levels(t1.hilbert, aon)
    f = x->findfirst(isequal(x), lvls)
    i1, j1 = f(t1.i), f(t1.j)
    i2, j2 = f(t2.i), f(t2.j)
    d1 = abs(i1 - j1)
    d2 = abs(i2 - j2)
    if d1==d2
        if i1==i2
            return j1 <= j2
        else
            m1 = min(i1, j1)
            m2 = min(i2, j2)
            if m1==m2
                return i1 > i2
            else
                return m1 < m2
            end
        end
    else
        return d1 < d2
    end
end

function _swap_aon_and_name(x::Average, aon1, aon2, names)
    _average(_swap_aon_and_name(x.arguments[1], aon1, aon2, names))
end
function _swap_aon_and_name(t::QMul, aon1, aon2, names)
    args = QSym[]
    for arg in t.args_nc
        idx = findfirst(isequal(acts_on(arg)), aon1)
        push!(args, _swap_aon_and_name(arg, aon1[idx], aon2[idx], names[idx]))
    end
    sort!(args; by=acts_on)
    return QMul(t.arg_c, args)
end
for T in (:Create, :Destroy)
    @eval _swap_aon_and_name(op::$(T), aon1, aon2, name::Symbol) = $(T)(
        op.hilbert, name, aon2; op.metadata
    )
end
function _swap_aon_and_name(op::Transition, aon1, aon2, name::Symbol)
    Transition(op.hilbert, name, op.i, op.j, aon2; op.metadata)
end
function _swap_aon_and_name(op::QSym, aon1, aon2, names)
    _swap_aon_and_name(op, aon1, aon2, _get_names(names, aon2))
end

_get_names(names, aon::Integer) = names[aon]
_get_names(names::Vector{<:Symbol}, aon::ClusterAon) = names[aon.j]
_get_names(names, aon::ClusterAon) = names[aon.i][aon.j]
function _get_names(names, aon)
    names_ = Symbol[]
    for a in aon
        n = [_get_names(names, a1) for a1 in a]
        append!(names_, n)
    end
end

## Complete
function filter_redundants!(
    missed, scale_aons_::Vector{<:Vector}, names, vhash=UInt[], vs′hash=UInt[]
)
    j = 1
    while j <= length(missed)
        for scale_aons in scale_aons_
            missed[j] = substitute_redundants(missed[j], scale_aons, names)
        end
        h = hash(missed[j])
        h′ = hash(_adjoint(missed[j]))
        if h∈vhash || h′∈vhash || h∈vs′hash || h′∈vs′hash
            deleteat!(missed, j)
        else
            j += 1
            push!(vhash, h)
        end
    end
    return missed
end
