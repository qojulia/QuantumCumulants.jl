function scale(he::HeisenbergEquation, scale_aons, N; kwargs...)
    M = length(scale_aons)
    rhs = []
    names = get_names(he)
    for i=1:length(he)
        rhs_ = _scale(he.lhs[i], he.rhs[i], scale_aons, N, M, names; kwargs...)
        push!(rhs, rhs_)
    end
    return ScaledHeisenbergEquation(he.lhs, rhs, he.hamiltonian, he.jumps, he.rates, ones(Bool, length(he)))
end

function _scale(lhs, rhs, scale_aons, N, M, names; kwargs...)
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
    rw = scaling_rewriter(N_, scale_aons, aon_lhs, names; kwargs...)
    rhs_ = qsimplify(rhs; rewriter=rw)
    rhs_ = substitute_redundants(rhs_, scale_aons, names)
    return rhs_
end

function scale(he::ScaledHeisenbergEquation, scale_aons, N; kwargs...)
    M = length(scale_aons)
    rhs = []
    names = get_names(he)
    for i=1:length(he)
        if he.was_scaled[i]
            push!(rhs, he.rhs[i])
        else
            rhs_ = _scale(he.lhs[i], he.rhs[i], scale_aons, N, M, names; kwargs...)
            push!(rhs, rhs_)
        end
    end
    return ScaledHeisenbergEquation(he.lhs, rhs, he.hamiltonian, he.jumps, he.rates, ones(Bool, length(he)))
end

# _undo_average(x::Number) = x
# function _undo_average(x::SymbolicUtils.Symbolic)
#     if SymbolicUtils.istree(x)
#         f = SymbolicUtils.operation(x)
#         if f === average
#             return SymbolicUtils.arguments(x)[1]
#         else
#             args = map(_undo_average, SymbolicUtils.arguments(x))
#             return SymbolicUtils.similarterm(x, f, args)
#         end
#     else
#         return x
#     end
# end

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

function scaling_rewriter(N, scale_aons, aon_lhs, names; kwargs...)

    was_scaled = []

    SCALE_RULES = [
        SymbolicUtils.@rule(~x::should_scale(scale_aons,aon_lhs,was_scaled) => N*~x)
        SymbolicUtils.@rule(~x::SymbolicUtils.isnotflat(*) => SymbolicUtils.flatten_term(*,~x))
    ]

    rule_tree = [
        SymbolicUtils.If(
            SymbolicUtils.is_operation(*), SymbolicUtils.Chain(SCALE_RULES)
            ),
    ] |> SymbolicUtils.Chain

    return SymbolicUtils.Postwalk(rule_tree)
end

function substitute_redundants(x,scale_aons,names;kwargs...)

    SUBS_RULES = [
        SymbolicUtils.@rule(~x::is_redundant(scale_aons) => _replace_redundant(~x,scale_aons,names))
    ]

    rule_tree = [
        SymbolicUtils.If(
            SymbolicUtils.is_operation(*), SymbolicUtils.Chain(SUBS_RULES)
        ),
        default_operator_simplifier(;kwargs...)
    ] |> SymbolicUtils.Chain

    rw = SymbolicUtils.Postwalk(rule_tree)

    return qsimplify(x;rewriter=rw)
end

should_scale(scale_aons,aon_lhs,was_scaled) = (x->should_scale(x, scale_aons, aon_lhs, was_scaled))
function should_scale(x, scale_aons, aon_lhs, was_scaled)
    h = hash(filter(!iscommutative, SymbolicUtils.arguments(x)))
    h ∈ was_scaled && return false # x was already scaled
    aon = acts_on(x)
    all(a ∈ aon_lhs for a ∈ aon) && return false # x acts only on things contained in lhs
    !any(a ∈ scale_aons for a ∈ aon) && return false # no scaling should occur if x is not part of a cluster
    should_scale_ = any(a ∈ scale_aons for a ∈ aon)
    should_scale_ && push!(was_scaled, h)
    return should_scale_
end

is_redundant(scale_aons) = x->is_redundant(x,scale_aons)
function is_redundant(x,scale_aons)
    aon = acts_on(x)
    idx = findall(in(scale_aons), aon)
    isempty(idx) && return false

    if length(idx)>1
        for i=1:length(idx)-1
            (idx[i]+1 == idx[i+1]) || return true
        end
    end

    return aon[idx[1]] != scale_aons[1]
end

function _replace_redundant(x,scale_aons,names)
    aon = acts_on(x)
    idx = findall(in(scale_aons), aon)
    if length(aon)==1
        aon_subs = scale_aons[idx[1]]
    else
        aon_subs = copy(aon)
        aon_subs[idx] .= scale_aons[1:length(idx)]
    end
    return _swap_aon_and_name(x, aon, aon_subs, names)
end

function _swap_aon_and_name(t::QTerm{<:typeof(*)}, aon1, aon2, names)
    args = []
    args_c = filter(iscommutative, SymbolicUtils.arguments(t))
    args_nc = filter(!iscommutative, SymbolicUtils.arguments(t))
    for arg in args_nc
        idx = findfirst(isequal(acts_on(arg)), aon1)
        push!(args, _swap_aon_and_name(arg, aon1[idx], aon2[idx], names[idx]))
    end
    sort!(args, by=acts_on)
    return QTerm(*, vcat(args_c, args))
end
_swap_aon_and_name(op::T, aon1, aon2, name::Symbol) where T<:QSym = T(op.hilbert, name, aon2)
function _swap_aon_and_name(op::Transition, aon1, aon2, name::Symbol)
    Transition(op.hilbert, name, op.i, op.j, aon2)
end
_swap_aon_and_name(op::QSym, aon1, aon2, names) = _swap_aon_and_name(op, aon1, aon2, names[aon2])
