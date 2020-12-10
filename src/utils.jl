"""
    find_missing(rhs::Vector, vs::Vector, vs_adj=adjoint.(vs), ps=[])

For a list of expressions contained in `rhs`, check whether all occurring symbols
are contained either in the variables given in `vs`. If a list of parameters `ps`
is provided, parameters that do not occur in the list `ps` are also added to the list.
Returns a list of missing symbols.
"""
function find_missing(rhs::Vector{<:Number}, vs::Vector{<:Number}; vs_adj::Vector=adjoint.(vs), ps=[])
    missed = Number[]
    for e=rhs
        append!(missed,get_symbolics(e))
    end
    unique!(missed)
    if isempty(ps)
        filter!(x->!isa(x,Parameter), missed)
    end
    filter!(x->!(x∈vs || x∈ps || x∈vs_adj),missed)
    isempty(ps) || (ps_adj = adjoint.(ps); filter!(x -> !(x∈ps_adj), missed))
    return missed
end
function find_missing(de::DifferentialEquation{<:Number,<:Number}; kwargs...)
    find_missing(de.rhs, de.lhs; kwargs...)
end

"""
    get_symbolics(ex)

Find all symbolic numbers occuring in `ex`.
"""
get_symbolics(x::Number) = SymbolicNumber[]
get_symbolics(x::SymbolicNumber) = [x]
function get_symbolics(t::NumberTerm)
    syms = SymbolicNumber[]
    for arg in t.arguments
        append!(syms, get_symbolics(arg))
    end
    return unique(syms)
end

"""
    complete(de::DifferentialEquation)

From a set of differential equation of averages, find all averages that are missing
and derive the corresponding equations of motion.
"""
function complete(de::DifferentialEquation{<:Number,<:Number};kwargs...)
    rhs_, lhs_ = complete(de.rhs,de.lhs,de.hamiltonian,de.jumps,de.rates;kwargs...)
    return DifferentialEquation(lhs_,rhs_,de.hamiltonian,de.jumps,de.rates)
end
function complete(rhs::Vector{<:Number}, vs::Vector{<:Number}, H, J, rates; order=nothing, filter_func=nothing, mix_choice=maximum, kwargs...)
    order_lhs = maximum(get_order.(vs))
    order_rhs = maximum(get_order.(rhs))
    if order isa Nothing
        order_ = max(order_lhs, order_rhs)
    else
        order_ = order
    end
    maximum(order_) >= order_lhs || error("Cannot form cumulant expansion of derivative; you may want to use a higher order!")

    vs_ = copy(vs)
    rhs_ = [cumulant_expansion(r, order_) for r in rhs]
    missed = unique_ops(find_missing(rhs_, vs_))
    filter!(x->isa(x,Average),missed)
    isnothing(filter_func) || filter!(filter_func, missed) # User-defined filter
    while !isempty(missed)
        ops = getfield.(missed, :operator)
        he = isempty(J) ? heisenberg(ops,H) : heisenberg(ops,H,J;rates=rates)
        he_avg = average(he,order_;mix_choice=mix_choice)
        rhs_ = [rhs_;he_avg.rhs]
        vs_ = [vs_;he_avg.lhs]
        missed = unique_ops(find_missing(rhs_,vs_))
        filter!(x->isa(x,Average),missed)
        isnothing(filter_func) || filter!(filter_func, missed) # User-defined filter
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
    return rhs_, vs_
end

"""
    find_operators(::HilbertSpace, order; names=nothing)

Find all operators that fully define a system up to the given `order`.
"""
function find_operators(h::HilbertSpace, order::Int; names=nothing, kwargs...)
    if names isa Nothing && (unique(typeof.(h.spaces))!=typeof.(h.spaces))
        alph = 'a':'z'
        names_ = Symbol.(alph[1:length(h.spaces)])
    else
        names_ = names
    end
    fund_ops = fundamental_operators(h;names=names_, kwargs...)
    fund_ops = unique([fund_ops;adjoint.(fund_ops)])
    ops = copy(fund_ops)
    for i=2:order
        ops = [ops;fund_ops]
    end

    all_ops = AbstractOperator[]
    for i=1:order
        for c in combinations(ops, i)
            push!(all_ops, prod(c))
        end
    end

    # Simplify and remove non-operators iteratively
    ops_1 = simplify_operators.(all_ops)
    ops_2 = all_ops
    while ops_1 != ops_2
        ops_2 = AbstractOperator[]
        for op in ops_1
            append!(ops_2, _get_operators(op))
        end
        ops_1 = simplify_operators.(ops_2)
    end

    return unique_ops(ops_2)
end
find_operators(op::AbstractOperator,args...) = find_operators(hilbert(op),args...)

"""
    hilbert(::AbstractOperator)

Return the Hilbert space of the operator.
"""
hilbert(op::BasicOperator) = op.hilbert
hilbert(t::OperatorTerm) = hilbert(t.arguments[findfirst(x->isa(x,AbstractOperator), t.arguments)])

"""
    fundamental_operators(::HilbertSpace)

Return all fundamental operators for a given Hilbertspace. For example,
a [`FockSpace`](@ref) only has one fundamental operator, `Destroy`.
"""
function fundamental_operators(h::FockSpace,aon::Int=1;names=nothing)
    name = names isa Nothing ? :a : names[aon]
    a = Destroy(h,name)
    return [a]
end
function fundamental_operators(h::NLevelSpace,aon::Int=1;names=nothing)
    sigmas = Transition[]
    lvls = levels(h)
    name = names isa Nothing ? :σ : names[aon]
    for i=1:length(lvls)
        for j=i:length(lvls)
            (i==j) && lvls[i]==ground_state(h) && continue
            s = Transition(h,name,lvls[i],lvls[j])
            push!(sigmas,s)
        end
    end
    return sigmas
end
function fundamental_operators(h::ProductSpace;kwargs...)
    ops = []
    for i=1:length(h.spaces)
        ops_ = fundamental_operators(h.spaces[i],i;kwargs...)
        ops_ = [embed(h,o,i) for o in ops_]
        append!(ops,ops_)
    end
    return ops
end


"""
    get_operators(::AbstractOperator)

Return a list of all [`BasicOperator`](@ref) in an expression.
"""
get_operators(x) = _get_operators(x)
function get_operators(t::OperatorTerm{<:typeof(*)})
    ops = AbstractOperator[]
    for arg in t.arguments
        append!(ops, get_operators(arg))
    end
    return ops
end

_get_operators(::Number) = []
_get_operators(op::BasicOperator) = [op]
_get_operators(op::OperatorTerm{<:typeof(^)}) = [op]
function _get_operators(op::OperatorTerm{<:typeof(*)})
    args = AbstractOperator[]
    for arg in op.arguments
        append!(args, _get_operators(arg))
    end
    isempty(args) && return args
    return [*(args...)]
end
function _get_operators(t::OperatorTerm)
    ops = AbstractOperator[]
    for arg in t.arguments
        append!(ops, _get_operators(arg))
    end
    return ops
end

"""
    unique_ops(ops)

For a given list of operators, return only unique ones taking into account
their adjoints.
"""
function unique_ops(ops)
    seen = eltype(ops)[]
    for op in ops
        if !(op in seen || op' in seen)
            push!(seen, op)
        end
    end
    return seen
end


_to_expression(x::Number) = x
function _to_expression(x::Complex) # For brackets when using latexify
    iszero(x) && return x
    if iszero(real(x))
        return :( $(imag(x))*im )
    elseif iszero(imag(x))
        return real(x)
    else
        return :( $(real(x)) + $(imag(x))*im )
    end
end
_to_expression(op::BasicOperator) = op.name
_to_expression(op::Create) = :(dagger($(op.name)))
_to_expression(op::Transition) = :(Transition($(op.name),$(op.i),$(op.j)) )
_to_expression(t::Union{OperatorTerm,NumberTerm}) = :( $(Symbol(t.f))($(_to_expression.(t.arguments)...)) )
_to_expression(p::Parameter) = p.name
function _to_expression(avg::Average)
    ex = _to_expression(avg.operator)
    return :(AVERAGE($ex))
end


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
function scale(de::DifferentialEquation, identical_ops::Vector, interaction_op::AbstractOperator, N::Number)
    identical_aons = unique(acts_on.(identical_ops))
    interaction_aon = acts_on(interaction_op)
    scale(de::DifferentialEquation, identical_aons, interaction_aon, N)
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
            N_sc = (N-length(avg_aon_lhs_all_id))/(length(identical_aons) - length(avg_aon_lhs_all_id)) #double sum extension to RHS??
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
    he_scale = DifferentialEquation(de_lhs, de_rhs, de.hamiltonian, de.jumps, de.rates)
    return filter_redundant_all_id_aon(he_scale, identical_aons, names)
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
function get_redundant_all_id_aon(lhs, all_id_aon, names)
    lhsf = filter(x->length(aon_in_all_id_aon(x, all_id_aon))>1, lhs)
    lhs_redundants = []
    dict_lhs_redunant = Dict()
    while length(lhsf) > 1
        len_lhsf = length(lhsf)
        lhsf1 = lhsf[1]
        id_avgs = get_all_id_permuted_avgs(lhsf1, all_id_aon, names)[2:end]
        id_avgs_adj = adjoint.(id_avgs)
        for it=2:len_lhsf
            if lhsf[it] ∈ [id_avgs; id_avgs_adj]
                push!(lhs_redundants, lhsf[it])
                [dict_lhs_redunant[id_avg] = lhsf1 for id_avg in id_avgs]
            end
        end
        filter!(x->x∉[lhsf1;lhs_redundants], lhsf)
    end
    return lhs_redundants, dict_lhs_redunant
end
function filter_redundant_all_id_aon(de::DifferentialEquation, all_id_aon, names)
    lhs_ = de.lhs; rhs_ = de.rhs
    lhsf = eltype(lhs_)[]; rhsf = eltype(rhs_)[]
    reds, dict_reds = get_redundant_all_id_aon(lhs_, all_id_aon, names)
    for it=1:length(lhs_)
        lhs_[it] ∉ reds && (push!(lhsf, lhs_[it]); push!(rhsf, substitute(rhs_[it], dict_reds)))
    end
    DifferentialEquation(lhsf, rhsf, de.hamiltonian, de.jumps, de.rates)
end

function get_names(de::DifferentialEquation)
    lhs_ops = [lhs_.operator  for lhs_ in de.lhs]
    ops = unique(collect(Iterators.flatten(get_operators.(lhs_ops))))
    sort!(ops; by=acts_on)
    indices = [findfirst(x->acts_on(x)==aon_, ops) for aon_=acts_on(ops[1]):acts_on(ops[end])]
    ops = ops[indices]
    names = [op.name for op in ops]
end
