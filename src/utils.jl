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
    vs_ = collect(Iterators.flatten(get_symbolics.(vs)))
    # vs_ = [(get_symbolics.(vs)...)...]
    filter!(x->isa(x,Average),vs_)
    vars = [vs_;adjoint.(vs_);ps]
    unique!(vars)
    if isempty(ps)
        filter!(x->!isa(x,Parameter), missed)
        filter!(x->!isa(x,IndexedParameter), missed)
    end
    filter!(x->!(x∈vars),missed)

    if any(.!isempty.(find_index.(missed)))
        filter!(x->!isa(x,Index), missed) #filter Indices (not indexed objects)
        missed = _filter_indexed(missed, vars)
    end
    return missed
end
function find_missing(de::DifferentialEquation{<:Number,<:Number}; kwargs...)
    find_missing(de.rhs, de.lhs; kwargs...)
end

function _filter_indexed(missed_, vars)
    missed = Number[]
    vars_and_already_missed = vars
    for m in missed_
        _in_without_index(m, vars_and_already_missed) || (push!(missed, m); push!(vars_and_already_missed, [m, adjoint(m)]...))
    end
    return missed
end

# Check if a symbolic with index is contained in a collection disregarding the respective indices
function _in_without_index(x,itr)
    idx = find_index(x)
    isempty(idx) && return false
    for i in itr
        _isequal_without_index(x, idx, i) && return true
    end
    return false
end
_isequal_without_index(x,idx_x,y) = false

# Compare two indexed objects disregarding their indices
function _isequal_without_index(x::T,idx_x,y::T) where T
    idx_y = find_index(y)
    _check_idx(idx_x, idx_y) || return false
    x_ = _construct_without_index(x)
    y_ = _construct_without_index(y)
    return isequal(x_,y_)
end

# Construct an Indexed object without index
_construct_without_index(p::IndexedParameter{T}) where T = Parameter{T}(p.name)
for T = [:Create,:Destroy]
    Name = Symbol(:Indexed,T)
    @eval _construct_without_index(op::$(Name)) = $(T)(op.hilbert,op.name,op.aon)
end
_construct_without_index(op::IndexedTransition) = Transition(op.hilbert,op.name,op.i,op.j,op.aon)
_construct_without_index(avg::Average) = Average(_construct_without_index(avg.operator))
function _construct_without_index(op::OperatorTerm)
    args = []
    for arg in op.arguments
        isempty(find_index(arg)) ? push!(args, arg) : push!(args, _construct_without_index(arg))
    end
    sort!(args, by=acts_on)
    return op.f(args...)
end

# Check whether two indices have the same properties besides their names
function _check_idx(idx1::Vector,idx2::Vector)
    length(idx1)==length(idx2) || return false
    for (i1,i2) in zip(idx1,idx2)
        _check_idx(i1,i2) || return false
    end
    return true
end
_check_idx(idx1::Index,idx2::Index) = isequal(idx1.count, idx2.count)

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

function help_f_missed_indexed(missed, LHS_idx_list) #just to avoid writing it twice
    for it = 1:length(missed)
        m = missed[it]
        m_idx = find_index(m)
        if !isempty(m_idx)
            LHS_idx_list_ = copy(LHS_idx_list) #needed to delete already use lhs indices
            new_m_idx = []
            for it_m = 1:length(m_idx)
                i_lhs = LHS_idx_list_[findfirst(x->_check_idx(x,m_idx[it_m]), LHS_idx_list_)] #corresponding index on the lhs
                push!(new_m_idx, i_lhs)
                filter!(x->!isequal(x,i_lhs), LHS_idx_list_)
            end
            missed[it] = swap_index(m, m_idx, new_m_idx)
        end
    end
end

function complete(rhs::Vector{<:Number}, vs::Vector{<:Number}, H, J, rates; order=nothing, mix_choice=maximum, LHS_idx_list=[], kwargs...)
    order_lhs = maximum(get_order.(vs))
    order_rhs = maximum(get_order.(rhs))
    if order isa Nothing
        order_ = max(order_lhs, order_rhs)
    else
        order_ = order
    end
    order_ >= order_lhs || error("Cannot form cumulant expansion of derivative; you may want to use a higher order!")

    vs_ = copy(vs)
    rhs_ = [cumulant_expansion(r, order_) for r in rhs]
    missed = find_missing(rhs_, vs_)
    filter!(x->isa(x,Average),missed) #keep only averages

    if any(has_indexed.(missed))
        isempty(LHS_idx_list) && error("LHS_idx_list is needed!")
        sub_ij0 = Dict()
        if length(LHS_idx_list) > 1 # for substituting i==j in nips on the LHS
            ijs = collect(combinations(LHS_idx_list, 2))
            for ij in ijs
                sub_ij0[ij[1]==ij[2]] = 0
                sub_ij0[ij[2]==ij[1]] = 0
            end
        end
        help_f_missed_indexed(missed,LHS_idx_list)
    end
    missed = unique_ops(missed)

    while !isempty(missed)
        ops = getfield.(missed, :operator)
        ops = simplify_operators.(ops)
        if length(LHS_idx_list) > 1
            ops = [substitute(op, sub_ij0) for op in ops]
        end
        he = isempty(J) ? heisenberg(ops,H) : heisenberg(ops,H,J;rates=rates)
        he_avg = average(he,order_;mix_choice=mix_choice)
        rhs_ = [rhs_;he_avg.rhs]
        vs_ = [vs_;he_avg.lhs]
        missed = find_missing(rhs_,vs_)
        filter!(x->isa(x,Average), missed)
        if any(has_indexed.(missed))
            help_f_missed_indexed(missed,LHS_idx_list)
        end
        missed = unique_ops(missed)
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

has_expr(f, x) = false
function has_expr(f::Function, x::Union{OperatorTerm,NumberTerm})
    (f===x.f) && return true
    for arg in x.arguments
        has_expr(f, arg) && return true
    end
    return false
end
function has_expr(y::Union{OperatorTerm,NumberTerm}, x::Union{OperatorTerm,NumberTerm})
    isequal(y, x) && return true
    for arg in x.arguments
        has_expr(y, arg) && return true
    end
    return false
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
function _to_expression(t::Union{OperatorTerm,NumberTerm})
    if t.f === Sum
        return :( Sum($(t.arguments[1]), $(t.arguments[2:end])) )
    else
        return :( $(Symbol(t.f))($(_to_expression.(t.arguments)...)) )
    end
end
_to_expression(p::Parameter) = p.name
function _to_expression(avg::Average)
    ex = _to_expression(avg.operator)
    return :(AVERAGE($ex))
end

_to_expression(p::IndexedParameter) = :(IndexedParameter($(p.name), $(_to_expression.(p.index))))
_to_expression(op::IndexedDestroy) = :(Indexed($(op.name), $(_to_expression(op.index))))
_to_expression(op::IndexedCreate) = :(Indexed(dagger($(op.name)), $(_to_expression(op.index))))
_to_expression(op::IndexedTransition) = :(Indexed(Transition($(op.name),$(op.i),$(op.j)), $(_to_expression(op.index))))
_to_expression(op::OperatorTerm{<:typeof(nip)}) = _to_expression(*(op.arguments...))
_to_expression(i::Index) = i.name
