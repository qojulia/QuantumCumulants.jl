"""
    find_missing(rhs::Vector, vs::Vector, vs_adj=adjoint.(vs) ps=[])

For a list of expressions contained in `rhs`, check whether all occurring symbols
are contained either in the variables given in `vs` or `ps`. Returns a list of
missing symbols.
"""
function find_missing(rhs::Vector{<:Number}, vs::Vector{<:Number}; vs_adj::Vector=adjoint.(vs), ps=[])
    missed = Number[]
    for e=rhs
        append!(missed,get_symbolics(e))
    end
    unique!(missed)
    filter!(x->!(x∈vs || x∈ps || x∈vs_adj),missed)
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
    complete(de::DifferentialEquation,H)
    complete(de::DifferentialEquation,H,J)

From a set of differential equation of averages, find all averages that are missing
and derive the corresponding equations of motion.
"""
function complete(de::DifferentialEquation{<:Number,<:Number},H::AbstractOperator;kwargs...)
    rhs_, lhs_ = complete(de.rhs,de.lhs, [H], [];kwargs...)
    return DifferentialEquation(lhs_,rhs_)
end
function complete(de::DifferentialEquation{<:Number,<:Number},H::AbstractOperator,J::Vector;kwargs...)
    rhs_, lhs_ = complete(de.rhs,de.lhs, [H,J]; kwargs...)
    return DifferentialEquation(lhs_,rhs_)
end
function complete(rhs::Vector{<:Number}, vs::Vector{<:Number}, he_args; order=nothing, mix_choice=maximum, kwargs...)
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
    missed = unique_ops(find_missing(rhs_, vs_))
    filter!(x->isa(x,Average),missed)
    while !isempty(missed)
        ops = getfield.(missed, :operator)
        he = heisenberg(ops,he_args...;kwargs...)
        he_avg = average(he,order_;mix_choice=mix_choice)
        rhs_ = [rhs_;he_avg.rhs]
        vs_ = [vs_;he_avg.lhs]
        missed = unique_ops(find_missing(rhs_,vs_))
        filter!(x->isa(x,Average),missed)
    end
    return rhs_, vs_
end

"""
    find_operators(::HilbertSpace, order; names=nothing)

Find all operators that fully define a system up to the given `order`.
"""
function find_operators(h::ProductSpace, order::Int; names=nothing, kwargs...)
    if names isa Nothing && (unique(typeof.(h.spaces))!=typeof.(h.spaces))
        alph = 'a':'z'
        names_ = Symbol.(alph[1:length(h.spaces)])
    else
        names_ = names
    end
    fund_ops = fundamental_operators(h;names=names_, kwargs...)
    fund_ops = unique([fund_ops;adjoint.(fund_ops)])
    return _find_operators(fund_ops,order;kwargs...)
end
function find_operators(h::HilbertSpace,order::Int;names=nothing,kwargs...)
    fund_ops = fundamental_operators(h;names=names)
    return _find_operators(fund_ops,order;kwargs...)
end
find_operators(op::AbstractOperator,args...;kwargs...) = find_operators(hilbert(op),args...;kwargs...)

function _find_operators(fund_ops,order; kwargs...)
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
function fundamental_operators(h::PositionSpace,aon::Int=1;names=nothing)
    names_ = names isa Nothing ? [:x,:p] : names[aon]
    x = Position(h,names_[1],aon)
    p = Momentum(h,names_[2],aon)
    return [x,p]
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
