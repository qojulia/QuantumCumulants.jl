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
function complete(rhs::Vector{<:Number}, vs::Vector{<:Number}, he_args; order=maximum(get_order.(vs)), mix_choice=maximum, kwargs...)
    rhs_ = copy(rhs)
    vs_ = copy(vs)
    missed = unique_ops(find_missing(rhs_, vs_))
    filter!(x->isa(x,Average),missed)
    while !isempty(missed)
        ops = getfield.(missed, :operator)
        he = heisenberg(ops,he_args...;kwargs...)
        he_avg = average(he,order;mix_choice=mix_choice)
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

"""
    swap_idx(x, i::Index, j::Index)

Swap all indices `i` in the expression `x` with `j`. If the index does not match,
the expression is left unchanged.
"""
function swap_idx(x::T, i::Index, j::Index) where T<:Union{BasicOperator,Parameter}
    if get_index(x)==i
        fields = [getfield(x, name) for name in fieldnames(T)]
        fields[end] = j
        return T(fields...)
    else
        return x
    end
end
function swap_idx(t::Union{OperatorTerm,NumberTerm}, i::Index, j::Index)
    args = Any[]
    for arg in t.arguments
        push!(args, swap_idx(arg, i, j))
    end
    return t.f(args...)
end
swap_idx(x::Number, args...) = x



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
function _to_expression(op::BasicOperator)
    n = op.name
    idx = get_index(op)
    if idx != default_index()
        return Expr(:ref, n, idx.i)
    else
        return n
    end
end
function _to_expression(op::Create)
    n = :(dagger($(op.name)))
    idx = get_index(op)
    if idx != default_index()
        return Expr(:ref, n, idx.i)
    else
        return n
    end
end
_to_expression(op::Transition) = :( Transition($(op.name),$(op.i),$(op.j),$(op.index)) )
_to_expression(t::Union{OperatorTerm,NumberTerm}) = :( $(Symbol(t.f))($(_to_expression.(t.arguments)...)) )
_to_expression(p::Parameter) = p.name
function _to_expression(avg::Average)
    ex = _to_expression(avg.operator)
    return :(AVERAGE($ex))
end
