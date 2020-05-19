using Combinatorics: partitions, combinations

struct Average{T<:Number,OP} <: SymbolicNumber
    operator::OP
end
Average(operator::OP) where OP = Average{Number,OP}(operator)
Base.:(==)(a1::Average,a2::Average) = (a1.operator==a2.operator)

Base.conj(a::Average) = Average(adjoint(a.operator))

average(op::BasicOperator) = Average(op)
function average(op::OperatorTerm)
    if op.f ∈ [+,-] # linearity
        avg = op.f(average.(op.arguments)...)
        return avg
    elseif op.f === (*)
        # Move constants out of average
        cs, ops = separate_constants(op)
        if isempty(cs)
            return Average(op)
        else
            return op.f(cs...)*Average(op.f(ops...))
        end
    else
        return Average(op)
    end
end
average(x::Number) = x

separate_constants(x::Number) = [x],[]
separate_constants(op::BasicOperator) = [],[op]
function separate_constants(op::OperatorTerm)
    cs = []
    ops = []
    for arg in op.arguments
        c_, op_ = separate_constants(arg)
        append!(cs, c_)
        append!(ops, op_)
    end
    return cs, ops
end

function average(de::DifferentialEquation)
    lhs = average(de.lhs)
    rhs = average(de.rhs)
    return DifferentialEquation(lhs,rhs)
end
function average(de::DifferentialEquationSet)
    lhs = average.(de.lhs)
    rhs = average.(de.rhs)
    return DifferentialEquationSet(lhs,rhs)
end
average(arg,order::Int) = cumulant_expansion(average(arg),order)

# Conversion to SymbolicUtils
_to_symbolic(a::Average{T}) where T<:Number = SymbolicUtils.term(average, _to_symbolic(a.operator); type=T)
function _to_qumulants(t::SymbolicUtils.Term{T}) where T<:Number
    if t.f===average
        return Average(_to_qumulants(t.arguments[1]))
    else
        return NumberTerm{T}(t.f, _to_qumulants.(t.arguments))
    end
end

# Cumulant expansion
"""
    cumulant_expansion(avg, order::Int)

For a product of operators whose constituents are given in `args`, expand it in terms
of moments up to `order` neglecting their joint cumulant.

See also: https://en.wikipedia.org/wiki/Cumulant#Joint_cumulants
"""
function cumulant_expansion(avg::Average,order::Int;simplify=true)
    @assert order > 0
    ord = get_order(avg)
    if ord <= order
        return avg
    else
        op = avg.operator
        @assert op.f === (*)
        if simplify
            return simplify_constants(_cumulant_expansion(op.arguments, order))
        else
            _cumulant_expansion(op.arguments, order)
        end
    end
end
cumulant_expansion(x::Number,order;kwargs...) = x
function cumulant_expansion(x::NumberTerm,order;kwargs...)
    cumulants = [cumulant_expansion(arg,order;simplify=false) for arg in x.arguments]
    return simplify_constants(x.f(cumulants...);kwargs...)
end
function cumulant_expansion(de::DifferentialEquation,order;kwargs...)
    return DifferentialEquation(cumulant_expansion(de.lhs,order;kwargs...),cumulant_expansion(de.rhs,order;kwargs...))
end
function cumulant_expansion(de::DifferentialEquationSet,order;kwargs...)
    lhs = Number[]
    rhs = Number[]
    for i=1:length(de.lhs)
        cl = cumulant_expansion(de.lhs[i],order;kwargs...)
        cr = cumulant_expansion(de.rhs[i],order;kwargs...)
        push!(lhs, cl)
        push!(rhs, cr)
    end
    return DifferentialEquationSet(lhs,rhs)
end

function _cumulant_expansion(args::Vector,order::Int)

    # Get all possible partitions; partitions(args,1) corresponds to the moment of order length(args)
    parts = [partitions(args,i) for i=2:length(args)]

    args_sum = Number[]
    for i=1:length(parts)
        p = collect(parts[i])
        for j=1:length(p) # Terms in the sum
            n = length(p[j])
            args_prod = Number[-factorial(n-1)*(-1)^(n-1)]
            for p_=p[j] # Product over partition blocks
                if length(p_) > order # If the encountered moment is larger than order, apply expansion
                    push!(args_prod, _cumulant_expansion(p_, order))
                else # Else, average and add its product
                    push!(args_prod, Average(*(p_...)))
                end
            end
            # Add terms in sum
            push!(args_sum, *(args_prod...))
        end
    end
    return average(+(args_sum...))
end

get_order(avg::Average) = get_order(avg.operator)
get_order(t::NumberTerm) = maximum(get_order.(t.arguments))
get_order(::Number) = 0
function get_order(t::OperatorTerm)
    if t.f in [+,-]
        return maximum(get_order.(t.arguments))
    else
        return length(t.arguments)
    end
end
get_order(::BasicOperator) = 1
