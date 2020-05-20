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
function find_missing(de::DifferentialEquationSet{<:Number,<:Number}; kwargs...)
    find_missing(de.rhs, de.lhs; kwargs...)
end
function find_missing(de::DifferentialEquation{<:Number,<:Number}; kwargs...)
    find_missing([de.rhs], [de.lhs]; kwargs...)
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
    complete(de::DifferentialEquationSet,order,H)
    complete(de::DifferentialEquationSet,order,H,J)

From a set of differential equation of averages, find all averages that are missing
and derive the corresponding equations of motion.
"""
function complete(de::DifferentialEquationSet{<:Number,<:Number},order::Int,H::AbstractOperator)
    rhs_, lhs_ = complete(de.rhs,de.lhs,order, [H], [])
    return DifferentialEquationSet(lhs_,rhs_)
end
function complete(de::DifferentialEquationSet{<:Number,<:Number},order::Int,H::AbstractOperator,J::Vector;kwargs...)
    rhs_, lhs_ = complete(de.rhs,de.lhs,order, [H,J]; kwargs...)
    return DifferentialEquationSet(lhs_,rhs_)
end
function complete(rhs::Vector{<:Number}, vs::Vector{<:Number}, order::Int, he_args; kwargs...)
    rhs_ = copy(rhs)
    vs_ = copy(vs)
    @assert !any(get_order.(vs_) .> order)
    missed = unique_ops(find_missing(rhs_, vs_))
    filter!(x->isa(x,Average),missed)
    while !isempty(missed)
        ops = getfield.(missed, :operator)
        he = heisenberg(ops,he_args...;kwargs...)
        he_avg = average(he,order)
        rhs_ = [rhs_;he_avg.rhs]
        vs_ = [vs_;he_avg.lhs]
        missed = unique_ops(find_missing(rhs_,vs_))
        filter!(x->isa(x,Average),missed)
    end
    return rhs_, vs_
end


function find_operators(h::HilbertSpace, order::Int)
    fund_ops = fundamental_operators(h)
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

    ops_simp = simplify_operators.(all_ops)
    ops_only = AbstractOperator[]
    for op in ops_simp
        append!(ops_only, get_operators(op))
    end
    return unique_ops(ops_only)
end
find_operators(op::AbstractOperator,args...) = find_operators(hilbert(op),args...)

hilbert(op::BasicOperator) = op.hilbert
hilbert(t::OperatorTerm) = hilbert(t.arguments[findfirst(x->isa(x,AbstractOperator), t.arguments)])

function fundamental_operators(h::FockSpace,aon::Int=1)
    a = Destroy(h,Symbol(:a,SubscriptMap[aon]))
    return [a]
end
function fundamental_operators(h::NLevelSpace,aon::Int=1)
    sigmas = Transition[]
    lvls = levels(h)
    for i=1:length(lvls)
        for j=i:length(lvls)
            (i==j) && lvls[i]==ground_state(h) && continue
            s = Transition(h,Symbol(:σ,SubscriptMap[aon]),lvls[i],lvls[j])
            push!(sigmas,s)
        end
    end
    return sigmas
end
function fundamental_operators(h::ProductSpace)
    ops = []
    for i=1:length(h.spaces)
        ops_ = fundamental_operators(h.spaces[i],i)
        ops_ = [embed(h,o,i) for o in ops_]
        append!(ops,ops_)
    end
    return ops
end

get_operators(::Number) = []
get_operators(op::BasicOperator) = [op]
function get_operators(op::OperatorTerm{<:typeof(*)})
    args = AbstractOperator[]
    for arg in op.arguments
        append!(args, get_operators(arg))
    end
    isempty(args) && return args
    return [*(args...)]
end
function get_operators(t::OperatorTerm)
    ops = AbstractOperator[]
    for arg in t.arguments
        append!(ops, get_operators(arg))
    end
    return ops
end

function unique_ops(ops)
    seen = eltype(ops)[]
    for op in ops
        if !(op in seen || op' in seen)
            push!(seen, op)
        end
    end
    return seen
end

const SubscriptMap = Dict{Int,Char}(
            0 => '₀',
            1 => '₁',
            2 => '₂',
            3 => '₃',
            4 => '₄',
            5 => '₅',
            6 => '₆',
            7 => '₇',
            8 => '₈',
            9 => '₉')
