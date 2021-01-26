"""
    Average <: SymbolicNumber

Symbolic number representing the average over an operator.
See also: [`average`](@ref)
"""
struct Average{T<:Number,OP} <: SymbolicUtils.Symbolic{T}
    operator::OP
end
Average(operator::OP) where OP = Average{Number,OP}(operator)

Base.hash(a::Average{T}, h::UInt) where T = hash(a.operator, hash(T, h))
Base.isequal(a1::Average,a2::Average) = isequal(a1.operator,a2.operator)
Base.isone(::Average) = false
Base.iszero(::Average) = false

SymbolicUtils.:<ₑ(a::Average, b::Number) = false
SymbolicUtils.:<ₑ(a::Number,   b::Average) = true
SymbolicUtils.:<ₑ(a::Average,   b::Average) = SymbolicUtils.:<ₑ(a.operator, b.operator)
# SymbolicUtils.:<ₑ(a::Average, b::SymbolicUtils.Symbolic{<:Number}) = false
# SymbolicUtils.:<ₑ(a::SymbolicUtils.Symbolic{<:Number}, b::Average) = true

# Symbolic type promotion
for f in [+,-,*,/,^]
    @eval SymbolicUtils.promote_symtype(::$(typeof(f)),
                   T::Type{<:Average},
                   S::Type{<:Number}) = S
    @eval SymbolicUtils.promote_symtype(::$(typeof(f)),
                   T::Type{<:Number},
                   S::Type{<:Average}) = T
   @eval SymbolicUtils.promote_symtype(::$(typeof(f)),
                  T::Type{<:Average{T1}},
                  S::Type{<:Average{T2}}) where {T1,T2} = promote_type(T1,T2)
end

SymbolicUtils.symtype(x::Average{T}) where T = T
SymbolicUtils.to_symbolic(x::Average) = x

SymbolicUtils.istree(::Average) = false
# SymbolicUtils.islike(::Average, ::Type{<:Number}) = true


acts_on(avg::Average) = acts_on(avg.operator)

# SymbolicUtils.@number_methods(Average, SymbolicUtils.term(f, a), SymbolicUtils.term(f, a, b), skipbasics)
# Methods
# Base.conj(a::Average) = Average(adjoint(a.operator))
# Base.:+(a::Number, b::Average) = SymbolicUtils.Add(SymbolicUtils.add_t(a,b), SymbolicUtils.makeadd(1, a, b)...)
# Base.:+(a::Average, b::Number) = SymbolicUtils.Add(SymbolicUtils.add_t(a,b), SymbolicUtils.makeadd(1, b, a)...)
# Base.:+(a::Average, b::Average) = SymbolicUtils.Add(SymbolicUtils.add_t(a,b), SymbolicUtils.makeadd(1, 0, a, b)...)
# Base.:+(a::SymbolicUtils.Symbolic, b::Average) = SymbolicUtils.Add(SymbolicUtils.add_t(a,b), SymbolicUtils.makeadd(1, 0, a, b)...)
# Base.:+(a::Average, b::SymbolicUtils.Symbolic) = SymbolicUtils.Add(SymbolicUtils.add_t(a,b), SymbolicUtils.makeadd(1, 0, a, b)...)

# Base.:^(avg::Average,n::Integer) = SymbolicUtils.Pow(avg, n)
# Base.:^(avg::Average,n::SymbolicUtils.Symbolic{<:Number}) = SymbolicUtils.Pow(avg, n)
# Base.:^(avg::Average,n::Number) = SymbolicUtils.Pow(avg, n)

"""
    average(::AbstractOperator)
    average(::AbstractOperator,order::Int)

Compute the average of an operator. If `order` is given, the [`cumulant_expansion`](@ref)
up to that order is computed immediately.
"""
average(op::BasicOperator) = Average(op)
function average(op::OperatorTerm)
    f = SymbolicUtils.operation(op)
    if f ∈ [+,-] # linearity
        avg = f(average.(op.arguments)...)
        return avg
    elseif f === (*)
        # Move constants out of average
        cs, ops = separate_constants(op)
        if isempty(cs)
            op_ = expand(op)
            if isequal(op_, op)
                return Average(op)
            else
                return average(op_)
            end
        else
            return f(cs...)*average(f(ops...))
        end
    elseif f === (^)
        arg, n = op.arguments
        op_ = *((arg for i=1:n)...)
        return average(op_)
    else
        return Average(op)
    end
end
average(x::Union{T,SymbolicUtils.Symbolic{T}}) where T<:Number = x

separate_constants(x::Union{T,SymbolicUtils.Symbolic{T}}) where T<:Number = [x],[]
separate_constants(op::T) where T<:AbstractOperator = [],[op]
function separate_constants(op::OperatorTerm{<:typeof(*)})
    cs = filter(x->isa(x,Number)||isa(x,SymbolicUtils.Symbolic{<:Number}), op.arguments)
    ops = filter(x->isa(x,AbstractOperator), op.arguments)
    return cs, ops
end

"""
    average(::DifferentialEquation;multithread=false)
    average(::DifferentialEquation,order::Int;multithread=false)

Compute the average of a [`DifferentialEquation`](@ref) (or a set of equations).
Returns a [`DifferentialEquation`](@ref) with containing the corresponding equations
for averages. If `order` is specified, the [`cumulant_expansion`](@ref) up to
that order is computed immediately. The keyword `multithread` specifies whether
the averaging (and [`cumulant_expansion`](@ref)) should be parallelized (defaults to `false`).
"""
function average(de::DifferentialEquation;multithread=false)
    lhs = Vector{Number}(undef, length(de.lhs))
    rhs = Vector{Number}(undef, length(de.lhs))
    if multithread
        Threads.@threads for i=1:length(de.lhs)
            lhs[i] = average(de.lhs[i])
            rhs[i] = average(de.rhs[i])
        end
    else
        for i=1:length(de.lhs)
            lhs[i] = average(de.lhs[i])
            rhs[i] = average(de.rhs[i])
        end
    end
    return DifferentialEquation(lhs,rhs,de.hamiltonian,de.jumps,de.rates)
end
average(arg,order;kwargs...) = cumulant_expansion(average(arg),order;kwargs...)

# Type promotion -- average(::Operator)::Number
# SymbolicUtils.promote_symtype(average, ::Type{<:AbstractOperator}) = Number

# Cumulant expansion
"""
    cumulant_expansion(avg, order::Int)

For an [`Average`](@ref) of an operator, expand it in terms
of moments up to `order` neglecting their joint cumulant.

See also: https://en.wikipedia.org/wiki/Cumulant#Joint_cumulants

Examples
=======
```
julia> avg = average(a*b)
⟨a*b⟩

julia> cumulant_expansion(avg,1)
(⟨a⟩*⟨b⟩)

julia> avg = average(a*b*c)
⟨a*b*c⟩

julia> cumulant_expansion(avg,2)
((⟨a*b⟩*⟨c⟩)+(⟨a*c⟩*⟨b⟩)+(⟨a⟩*⟨b*c⟩)+(-2*⟨a⟩*⟨b⟩*⟨c⟩))
```

Optional arguments
=================
*simplify=true: Specify whether the result should be simplified.
*kwargs...: Further keyword arguments being passed to [`simplify_constants`](@ref)
"""
function cumulant_expansion(avg::Average,order::Int;simplify=true,kwargs...)
    @assert order > 0
    ord = get_order(avg)
    if ord <= order
        return avg
    else
        op = avg.operator
        if simplify
            avg_ = average(simplify_operators(op))
        else
            avg_ = average
        end
        if simplify && !isequal(avg, avg_) # TODO: better strategy to get proper ordering
            return cumulant_expansion(avg_,order;simplify=simplify,kwargs...)
        else
            @assert SymbolicUtils.operation(op) === (*)
            if simplify
                return simplify_operators(_cumulant_expansion(op.arguments, order), kwargs...)
            else
                _cumulant_expansion(op.arguments, order)
            end
        end
    end
end
function cumulant_expansion(avg::Average,order::Vector;mix_choice=maximum,kwargs...)
    aon = acts_on(avg.operator)
    order_ = mix_choice(order[i] for i in aon)
    return cumulant_expansion(avg,order_;kwargs...)
end
cumulant_expansion(x::Union{SymbolicUtils.Symbolic{T},T},order;kwargs...) where T<:Number = x
function cumulant_expansion(x::SymbolicUtils.Symbolic,order;mix_choice=maximum, simplify=false, kwargs...)
    if SymbolicUtils.istree(x)
        cumulants = [cumulant_expansion(arg,order;simplify=false,mix_choice=mix_choice) for arg in x.arguments]
        return simplify_operators(x.f(cumulants...);kwargs...)
    else
        error()
    end
end
function cumulant_expansion(de::DifferentialEquation{<:Number,<:Number},order;multithread=false,mix_choice=maximum,kwargs...)
    rhs = Vector{Number}(undef, length(de.lhs))
    if multithread
        Threads.@threads for i=1:length(de.lhs)
            check_lhs(de.lhs[i],order;mix_choice=mix_choice)
            cr = cumulant_expansion(de.rhs[i],order;mix_choice=mix_choice,kwargs...)
            rhs[i] = cr
        end
    else
        for i=1:length(de.lhs)
            check_lhs(de.lhs[i],order;mix_choice=mix_choice)
            cr = cumulant_expansion(de.rhs[i],order;mix_choice=mix_choice,kwargs...)
            rhs[i] = cr
        end
    end
    return DifferentialEquation(de.lhs,rhs,de.hamiltonian,de.jumps,de.rates)
end

function check_lhs(lhs,order::Int;kwargs...)
    (get_order(lhs) > order) && error("Cannot form cumulant expansion of derivative! Check the left-hand-side of your equations; you may want to use a higher order!")
    return nothing
end
function check_lhs(lhs,order::Vector;mix_choice=maximum)
    aon = acts_on(lhs.operator)
    order_ = mix_choice(order[i] for i in aon)
    (get_order(lhs) > order_) && error("Cannot form cumulant expansion of derivative! Check the left-hand-side of your equations; you may want to use a higher order!")
    return nothing
end

function _cumulant_expansion(args::Vector,order::Int)

    # Get all possible partitions; partitions(args,1) corresponds to the moment of order length(args)
    parts = [partitions(args,i) for i=2:length(args)]

    args_sum = Any[]
    for p in parts
        for pj in p
            n = length(pj)
            args_prod = Any[-factorial(n-1)*(-1)^(n-1)]
            for p_=pj # Product over partition blocks
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

"""
    cumulant(x,n=get_order(x);simplify=true,kwargs...)

Compute the `n`th cumulant of `x` (either an operator or an average).
The output is simplified when `simplify=true`. Further keyword arguments are
passed on to simplification.

Examples
========
```
julia> cumulant(a)
⟨a⟩

julia> cumulant(a*b)
(⟨a*b⟩+(-1*⟨a⟩*⟨b⟩))

julia> cumulant(a*b,1)
⟨a*b⟩

julia> cumulant(a*b,3)
0
```
"""
function cumulant(op::OperatorTerm,n::Int=get_order(op);simplify=true,kwargs...)
    order = get_order(op)
    if order < n
        return zero(op)
    else
        if simplify
            avg_ = average(simplify_operators(op))
        else
            avg_ = average(op)
        end
        if simplify && !isequal(average(op), avg_) # TODO: better strategy to get proper ordering
            return cumulant(avg_.operator,order;simplify=simplify,kwargs...)
        else
            @assert SymbolicUtils.operation(op) === (*)
            if simplify
                return simplify_operators(_cumulant(op.arguments, n), kwargs...)
            else
                return _cumulant(op.arguments, n)
            end
        end
    end
end
cumulant(avg::Average,args...;kwargs...) = cumulant(avg.operator,args...;kwargs...)
cumulant(op::BasicOperator,n::Int=1;kwargs...) = isone(n) ? average(op) : zero(op)

function _cumulant(args::Vector,m::Int=length(args))
    parts = [partitions(args,i) for i=1:m]
    args_sum = Any[]
    for i=1:length(parts)
        p = collect(parts[i])
        for j=1:length(p) # Terms in the sum
            n = length(p[j])
            args_prod = Any[factorial(n-1)*(-1)^(n-1)]
            for p_=p[j] # Product over partition blocks
                push!(args_prod, Average(*(p_...)))
            end
            # Add terms in sum
            push!(args_sum, *(args_prod...))
        end
    end
    return average(+(args_sum...))
end

"""
    get_order(arg)

Compute the order of a given argument. This is the order used to decide whether
something should be expanded using a [`cumulant_expansion`](@ref) method.

Examples
=======
```
julia> get_order(a)
1

julia> get_order(a*b)
2

julia> get_order(1)
0
```
"""
get_order(avg::Average) = get_order(avg.operator)
function get_order(t::SymbolicUtils.Symbolic)
    if SymbolicUtils.istree(t)
        return maximum(map(get_order, SymbolicUtils.arguments(t)))
    else
        return 0
    end
end
get_order(::Number) = 0
function get_order(t::OperatorTerm)
    if t.f in [+,-]
        return maximum(get_order.(t.arguments))
    elseif t.f === (*)
        return length(t.arguments)
    elseif t.f === (^)
        n = t.arguments[end]
        @assert n isa Integer
        return n
    end
    error("Unknown function $(t.f)")
end
get_order(::BasicOperator) = 1
