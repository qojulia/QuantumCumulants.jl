"""
    Average <: SymbolicNumber

Symbolic number representing the average over an operator.
See also: [`average`](@ref)
"""
struct Average{T<:Number,OP} <: SymbolicNumber
    operator::OP
end
Average(operator::OP) where OP = Average{Number,OP}(operator)
Base.:(==)(a1::Average,a2::Average) = (a1.operator==a2.operator)
Base.hash(a::Average{T}, h::UInt) where T = hash(a.operator, hash(T, h))

Base.conj(a::Average) = Average(adjoint(a.operator))

"""
    average(::AbstractOperator)
    average(::AbstractOperator,order::Int)

Compute the average of an operator. If `order` is given, the [`cumulant_expansion`](@ref)
up to that order is computed immediately.
"""
average(op::BasicOperator) = Average(op)
function average(op::OperatorTerm)
    if op.f ∈ [+,-] # linearity
        avg = op.f(average.(op.arguments)...)
        return avg
    elseif op.f === (*)
        # Move constants out of average
        cs, ops = separate_constants(op)
        if isempty(cs)
            op_ = expand_trigs(op)
            op_ = expand(op_)
            if isequal(op_, op)
                return Average(op)
            else
                return average(op_)
            end
        else
            return op.f(cs...)*average(op.f(ops...))
        end
    elseif op.f === (^)
        arg, n = op.arguments
        op_ = *((arg for i=1:n)...)
        return average(op_)
    elseif op.f in [cos,sin]
        op_ = expand_trigs(op)
        return average(op_)
    else
        return Average(op)
    end
end
average(x::Number) = x

separate_constants(x::Number) = [x],[]
separate_constants(op::AbstractOperator) = [],[op]
function separate_constants(op::OperatorTerm{<:typeof(*)})
    cs = filter(x->isa(x,Number), op.arguments)
    ops = filter(x->isa(x,AbstractOperator), op.arguments)
    return cs, ops
end

function expand_trigs(op::OperatorTerm)
    if op.f === (cos)
        @assert length(op.arguments)==1
        arg = op.arguments[1]
        return 0.5*(exp(-im*arg) + exp(im*arg))
    elseif op.f === (sin)
        @assert length(op.arguments)==1
        arg = op.arguments[1]
        return 0.5im*(exp(-im*arg) + -1*exp(im*arg))
    else
        return op.f(expand_trigs.(op.arguments)...)
    end
end
expand_trigs(x) = x



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
    return DifferentialEquation(lhs,rhs)
end
average(arg,order;kwargs...) = cumulant_expansion(average(arg),order;kwargs...)

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
    elseif ord == Inf
        op = avg.operator
        out = _transcendent_cumulant_expansion(op,order)
        return simplify ? simplify_constants(out, kwargs...) : out
    else
        op = avg.operator
        @assert op.f === (*)
        if simplify
            return simplify_constants(_cumulant_expansion(op.arguments, order), kwargs...)
        else
            _cumulant_expansion(op.arguments, order)
        end
    end
end
function cumulant_expansion(avg::Average,order::Vector;mix_choice=maximum,kwargs...)
    aon = acts_on(avg.operator)
    order_ = mix_choice(order[i] for i in aon)
    return cumulant_expansion(avg,order_;kwargs...)
end
cumulant_expansion(x::Number,order;kwargs...) = x
function cumulant_expansion(x::NumberTerm,order;mix_choice=maximum, simplify=false, kwargs...)
    cumulants = [cumulant_expansion(arg,order;simplify=false,mix_choice=mix_choice) for arg in x.arguments]
    return simplify_constants(x.f(cumulants...);kwargs...)
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
    return DifferentialEquation(de.lhs,rhs)
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

    args_sum = Number[]
    for i=1:length(parts)
        p = collect(parts[i])
        for j=1:length(p) # Terms in the sum
            n = length(p[j])
            args_prod = Number[factorial(n-1)*(-1)^n]
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

# Compute the (joint) cumulant of a list of operators in terms of their raw moments
function cumulant(args::Vector)
    # Get all possible partitions
    parts = [partitions(args,i) for i=1:length(args)]

    args_sum = Number[]
    for i=1:length(parts)
        p = collect(parts[i])
        for j=1:length(p) # Terms in the sum
            n = length(p[j])
            args_prod = Number[factorial(n-1)*(-1)^(n-1)]
            for p_=p[j] # Product over partition blocks
                push!(args_prod, average(*(p_...)))
            end
            # Add terms in sum
            push!(args_sum, *(args_prod...))
        end
    end
    return +(args_sum...)
end

function _transcendent_cumulant_expansion(op::OperatorTerm,args...;kwargs...)
    if op.f === (*)
        return _transcendent_cumulant_expansion(op.arguments,args...;kwargs...)
    else
        error("Unknown function to expand, $(op.f)")
    end
end
function _transcendent_cumulant_expansion(op::OperatorTerm{<:typeof(exp)},order::Int;kwargs...)
    @assert length(op.arguments)==1
    arg = op.arguments[1]
    args = [arg]
    kn = Number[average(arg)]
    fac = 1.0
    for n=2:order
        fac /= n # factorial
        push!(args,arg)
        k = cumulant(args)
        push!(kn,fac*k)
    end
    return *((exp(k) for k in kn)...)
end
function _transcendent_cumulant_expansion(args::Vector,order::Int;kwargs...)
    inds = findall(x->(x isa OperatorTerm && x.f==(exp)), args)
    @assert length(inds)==1 # TODO: generalize
    idx = inds[1]
    expf = args[idx]
    @assert length(expf.arguments)==1
    x = expf.arguments[1]
    inds_comp = filter(!isequal(idx),1:length(args))
    inds2 = findall(x->(x isa AbstractOperator), args[inds_comp])
    if length(inds2) > 1
        error("Composite expression with arguments $args")
    end
    idx2 = inds2[1]
    y = args[inds_comp][idx2]
    c = commutator(y,x)
    (c isa Number) || error("Cumulant expansion not implemented for non-central operators $y and $x")

    exp_expanded = _transcendent_cumulant_expansion(expf,order;kwargs...)

    kn = Number[average(y)]
    fac = 1.0

    for n=2:order
        fac /= n
        xs = [x for i=1:n-1]
        ks_ = Number[]
        for k=1:n
            args_ = [xs[1:k-1];y;xs[1:n-k]]
            push!(ks_,cumulant(args_))
        end
        kn_ = +(ks_...)
        push!(kn, fac*kn_)
    end
    arg_expanded = +(kn...)

    if iszero(c)
        return arg_expanded*exp_expanded
    elseif idx2<idx
        return (arg_expanded + 0.5c)*exp_expanded
    else
        return (arg_expanded - 0.5c)*exp_expanded
    end
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
get_order(t::NumberTerm) = maximum(get_order.(t.arguments))
get_order(::Number) = 0
function get_order(t::OperatorTerm)
    if t.f in [+,-]
        return maximum(get_order.(t.arguments))
    elseif t.f === (*)
        return sum(get_order.(t.arguments))
    elseif t.f === (^)
        n = t.arguments[end]
        @assert n isa Integer
        return n
    elseif t.f in [cos,sin,exp] # TODO: add more functions
        if any(get_order.(t.arguments) .> 0) # Operator arguments
            return Inf
        else # Only number arguments
            return 0
        end
    end
    error("Couldn't determine order of expression $t")
end
get_order(::BasicOperator) = 1
