function average end

"""
    AvgSym <: CNumber

Symbolic number representing the average over an operator.
See also: [`average`](@ref)
"""
struct AvgSym <: CNumber end

const Average = SymbolicUtils.Term{<:AvgSym}

# Direct construction of average symbolic expression
function _average(operator)
    return SymbolicUtils.Term{AvgSym}(average, [operator])
end

# Type promotion -- average(::Operator)::Number
SymbolicUtils.promote_symtype(average, ::Type{<:QNumber}) = AvgSym

function acts_on(s::SymbolicUtils.Symbolic)
    if SymbolicUtils.istree(s)
        f = SymbolicUtils.operation(s)
        if f === average
            return acts_on(SymbolicUtils.arguments(s)[1])
        else
            aon = Int[]
            for arg in SymbolicUtils.arguments(s)
                append!(aon, acts_on(arg))
            end
            unique!(aon)
            sort!(aon)
            return aon
        end
    else
        return Int[]
    end
end

"""
    average(::QNumber)
    average(::QNumber,order::Int)

Compute the average of an operator. If `order` is given, the [`cumulant_expansion`](@ref)
up to that order is computed immediately.
"""
average(op::QSym) = _average(op)
function average(op::QTerm)
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
                return _average(op)
            else
                return average(op_)
            end
        else
            return f(cs...)*average(f(ops...))
        end
    elseif f === (^)
        arg, n = op.arguments
        op_ = QTerm(*, [arg for i=1:n])
        return average(op_)
    else
        return _average(op)
    end
end
average(x::Union{T,SymbolicUtils.Symbolic{T}}) where T<:Number = x

separate_constants(x::Union{T,SymbolicUtils.Symbolic{T}}) where T<:Number = [x],[]
separate_constants(op::T) where T<:QNumber = [],[op]
function separate_constants(op::QTerm{<:typeof(*)})
    cs = filter(x->isa(x,Number)||isa(x,SymbolicUtils.Symbolic{<:Number}), op.arguments)
    ops = filter(x->isa(x,QNumber), op.arguments)
    return cs, ops
end

"""
    average(::HeisenbergEquation;multithread=false)
    average(::HeisenbergEquation,order::Int;multithread=false)

Compute the average of a [`HeisenbergEquation`](@ref) (or a set of equations).
Returns a [`HeisenbergEquation`](@ref) with containing the corresponding equations
for averages. If `order` is specified, the [`cumulant_expansion`](@ref) up to
that order is computed immediately. The keyword `multithread` specifies whether
the averaging (and [`cumulant_expansion`](@ref)) should be parallelized (defaults to `false`).
"""
function average(de::HeisenbergEquation;multithread=false)
    lhs = Vector{Any}(undef, length(de.lhs))
    rhs = Vector{Any}(undef, length(de.lhs))
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
    return HeisenbergEquation(lhs,rhs,de.hamiltonian,de.jumps,de.rates)
end
average(arg,order;kwargs...) = cumulant_expansion(average(arg),order;kwargs...)


# Cumulant expansion
"""
    cumulant_expansion(avg, order::Int)

For an [`AvgSym`](@ref) of an operator, expand it in terms
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
*kwargs...: Further keyword arguments being passed to [`qsimplify`](@ref)
"""
function cumulant_expansion(avg::SymbolicUtils.Term{<:AvgSym},order::Int;simplify_input=true,simplify_output=true,kwargs...)
    @assert order > 0
    ord = get_order(avg)
    if ord <= order
        return avg
    else
        op = SymbolicUtils.arguments(avg)[1]
        if simplify_input
            avg_ = average(qsimplify(op))
        else
            avg_ = avg
        end
        if simplify_input && !isequal(avg, avg_)
            return cumulant_expansion(avg_,order;simplify_input=simplify_input,simplify_output=simplify_output,kwargs...)
        else
            @assert SymbolicUtils.operation(op) === (*)
            if simplify_output
                return qsimplify(_cumulant_expansion(op.arguments, order); kwargs...)
            else
                _cumulant_expansion(op.arguments, order)
            end
        end
    end
end
function cumulant_expansion(avg::SymbolicUtils.Term{<:AvgSym},order::Vector;mix_choice=maximum,kwargs...)
    aon = acts_on(avg)
    order_ = mix_choice(order[i] for i in aon)
    return cumulant_expansion(avg,order_;kwargs...)
end
cumulant_expansion(x::Number,order;kwargs...) = x
function cumulant_expansion(x::SymbolicUtils.Symbolic,order;mix_choice=maximum,simplify_input=false,simplify_output=true,kwargs...)
    if SymbolicUtils.istree(x)
        if order isa Int
            get_order(x) <= order && return x
        end
        f = SymbolicUtils.operation(x)
        args = SymbolicUtils.arguments(x)
        cumulants = [cumulant_expansion(arg,order;simplify_input=simplify_input,simplify_output=false,mix_choice=mix_choice) for arg in args]
        if simplify_output
            return qsimplify(f(cumulants...);kwargs...)
        else
            return f(cumulants...)
        end
    else
        return x
    end
end
function cumulant_expansion(de::HeisenbergEquation,order;multithread=false,mix_choice=maximum,kwargs...)
    rhs = Vector{Any}(undef, length(de.lhs))
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
    return HeisenbergEquation(de.lhs,rhs,de.hamiltonian,de.jumps,de.rates)
end

function check_lhs(lhs,order::Int;kwargs...)
    (get_order(lhs) > order) && error("Cannot form cumulant expansion of derivative! Check the left-hand-side of your equations; you may want to use a higher order!")
    return nothing
end
function check_lhs(lhs,order::Vector;mix_choice=maximum)
    aon = acts_on(lhs)
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
                    if length(p_)==1
                        op_ = p_[1]
                    else
                        op_ = QTerm(*, p_)
                    end
                    push!(args_prod, _average(op_))
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
function cumulant(op::QTerm,n::Int=get_order(op);simplify=true,kwargs...)
    order = get_order(op)
    if order < n
        return zero(op)
    else
        if simplify
            avg_ = average(qsimplify(op))
        else
            avg_ = average(op)
        end
        if simplify && !isequal(average(op), avg_) # TODO: better strategy to get proper ordering
            return cumulant(avg_.operator,order;simplify=simplify,kwargs...)
        else
            @assert SymbolicUtils.operation(op) === (*)
            if simplify
                return qsimplify(_cumulant(op.arguments, n), kwargs...)
            else
                return _cumulant(op.arguments, n)
            end
        end
    end
end
cumulant(avg::SymbolicUtils.Term{<:AvgSym},args...;kwargs...) = cumulant(SymbolicUtils.arguments(avg)[1],args...;kwargs...)
cumulant(op::QSym,n::Int=1;kwargs...) = isone(n) ? average(op) : zero(op)

function _cumulant(args::Vector,m::Int=length(args))
    parts = [partitions(args,i) for i=1:m]
    args_sum = Any[]
    for i=1:length(parts)
        p = collect(parts[i])
        for j=1:length(p) # Terms in the sum
            n = length(p[j])
            args_prod = Any[factorial(n-1)*(-1)^(n-1)]
            for p_=p[j] # Product over partition blocks
                push!(args_prod, _average(*(p_...)))
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
get_order(avg::SymbolicUtils.Term{<:AvgSym}) = get_order(SymbolicUtils.arguments(avg)[1])
function get_order(t::SymbolicUtils.Symbolic)
    if SymbolicUtils.istree(t)
        return maximum(map(get_order, SymbolicUtils.arguments(t)))
    else
        return 0
    end
end
get_order(::Number) = 0
function get_order(t::QTerm)
    if t.f in [+,-]
        order = Int[get_order(arg) for arg in t.arguments]
        return maximum(order)
    elseif t.f === (*)
        order = Int[get_order(arg) for arg in t.arguments]
        return sum(order)
    elseif t.f === (^)
        n = t.arguments[end]::Int
        return n
    end
    error("Unknown function $(t.f)")
end
get_order(::QSym) = 1
