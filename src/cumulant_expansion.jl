## Cumulant expansion

"""
    cumulant_expansion(avg, order::Int)

For an [`average`](@ref) of an operator, expand it in terms
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
*kwargs...: Further keyword arguments being passed to simplification.
"""
function cumulant_expansion(
    x::SymbolicUtils.Symbolic,
    order::Integer;
    simplify = true,
    kwargs...,
)
    if SymbolicUtils.iscall(x)
        get_order(x) <= order && return x
        f = SymbolicUtils.operation(x)
        if f===sym_average
            op = SymbolicUtils.arguments(x)[1]
            return _cumulant_expansion(op.args_nc, order)
        else
            args = SymbolicUtils.arguments(x)
            cumulants = [cumulant_expansion(arg, order; kwargs...) for arg in args]
            return f(cumulants...)
        end
    elseif x isa AvgSums
        return _cumulant_expansion(x, order; simplify, kwargs...) # basically just another cumulant_expansion dispatch
    else
        return x
    end
end
function cumulant_expansion(avg::Average, order::Vector; mix_choice = maximum, kwargs...)
    aon = acts_on(avg)
    order_ = mix_choice(order[get_i(i)] for i in aon)
    return cumulant_expansion(avg, order_; kwargs...)
end
cumulant_expansion(x::Number, order; kwargs...) = x

"""
    average(::QNumber, order)

Compute the average of an operator. If `order` is given, the [`cumulant_expansion`](@ref)
up to that order is computed immediately.
"""
SQA.average(x, order; kwargs...) = cumulant_expansion(average(x), order; kwargs...)

function cumulant_expansion(
    x::SymbolicUtils.Symbolic,
    order;
    mix_choice = maximum,
    simplify = true,
    kwargs...,
)
    if SymbolicUtils.iscall(x)
        f = SymbolicUtils.operation(x)
        args = SymbolicUtils.arguments(x)
        cumulants = [
            cumulant_expansion(arg, order; simplify = simplify, mix_choice = mix_choice) for arg in args
        ]
        if simplify
            return SymbolicUtils.simplify(f(cumulants...); kwargs...)
        else
            return f(cumulants...)
        end
    elseif x isa AvgSums
        return _cumulant_expansion(x, order; mix_choice, simplify, kwargs...) # basically just another cumulant_expansion dispatch
    else
        return x
    end
end
function cumulant_expansion(
    de::AbstractMeanfieldEquations,
    order;
    multithread = false,
    mix_choice = maximum,
    kwargs...,
)
    order==de.order && return de
    eqs = de.equations
    eqs_out = Vector{Symbolics.Equation}(undef, length(eqs))
    if multithread
        Threads.@threads for i = 1:length(eqs)
            cr = cumulant_expansion(eqs[i].rhs, order; mix_choice = mix_choice, kwargs...)
            eqs_out[i] = Symbolics.Equation(eqs[i].lhs, cr)
        end
    else
        for i = 1:length(eqs)
            cr = cumulant_expansion(eqs[i].rhs, order; mix_choice = mix_choice, kwargs...)
            eqs_out[i] = Symbolics.Equation(eqs[i].lhs, cr)
        end
    end
    return MeanfieldEquations(
        eqs_out,
        de.operator_equations,
        de.states,
        de.operators,
        de.hamiltonian,
        de.jumps,
        de.jumps_dagger,
        de.rates,
        de.iv,
        de.varmap,
        order,
    )
end
function cumulant_expansion(
    de::ScaledMeanfieldEquations,
    order;
    multithread = false,
    mix_choice = maximum,
    kwargs...,
)
    order==de.order && return de
    eqs = de.equations
    eqs_out = Vector{Symbolics.Equation}(undef, length(eqs))
    if multithread
        Threads.@threads for i = 1:length(eqs)
            cr = cumulant_expansion(eqs[i].rhs, order; mix_choice = mix_choice, kwargs...)
            cr = substitute_redundants(cr, de.scale_aons, de.names)
            eqs_out[i] = Symbolics.Equation(eqs[i].lhs, cr)
        end
    else
        for i = 1:length(eqs)
            cr = cumulant_expansion(eqs[i].rhs, order; mix_choice = mix_choice, kwargs...)
            cr = substitute_redundants(cr, de.scale_aons, de.names)
            eqs_out[i] = Symbolics.Equation(eqs[i].lhs, cr)
        end
    end

    return ScaledMeanfieldEquations(
        eqs_out,
        de.operator_equations,
        de.states,
        de.operators,
        de.hamiltonian,
        de.jumps,
        de.jumps_dagger,
        de.rates,
        de.iv,
        de.varmap,
        order,
        de.scale_aons,
        de.names,
        de.was_scaled,
    )
end

function _cumulant_expansion(args::Vector, order::Int)
    # Get all possible partitions; partitions(args,1) corresponds to the moment of order length(args)
    parts = [partitions(args, i) for i = 2:length(args)]

    args_sum = Any[]
    for p in parts
        for pj in p
            n = length(pj)
            args_prod = Any[-factorial(n-1)*(-1)^(n-1)]
            for p_ in pj # Product over partition blocks
                if length(p_) > order # If the encountered moment is larger than order, apply expansion
                    push!(args_prod, _cumulant_expansion(p_, order))
                else # Else, average and add its product
                    op_ = QMul(1, p_)
                    push!(args_prod, _average(op_))
                end
            end
            # Add terms in sum
            push!(args_sum, *(args_prod...))
        end
    end
    return average(+(args_sum...))
end

function _cumulant_expansion(x::IndexedAverageSum, order; kwargs...)
    return IndexedAverageSum(
        simplifyMultiplication(cumulant_expansion(x.term, order; kwargs...)),
        x.sum_index,
        x.non_equal_indices,
    )
end
function _cumulant_expansion(x::IndexedAverageDoubleSum, order; kwargs...)
    inner = _cumulant_expansion(x.innerSum, order; kwargs...)
    return IndexedAverageDoubleSum(inner, x.sum_index, x.non_equal_indices)
end
function _cumulant_expansion(a::BasicSymbolic{IndexedAverageSum}, order; kwargs...)
    if SymbolicUtils.hasmetadata(a, IndexedAverageSum)
        meta = TermInterface.metadata(a)[IndexedAverageSum]
        return _cumulant_expansion(meta, order; kwargs...)
    end
end
function _cumulant_expansion(a::BasicSymbolic{IndexedAverageDoubleSum}, order; kwargs...)
    if SymbolicUtils.hasmetadata(a, IndexedAverageDoubleSum)
        meta = TermInterface.metadata(a)[IndexedAverageDoubleSum]
        return _cumulant_expansion(meta, order; kwargs...)
    end
end
function _cumulant_expansion(a::BasicSymbolic{SpecialIndexedAverage}, order; kwargs...)
    if SymbolicUtils.hasmetadata(a, SpecialIndexedAverage)
        meta = TermInterface.metadata(a)[SpecialIndexedAverage]
        return SpecialIndexedAverage(
            cumulant_expansion(meta.term, order; kwargs...),
            meta.indexMapping,
        )
    end
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
function cumulant(op::QMul, n::Int = get_order(op); simplify = true, kwargs...)
    order = get_order(op)
    if order < n
        return zero(op)
    else
        c = _cumulant(op.args_nc, n)
        if simplify
            return SymbolicUtils.simplify(op.arg_c*c)
        else
            return op.arg_c*c
        end
    end
end
cumulant(avg::Average, args...; kwargs...) =
    cumulant(SymbolicUtils.arguments(avg)[1], args...; kwargs...)
cumulant(op::QSym, n::Int = 1; kwargs...) = isone(n) ? average(op) : zero(op)

function _cumulant(args::Vector, m::Int = length(args))
    parts = [partitions(args, i) for i = 1:m]
    args_sum = Any[]
    for i = 1:length(parts)
        p = collect(parts[i])
        for j = 1:length(p) # Terms in the sum
            n = length(p[j])
            args_prod = Any[factorial(n-1)*(-1)^(n-1)]
            for p_ in p[j] # Product over partition blocks
                push!(args_prod, average(QMul(1, p_))) #_average before
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
get_order(avg::Average) = get_order(SymbolicUtils.arguments(avg)[1])
function get_order(t::SymbolicUtils.Symbolic)
    if SymbolicUtils.iscall(t)
        return maximum(map(get_order, SymbolicUtils.arguments(t)))
    else
        return 0
    end
end
get_order(::Number) = 0
get_order(q::QMul) = length(q.args_nc)
function get_order(q::QAdd)
    order = Int[get_order(arg) for arg in SymbolicUtils.arguments(q)]
    return maximum(order)
end
get_order(::QSym) = 1

get_order(::IndexedOperator) = 1

get_order(x::SingleSum) = get_order(x.term)
get_order(x::SpecialIndexedTerm) = get_order(x.term)

get_order(a::IndexedAverageSum) = get_order(a.term)
get_order(a::IndexedAverageDoubleSum) = get_order(a.innerSum)
get_order(a::SpecialIndexedAverage) = get_order(a.term)
function get_order(a::BasicSymbolic{IndexedAverageSum})
    if SymbolicUtils.hasmetadata(a, IndexedAverageSum)
        meta = TermInterface.metadata(a)[IndexedAverageSum]
        return get_order(meta)
    end
end
function get_order(a::BasicSymbolic{IndexedAverageDoubleSum})
    if SymbolicUtils.hasmetadata(a, IndexedAverageDoubleSum)
        meta = TermInterface.metadata(a)[IndexedAverageDoubleSum]
        return get_order(meta)
    end
end
function get_order(a::BasicSymbolic{SpecialIndexedAverage})
    if SymbolicUtils.hasmetadata(a, SpecialIndexedAverage)
        meta = TermInterface.metadata(a)[SpecialIndexedAverage]
        return get_order(meta)
    end
end
get_order(x::NumberedOperator) = get_order(x.op)
