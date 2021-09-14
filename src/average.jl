function average end

"""
    AvgSym <: CNumber

Symbolic number representing the average over an operator.
See also: [`average`](@ref)
"""
struct AvgSym <: CNumber end

const Average = SymbolicUtils.Term{<:AvgSym}

const sym_average = begin # Symbolic function for averages
    T = SymbolicUtils.FnType{Tuple{QNumber}, AvgSym}
    SymbolicUtils.Sym{T}(:avg)
end

# Type promotion -- average(::QNumber)::Number
SymbolicUtils.promote_symtype(::typeof(sym_average), ::Type{<:QNumber}) = AvgSym

# Direct construction of average symbolic expression
function _average(operator)
    return SymbolicUtils.Term{AvgSym}(sym_average, [operator])
end

function acts_on(s::SymbolicUtils.Symbolic)
    if SymbolicUtils.istree(s)
        f = SymbolicUtils.operation(s)
        if f === sym_average
            return acts_on(SymbolicUtils.arguments(s)[1])
        else
            aon = []
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
    average(::QNumber,order)

Compute the average of an operator. If `order` is given, the [`cumulant_expansion`](@ref)
up to that order is computed immediately.
"""
average(op::QSym) = _average(op)
function average(op::QTerm)
    f = SymbolicUtils.operation(op)
    if f===(+) || f===(-) # linearity
        args = map(average, SymbolicUtils.arguments(op))
        return f(args...)
    elseif f === (*)
        # Move constants out of average
        c = op.arg_c
        op_ = QMul(1,op.args_nc)
        return c*_average(op_)
    else
        error("Unknown function $f")
    end
end
average(x::SNuN) = x
average(x,order;kwargs...) = cumulant_expansion(average(x),order;kwargs...)

function undo_average(t)
    if SymbolicUtils.istree(t)
        f = SymbolicUtils.operation(t)
        if f === sym_average
            return SymbolicUtils.arguments(t)[1]
        else
            args = map(undo_average, SymbolicUtils.arguments(t))
            return f(args...)
        end
    else
        return t
    end
end

function undo_average(eq::Symbolics.Equation)
    lhs = undo_average(eq.lhs)
    rhs = undo_average(eq.rhs)
    return Symbolics.Equation(lhs,rhs)
end


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
function cumulant_expansion(x::SymbolicUtils.Symbolic,order::Integer;simplify=true,kwargs...)
    if SymbolicUtils.istree(x)
        get_order(x) <= order && return x
        f = SymbolicUtils.operation(x)
        if f===sym_average
            op = SymbolicUtils.arguments(x)[1]
            return _cumulant_expansion(op.args_nc, order)
        else
            args = SymbolicUtils.arguments(x)
            cumulants = [cumulant_expansion(arg,order;kwargs...) for arg in args]
            return f(cumulants...)
        end
    else
        return x
    end
end
function cumulant_expansion(avg::SymbolicUtils.Term{<:AvgSym},order::Vector;mix_choice=maximum,kwargs...)
    aon = acts_on(avg)
    order_ = mix_choice(order[get_i(i)] for i in aon)
    return cumulant_expansion(avg,order_;kwargs...)
end
cumulant_expansion(x::Number,order;kwargs...) = x

function cumulant_expansion(x::SymbolicUtils.Symbolic,order;mix_choice=maximum,simplify=true,kwargs...)
    if SymbolicUtils.istree(x)
        f = SymbolicUtils.operation(x)
        args = SymbolicUtils.arguments(x)
        cumulants = [cumulant_expansion(arg,order;simplify=simplify,mix_choice=mix_choice) for arg in args]
        if simplify
            return SymbolicUtils.simplify(f(cumulants...);kwargs...)
        else
            return f(cumulants...)
        end
    else
        return x
    end
end
function cumulant_expansion(de::MeanfieldEquations,order;multithread=false,mix_choice=maximum,kwargs...)
    order==de.order && return de
    eqs = de.equations
    eqs_out = Vector{Symbolics.Equation}(undef, length(eqs))
    if multithread
        Threads.@threads for i=1:length(eqs)
            cr = cumulant_expansion(eqs[i].rhs,order;mix_choice=mix_choice,kwargs...)
            eqs_out[i] = Symbolics.Equation(eqs[i].lhs, cr)
        end
    else
        for i=1:length(eqs)
            cr = cumulant_expansion(eqs[i].rhs,order;mix_choice=mix_choice,kwargs...)
            eqs_out[i] = Symbolics.Equation(eqs[i].lhs, cr)
        end
    end
    return MeanfieldEquations(eqs_out,de.operator_equations,de.states,de.operators,
                            de.hamiltonian,de.jumps,de.jumps_dagger,de.rates,de.iv,de.varmap,
                            order)
end
function cumulant_expansion(de::ScaledMeanfieldEquations,order;multithread=false,mix_choice=maximum,kwargs...)
    order==de.order && return de
    eqs = de.equations
    eqs_out = Vector{Symbolics.Equation}(undef, length(eqs))
    if multithread
        Threads.@threads for i=1:length(eqs)
            cr = cumulant_expansion(eqs[i].rhs,order;mix_choice=mix_choice,kwargs...)
            cr = substitute_redundants(cr, de.scale_aons, de.names)
            eqs_out[i] = Symbolics.Equation(eqs[i].lhs, cr)
        end
    else
        for i=1:length(eqs)
            cr = cumulant_expansion(eqs[i].rhs,order;mix_choice=mix_choice,kwargs...)
            cr = substitute_redundants(cr, de.scale_aons, de.names)
            eqs_out[i] = Symbolics.Equation(eqs[i].lhs, cr)
        end
    end

    return ScaledMeanfieldEquations(eqs_out,de.operator_equations,de.states,de.operators,
                                    de.hamiltonian,de.jumps,de.jumps_dagger,de.rates,de.iv,
                                    de.varmap,order,
                                    de.scale_aons,de.names,de.was_scaled
                                    )
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
function cumulant(op::QMul,n::Int=get_order(op);simplify=true,kwargs...)
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
                push!(args_prod, _average(QMul(1, p_)))
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
get_order(q::QMul) = length(q.args_nc)
function get_order(q::QAdd)
    order = Int[get_order(arg) for arg in SymbolicUtils.arguments(t)]
    return maximum(order)
end
get_order(::QSym) = 1
