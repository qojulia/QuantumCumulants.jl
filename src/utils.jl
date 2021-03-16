"""
    find_missing(rhs::Vector, vs::Vector, vs_adj=_conj(vs), ps=[])

For a list of expressions contained in `rhs`, check whether all occurring symbols
are contained either in the variables given in `vs`. If a list of parameters `ps`
is provided, parameters that do not occur in the list `ps` are also added to the list.
Returns a list of missing symbols.
"""
function find_missing(rhs::Vector, vs::Vector; vs_adj::Vector=_conj.(vs), ps=[])
    missed = []
    for e=rhs
        append!(missed,get_symbolics(e))
    end
    unique!(missed)
    if isempty(ps)
        filter!(!SymbolicUtils.sym_isa(Parameter), missed)
    end
    filter!(x->!(_in(x, vs) || _in(x, ps) || _in(x, vs_adj)),missed)
    isempty(ps) || (ps_adj = _conj.(ps); filter!(x -> !_in(x,ps_adj), missed))
    return missed
end
function find_missing(de::HeisenbergEquation; kwargs...)
    find_missing(de.rhs, de.lhs; kwargs...)
end

"""
    _in(x, itr)

Same as `Base.in` but uses `isequal` instead of `==`.
"""
function _in(x, itr)
    anymissing = false
    for y in itr
        v = isequal(y, x)
        if ismissing(v)
            anymissing = true
        elseif v
            return true
        end
    end
    return anymissing ? missing : false
end

"""
    get_symbolics(ex)

Find all symbolic numbers occuring in `ex`.
"""
get_symbolics(x::Number) = []
function get_symbolics(t::SymbolicUtils.Symbolic)
    if SymbolicUtils.istree(t)
        if SymbolicUtils.is_operation(average)(t)
            return [t]
        else
            syms = []
            for arg in SymbolicUtils.arguments(t)
                append!(syms, get_symbolics(arg))
            end
            return unique(syms)
        end
    else
        return [t]
    end
end

"""
    complete(de::HeisenbergEquation)

From a set of differential equation of averages, find all averages that are missing
and derive the corresponding equations of motion.
"""
function complete(de::HeisenbergEquation;kwargs...)
    rhs_, lhs_ = complete(de.rhs,de.lhs,de.hamiltonian,de.jumps,de.rates;kwargs...)
    de_ = HeisenbergEquation(lhs_,rhs_,de.hamiltonian,de.jumps,de.rates,de.iv,copy(de.varmap))
    add_vars!(de_.varmap, lhs_, de.iv)
    return de_
end
function complete(rhs::Vector, vs::Vector, H, J, rates; order=nothing, filter_func=nothing, mix_choice=maximum, kwargs...)
    order_lhs = maximum(get_order.(vs))
    order_rhs = maximum(get_order.(rhs))
    if order isa Nothing
        order_ = max(order_lhs, order_rhs)
    else
        order_ = order
    end
    maximum(order_) >= order_lhs || error("Cannot form cumulant expansion of derivative; you may want to use a higher order!")

    vs_ = copy(vs)
    rhs_ = [cumulant_expansion(r, order_) for r in rhs]
    missed = unique_ops(find_missing(rhs_, vs_))
    filter!(SymbolicUtils.sym_isa(AvgSym),missed)
    isnothing(filter_func) || filter!(filter_func, missed) # User-defined filter
    while !isempty(missed)
        ops = [SymbolicUtils.arguments(m)[1] for m in missed]
        he = isempty(J) ? heisenberg(ops,H; kwargs...) : heisenberg(ops,H,J;rates=rates, kwargs...)
        he_avg = average(he,order_;mix_choice=mix_choice, kwargs...)
        rhs_ = [rhs_;he_avg.rhs]
        vs_ = [vs_;he_avg.lhs]
        missed = unique_ops(find_missing(rhs_,vs_))
        filter!(SymbolicUtils.sym_isa(AvgSym),missed)
        isnothing(filter_func) || filter!(filter_func, missed) # User-defined filter
    end

    if !isnothing(filter_func)
        # Find missing values that are filtered by the custom filter function,
        # but still occur on the RHS; set those to 0
        missed = unique_ops(find_missing(rhs_, vs_))
        filter!(SymbolicUtils.sym_isa(AvgSym),missed)
        filter!(!filter_func, missed)
        missed_adj = map(_adjoint, missed)
        subs = Dict(vcat(missed, missed_adj) .=> 0)
        rhs_ = [substitute(r, subs) for r in rhs_]
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

    all_ops = QSymbolic[]
    for i=1:order
        for c in combinations(ops, i)
            push!(all_ops, prod(c))
        end
    end

    # Simplify and remove non-operators iteratively
    ops_1 = map(qsimplify, all_ops)
    ops_2 = all_ops
    while !isequal(ops_1,ops_2)
        ops_2 = QSymbolic[]
        for op in ops_1
            append!(ops_2, _get_operators(op))
        end
        ops_1 = map(qsimplify, ops_2)
    end

    return unique_ops(ops_2)
end
find_operators(op::QSymbolic,args...) = find_operators(hilbert(op),args...)

"""
    hilbert(::QNumber)

Return the Hilbert space of the operator.
"""
hilbert(op::QSym) = op.hilbert
hilbert(t::QTerm) = hilbert(t.arguments[findfirst(x->isa(x,QSymbolic), t.arguments)])

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
    get_operators(::QNumber)

Return a list of all [`QSym`](@ref) in an expression.
"""
get_operators(x) = _get_operators(x)
function get_operators(t::QTerm)
    ops = QSymbolic[]
    for arg in SymbolicUtils.arguments(t)
        append!(ops, get_operators(arg))
    end
    return ops
end

_get_operators(::Number) = []
_get_operators(op::QSym) = [op]
function _get_operators(t::QTerm)
    f = SymbolicUtils.operation(t)
    if f===(*)
        args = QSymbolic[]
        for arg in SymbolicUtils.arguments(t)
            append!(args, _get_operators(arg))
        end
        isempty(args) && return args
        return [*(args...)]
    elseif f===(^)
        return [t]
    else
        ops = QSymbolic[]
        for arg in SymbolicUtils.arguments(t)
            append!(ops, _get_operators(arg))
        end
        return ops
    end
end

"""
    unique_ops(ops)

For a given list of operators, return only unique ones taking into account
their adjoints.
"""
function unique_ops(ops)
    seen = eltype(ops)[]
    ops_adj = _adjoint.(ops)
    for (op,op′) in zip(ops,ops_adj)
        if !(_in(op, seen) || _in(op′, seen))
            push!(seen, op)
        end
    end
    return seen
end

# Overload getindex to obtain solutions with averages
function Base.getindex(sol::SciMLBase.AbstractTimeseriesSolution, avg::SymbolicUtils.Term{<:AvgSym})
    tsym = sol.prob.f.indepsym # This is a bit hacky
    t = SymbolicUtils.Sym{Real}(tsym)
    syms = SciMLBase.getsyms(sol)
    var = _make_var(avg, t)
    sym = Symbolics.tosymbol(var)
    if sym∈syms
        return getindex(sol, var)
    else
        var_ = _make_var(_conj(avg), t)
        return map(conj, getindex(sol, var_))
    end
end
Base.getindex(sol::SciMLBase.AbstractTimeseriesSolution, op::QSymbolic) = getindex(sol, average(op))


# Internal functions
_conj(v::SymbolicUtils.Term{<:AvgSym}) = _average(adjoint(v.arguments[1]))
function _conj(v::SymbolicUtils.Symbolic)
    if SymbolicUtils.istree(v)
        f = SymbolicUtils.operation(v)
        args = map(_conj, SymbolicUtils.arguments(v))
        return SymbolicUtils.similarterm(v, f, args)
    else
        return conj(v)
    end
end
_conj(x::Number) = conj(x)

_adjoint(op::QSymbolic) = adjoint(op)
_adjoint(s::SymbolicUtils.Symbolic{<:Number}) = _conj(s)
_adjoint(x) = adjoint(x)

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
_to_expression(op::QSym) = op.name
_to_expression(op::Create) = :(dagger($(op.name)))
_to_expression(op::Transition) = :(Transition($(op.name),$(op.i),$(op.j)) )
_to_expression(t::QTerm) = :( $(Symbol(t.f))($(_to_expression.(t.arguments)...)) )
_to_expression(p::Parameter) = p.name
function _to_expression(s::SymbolicUtils.Symbolic)
    if SymbolicUtils.istree(s)
        f = SymbolicUtils.operation(s)
        args = map(_to_expression, SymbolicUtils.arguments(s))
        return :( $(Symbol(f))($(args...)) )
    else
        return nameof(s)
    end
end
