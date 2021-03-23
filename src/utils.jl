"""
    find_missing(rhs::Vector, vs::Vector, vs_adj=_conj(vs), ps=[])

For a list of expressions contained in `rhs`, check whether all occurring symbols
are contained either in the variables given in `vs`. If a list of parameters `ps`
is provided, parameters that do not occur in the list `ps` are also added to the list.
Returns a list of missing symbols.
"""
function find_missing(rhs::Vector, vs::Vector; vs_adj::Vector=_conj.(vs), get_adjoints=true)
    missed = []
    missed_hashes = UInt[]
    vhash = map(hash, vs)
    v′hash = map(hash, vs_adj)
    for r∈rhs
        find_missing!(missed, missed_hashes, r, vhash, v′hash; get_adjoints=get_adjoints)
    end
    return missed
end
function find_missing!(missed, missed_hashes, r::SymbolicUtils.Symbolic, vhash, vs′hash; get_adjoints=true)
    if SymbolicUtils.istree(r)
        for arg∈SymbolicUtils.arguments(r)
            find_missing!(missed, missed_hashes, arg, vhash, vs′hash; get_adjoints=get_adjoints)
        end
    end
    return missed
end
function find_missing!(missed, missed_hashes, r::Average, vhash, vs′hash; get_adjoints=true)
    rhash = hash(r)
    if !(rhash ∈ vhash) && !(rhash ∈ vs′hash) && !(rhash ∈ missed_hashes)
        push!(missed, r)
        push!(missed_hashes, rhash)
        if !get_adjoints
            # To avoid collecting adjoints as missing variables,
            # collect the hash of the adjoint right away
            r′ = _conj(r)
            r′hash = hash(r′)
            if !(r′hash ∈ missed_hashes)
                push!(missed_hashes, r′hash)
            end
        end
    end
    return missed
end
find_missing!(missed, missed_hashes, r::Number, vs, vs′hash; kwargs...) = missed

function find_missing(eqs::Vector, vhash::Vector{UInt}, vs′hash::Vector{UInt}; get_adjoints=true)
    missed = []
    missed_hashes = UInt[]
    for i=1:length(eqs)
        find_missing!(missed, missed_hashes, eqs[i].rhs, vhash, vs′hash; get_adjoints=get_adjoints)
    end
    return missed
end

function find_missing(he::HeisenbergEquation; vs_adj=nothing, get_adjoints=true)
    vs = he.states
    vhash = map(hash, vs)
    vs′ = if vs_adj===nothing
        map(_conj, vs)
    else
        vs_adj
    end
    vs′hash = map(hash, vs′)
    filter!(!in(vhash), vs′hash)

    missed = []
    missed_hashes = UInt[]

    eqs = he.equations
    for i=1:length(eqs)
        find_missing!(missed, missed_hashes, eqs[i].rhs, vhash, vs′hash; get_adjoints=get_adjoints)
    end
    return missed
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
    complete!(de::HeisenbergEquation)

From a set of differential equation of averages, find all averages that are missing
and derive the corresponding equations of motion.
"""
function complete(de::HeisenbergEquation;kwargs...)
    de_ = deepcopy(de)
    complete!(de_;kwargs...)
    return de_
end
function complete!(de::HeisenbergEquation;
                                order=nothing,
                                filter_func=nothing,
                                mix_choice=maximum,
                                simplify=true,
                                kwargs...)
    vs = de.states
    order_lhs = maximum(get_order.(vs))
    order_rhs = 0
    for i=1:length(de.equations)
        k = get_order(de.equations[i].rhs)
        k > order_rhs && (order_rhs = k)
    end
    if order isa Nothing
        order_ = max(order_lhs, order_rhs)
    else
        order_ = order
    end
    maximum(order_) >= order_lhs || error("Cannot form cumulant expansion of derivative; you may want to use a higher order!")

    if order_ != de.order
        for i=1:length(de.equations)
            lhs = de.equations[i].lhs
            rhs = cumulant_expansion(de.equations[i].rhs, order_; simplify=simplify)
            de.equations[i] = Symbolics.Equation(lhs, rhs)
        end
    end

    vhash = map(hash, vs)
    vs′ = map(_conj, vs)
    vs′hash = map(hash, vs′)
    filter!(!in(vhash), vs′hash)
    missed = find_missing(de.equations, vhash, vs′hash; get_adjoints=false)
    isnothing(filter_func) || filter!(filter_func, missed) # User-defined filter

    while !isempty(missed)
        ops_ = [SymbolicUtils.arguments(m)[1] for m in missed]
        he = heisenberg(ops_,de.hamiltonian,de.jumps;
                                rates=de.rates,
                                simplify=simplify,
                                expand=true,
                                order=order_,
                                kwargs...)

        _append!(de, he)

        vhash_ = hash.(he.states)
        vs′hash_ = hash.(_conj.(he.states))
        append!(vhash, vhash_)
        for i=1:length(vhash_)
            vs′hash_[i] ∈ vhash_ || push!(vs′hash, vs′hash_[i])
        end

        missed = find_missing(he.equations, vhash, vs′hash; get_adjoints=false)
        isnothing(filter_func) || filter!(filter_func, missed) # User-defined filter
    end

    if !isnothing(filter_func)
        # Find missing values that are filtered by the custom filter function,
        # but still occur on the RHS; set those to 0
        missed = find_missing(de.equations, vhash, vs′hash; get_adjoints=false)
        filter!(!filter_func, missed)
        missed_adj = map(_adjoint, missed)
        subs = Dict(vcat(missed, missed_adj) .=> 0)
        for i=1:length(de.equations)
            de.equations[i] = substitute(de.equations[i], subs)
            de.states[i] = de.equations[i].lhs
        end
        if false
            for i=1:length(de.equations)
                de.equations[i] = SymbolicUtils.simplify(de.equations[i])
            end
        end
    end
    return de
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

    all_ops = QNumber[]
    for i=1:order
        for c in combinations(ops, i)
            push!(all_ops, prod(c))
        end
    end

    # Simplify and remove non-operators iteratively
    ops_1 = map(qsimplify, all_ops)
    ops_2 = all_ops
    while !isequal(ops_1,ops_2)
        ops_2 = QNumber[]
        for op in ops_1
            append!(ops_2, _get_operators(op))
        end
        ops_1 = map(qsimplify, ops_2)
    end

    return unique_ops(ops_2)
end
find_operators(op::QNumber,args...) = find_operators(hilbert(op),args...)

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
    ops = QNumber[]
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
        args = QNumber[]
        for arg in SymbolicUtils.arguments(t)
            append!(args, _get_operators(arg))
        end
        isempty(args) && return args
        return [*(args...)]
    elseif f===(^)
        return [t]
    else
        ops = QNumber[]
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
    var = make_var(avg, t)
    sym = Symbolics.tosymbol(var)
    if sym∈syms
        return getindex(sol, var)
    else
        var_ = make_var(_conj(avg), t)
        return map(conj, getindex(sol, var_))
    end
end
Base.getindex(sol::SciMLBase.AbstractTimeseriesSolution, op::QNumber) = getindex(sol, average(op))


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

_adjoint(op::QNumber) = adjoint(op)
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
