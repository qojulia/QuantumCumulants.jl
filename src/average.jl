import SymPy
using Combinatorics: partitions, combinations

"""
    get_order(::AbstractOperator)

Returns the order of the given operator or expression.

# Arguments
*`op::AbstractOperator`: The operator or expression of which the order is to
    be computed.
"""
get_order(::BasicOperator) = 1
function get_order(a::Prod)
    check = Bool[isa(a1,AbstractOperator) for a1=a.args]
    check0 = Bool[isa(a1,Identity) || isa(a1, Zero) for a1=a.args[check]]
    return sum(check) - sum(check0)
end
function get_order(a::TensorProd)
    args_ = []
    for a1=a.args
        if isa(a1,AbstractOperator)
            isa(a1,Identity) || isa(a1,Zero) || push!(args_,a1)
        end
    end
    isempty(args_) && return 0
    return sum(get_order.(args_))::Int
end
get_order(a::Add) = maximum(get_order.(a.args))::Int
get_order(s::SumType) = get_order(s.args[1])::Int

"""
    average(op::AbstractOperator,order::Int=get_order(op))

Compute the average value (moment) of the operator or expression `op`.
Moments that are of an order higher than the given `order`, i.e. average values
of operators that occur in the expression with larger order, are expanded
using the `cumulant_expansion` function. The resulting operator/expression
will hence be of order ≦ `order`.

# Arguments
*`op::AbstractOperator`: The operator or expression of which the average is to be
    computed.
*`order::Int=get_order(op)`: The order to which encountered moments should be reduced.
"""
function average(a::BasicOperator,order::Int=1)
    return SymPy.symbols("⟨"*gen_label(a)*"⟩")
end
function average(a::Transition,order::Int=1)
    label = gen_label(a)
    if a.i==a.j==a.GS
        return average(simplify_operators(a),order)
    else
        return SymPy.symbols("⟨"*label*"⟩",real=ishermitian(a))
    end
end
average(::Identity,order::Int=1) = 1
average(::Zero,order::Int=1) = 0
function average(a::IndexedOperator,order::Int=1)
    avg = average(a.operator,order)
    i = sympify(a.index)
    b = SymPy.sympy.IndexedBase(avg)
    return b[i]
end

# Expression averaging
function average(a::Prod,order::Int=get_order(a))
    if get_order(a) <= order
        return _average_proper_order(a)
    else
        if isa(a.args[1],Number)
            return a.args[1]*cumulant_expansion(a.f,a.args[2:end],order)
        else
            return cumulant_expansion(a.f,a.args,order)
        end
    end
end
function average(a::TensorProd,order::Int=get_order(a))
    if get_order(a) <= order
        return _average_proper_order(a)
    else
        c, args_ = _flatten_prods(a.args)
        return c*cumulant_expansion(a.f,args_,order)
    end
end
average(a::Add,order::Int=get_order(a)) = sum(average(a1,order) for a1=a.args)
function average(s::SumType,order::Int=get_order(s))
    arg = average(s.args[1],order)
    i = sympify(s.f.index)
    l = i.__pyobject__.lower
    u = i.__pyobject__.upper
    return SymPy.sympy.Sum(arg,(i,l,u))
end

"""
    average(op::AbstractOperator,order::Vector{<:Int};mix_choice::Function=maximum)

Average an operator or expression, where `order` is a vector specifying the
order to which each subspace on which the operator acts nontrivially is to be
expanded via cumulant expansion. For example, the average of a tensor product of
the form `A⊗𝟙` with `order=[2,1]` will be expanded to second order, whereas
`𝟙⊗A` would be expanded to first order.
If an operator acts nontrivially on multiple subspaces, the order is chosen according
to the optional argument `mix_choice`, which defaults to `maximum`. If you want
to reduce the order here as well, set `mix_choice=minimum`.

# Arguments
*`op::AbstractOperator`: The operator or expression of which the average is to be
    computed.
*`order::Vector{<:Int}``: A vector specifying to which order each subspace should be
    reduced.

# Optional arguments:
*`mix_choice::Function=maximum`: Choose by which function the order to which an operator
    that acts on different subspaces should be reduced.
"""
function average(a::AbstractOperator,order::Vector{<:Int};mix_choice::Function=maximum)
    aon = acts_on(a)
    isempty(aon) && return average(a)
    order_ = mix_choice(order[aon])
    return average(a,order_)
end
average(a::Add,order::Vector{Int};kwargs...) = sum(average(a1,order;kwargs...) for a1=a.args)

"""
    average(de::DifferentialEquation,order::Int=get_order(de.lhs))

Compute the average up to `order` for a given `DifferentialEquation`.
See also `average(::AbstractOperator,::Int)`. Returns a `DifferentialEquation`
of averages.

# Arguments
*`de::DifferentialEquation`: The equation of which the average is to be computed.
*`order::Int=get_order(de.lhs)`: The order to which encountered moments should be
    reduced. Must be larger or equal than the order of the operator on the
    left-hand side of the equation.
"""
function average(de::DifferentialEquation,order::Int=get_order(de.lhs))
    get_order(de.lhs) <= order || error("The left-hand-side operator is larger than the given order!")
    return DifferentialEquation(average(de.lhs,order),average(de.rhs,order))
end

"""
    average(de::DifferentialEquation,order::Vector{<:Int};kwargs...)

Compute the mixed-order average of a `DifferentialEquation`. See also
`average(de::DifferentialEquation,order::Int)` and
`average(op::AbstractOperator,order::Vector{<:Int}`.
"""
function average(de::DifferentialEquation,order::Vector{<:Int};kwargs...)
    return DifferentialEquation(average(de.lhs,order;kwargs...),average(de.rhs,order;kwargs...))
end

"""
    average(de::DifferentialEquationSet,order::Int=maximum(get_order.(de.lhs)))

Compute average of a set of equations provided as `DifferentialEquationSet`.
The advantage of computing the average directly from a set is that a closed
system can be more easily found since averages of adjoint operators can be
replaced by the respective complex conjugates.
Returns a `DifferentialEquationSet` of averages.

# Arguments
*`de::DifferentialEquationSet`: The set of equations of which the averages are to
    be computed.
*`order::Int=maximum(get_order.(lhs))`: The order to which encountered moments
    should be reduced. Must be larger or equal than the order of any operator on
    the left-hand side of the equations.
"""
function average(de::DifferentialEquationSet,order::Int=maximum(get_order.(de.lhs)))
    maximum(get_order.(de.lhs)) <= order || error("The left-hand-side operator is larger than the given order!")

    # Average over equations
    lhs_avg = average.(de.lhs,order)
    rhs_avg = average.(de.rhs,order)
    rhs_ = replace_adjoints(rhs_avg,de.lhs,order)
    return DifferentialEquationSet(lhs_avg,rhs_)
end

"""
    average(de::DifferentialEquationSet,order::Int=maximum(get_order.(de.lhs)))

Compute mixed-order average of a set of equations provided as
`DifferentialEquationSet`. See also
`average(de::DifferentialEquationSet,order::Int)` and
`average(op::AbstractOperator,order::Vector{<:Int}`.
"""
function average(de::DifferentialEquationSet,order::Vector{<:Int};kwargs...)
    # Average over equations
    lhs_avg = [average(de.lhs[i],order;kwargs...) for i=1:length(de.lhs)]
    rhs_avg = [average(de.rhs[i],order;kwargs...) for i=1:length(de.lhs)]
    rhs_ = replace_adjoints(rhs_avg,de.lhs,order)
    return DifferentialEquationSet(lhs_avg,rhs_)
end

"""
    replace_adjoints(exprs::Vector{<:SymPy.Sym},ops::Vector{<:AbstractOperator},order::Int)

For a set of expressions containing averages of `ops` up to `order`, replace
possibly occurring averages that correspond to the hermitian conjugate of an
operator contained in `ops` by the complex conjugate of the average of the
operator. Returns a vector of the new expressions. Note: this function is used
in order to obtain a closed set of equations when averaging over a
`DifferentialEquationSet`.

# Arguments
*`exprs::Vector{<:SymPy.Sym}`: A vector of averages given as `SymPy.Sym`.
*`ops::Vector{<:AbstractOperator}`: Vector of operator
*`order::Int`
"""
function replace_adjoints(exprs::Vector{<:SymPy.Sym},ops::Vector{<:AbstractOperator},order::Int=maximum(get_order.(ops)))
    check_adj = [!ishermitian(op) for op=ops]
    if any(check_adj)
        ops_adj = ops[check_adj]
        op_avg = average.(ops_adj,order)
        adj = adjoint.(ops_adj)
        adj_avg = average.(adj,order)
        subs_adj = Dict{SymPy.Sym,SymPy.Sym}()
        for i=1:length(adj_avg)
            subs_adj[adj_avg[i]] = op_avg[i]'
        end
        exprs_ = [ex(subs_adj) for ex=exprs]
        # TODO: find cleaner solution to substitute in sums
        return replace_adjoints_indexed(exprs_,subs_adj)
    else
        return exprs
    end
end
function replace_adjoints(exprs::Vector{<:SymPy.Sym},ops::Vector{<:AbstractOperator},order::Vector{<:Int};kwargs...)
    check_adj = [!ishermitian(op) for op=ops]
    if any(check_adj)
        ops_adj = ops[check_adj]
        op_avg = [average(ops_adj[i],order;kwargs...) for i=1:length(ops_adj)]
        adj = adjoint.(ops_adj)
        adj_avg = [average(adj[i],order;kwargs...) for i=1:length(ops_adj)]
        subs_adj = Dict{SymPy.Sym,SymPy.Sym}()
        for i=1:length(adj_avg)
            subs_adj[adj_avg[i]] = op_avg[i]'
        end
        exprs_ = [ex(subs_adj) for ex=exprs]
        # TODO: find cleaner solution to substitute in sums
        return replace_adjoints_indexed(exprs_,subs_adj)
    else
        return exprs
    end
end
replace_adjoints(expr::SymPy.Sym,ops::Vector,order;kwargs...) = replace_adjoints([expr],ops,order;kwargs...)
replace_adjoints(expr::SymPy.Sym,ops::AbstractOperator,order;kwargs...) = replace_adjoints([expr],[ops],order;kwargs...)

function replace_adjoints_indexed(exprs, subs)
    exprs_ = SymPy.Sym[]
    for i=1:length(exprs)
        push!(exprs_, exprs[i])
        for (keys,vals)=(subs)
            # Check for indexed objects; others have already been replaced
            if classname(keys)==classname(vals')=="Indexed"
                k_inds = keys.__pyobject__.indices
                m_inds = vals'.__pyobject__.indices
                length(k_inds)==length(m_inds)==1 || continue# TODO: other cases
                k_ = IndexOrder[findfirst(x->sympify(x)==k_inds[1],IndexOrder)]
                m_ = IndexOrder[findfirst(x->sympify(x)==m_inds[1],IndexOrder)]
                # Replace indices by any other known index and try to substitute adjoint in expression
                for j=IndexOrder
                    key_ = swap_index(keys, k_, j)
                    val_ = swap_index(vals, m_, j)
                    exprs_[i] = exprs_[i].__pyobject__.replace(key_,val_)
                end
            end
        end
    end
    return exprs_
end

function _flatten_prods(args)
    c = 1
    args_ = []
    for a1=args
        if isa(a1,Prod)
            if isa(a1.args[1],Number)
                c *= a1.args[1]
                for a2=a1.args[2:end]
                    iszero(a2) && (c *= 0)
                    isa(a2,Identity) || push!(args_,a2)
                end
            else
                for a2=a1.args
                    iszero(a2) && (c *= 0)
                    isa(a2,Identity) || push!(args_,a2)
                end
            end
        elseif isa(a1,BasicOperator)
            iszero(a1) && (c *= 0)
            isa(a1,Identity) || push!(args_,a1)
        elseif isa(a1,Number)
            c *= a1
        else
            error()
        end
    end
    return c, args_
end

function _average_proper_order(a::Prod)
    if isa(a.args[1],Number)
        return a.args[1]*_average_proper_order(prod(a.args[2:end]))
    else
        label = "⟨"*prod(gen_label.(a.args))*"⟩"
        s = SymPy.symbols(label)
        if any([isa(a1,IndexedOperator) for a1=a.args])
            _, inds = find_index(a)
            av = SymPy.sympy.IndexedBase(s)
            return av[sympify.(inds)...]
        else
            return s
        end
    end
end
_average_proper_order(a::BasicOperator) = average(a)
function _average_proper_order(a::TensorProd)
    c, args_ = _flatten_prods(a.args)
    if isempty(args_)
        out = 1
    elseif length(args_)==1
        out = _average_proper_order(args_[1])
    else
        label = "⟨"*prod(gen_label.(args_))*"⟩"
        out_ = SymPy.symbols(label)
        _, inds = find_index(a)
        if isempty(inds)
            out = out_
        else
            av = SymPy.sympy.IndexedBase(out_)
            out = av[sympify.(inds)...]
        end
    end
    return c*out
end
_average_proper_order(a::Add) = error("Something went wrong here!")

gen_label(a::BasicOperator) = string(a.label)
gen_label(a::Create) = string(a.label)*"ᵗ"
gen_label(a::Transition) = string(a.label)*"^{"*string(a.i)*string(a.j)*"}"
gen_label(a::IndexedOperator) = gen_label(a.operator)


"""
    cumulant_expansion(f, args, order::Int)

For a product of operators whose constituents are given in `args`, expand it in terms
of moments up to `order` neglecting their joint cumulant.

See also: https://en.wikipedia.org/wiki/Cumulant#Joint_cumulants
"""
function cumulant_expansion(f::Function,args::Vector,order::Int)
    length(args) > order  || error("Something went wrong here!")

    # Get all possible partitions; partitions(args,1) corresponds to the moment of order length(args)
    parts = [collect(partitions(args,i)) for i=2:length(args)]

    c = 0.0
    for i=1:length(parts)
        p = parts[i]
        for j=1:length(p) # Terms in the sum
            n = length(p[j])
            c_ = -factorial(n-1)*(-1)^(n-1)
            for p_=p[j] # Product over partition blocks
                if length(p_) > order # If the encountered moment is larger than order, apply expansion
                    c_ *= cumulant_expansion(f,p_, order)
                else # Else, average and add its product
                    # op_ = apply_comms(f(p_...)) # Can be necessary for order > 2
                    # c_ *= average(op_,order)
                    c_ *= _average_proper_order(f(p_...))
                end
            end
            # Add terms in sum
            c += c_
        end
    end
    return SymPy.expand(c)
end

function combine_operators(ops::Vector{<:AbstractOperator}, order::Int)
    @assert sum(get_order.(ops)) == length(ops)
    if order==1
        return ops
    else
        ops_tot = [adjoint.(ops); ops]
        ops_combs = combinations(ops_tot,order)
        ops_out = AbstractOperator[]
        for op=ops_combs
            op_ = simplify_operators(prod(op))
            (op_ ∈ ops_out) || push!(ops_out, op_)
        end
        return ops_out
    end
end
function combine_operators(ops::Vector{<:AbstractOperator}, order::Vector{<:Int})
    ops_ = AbstractOperator[]
    for ord=order
        append!(ops_, combine_operators(ops, ord))
    end
    return ops_
end
