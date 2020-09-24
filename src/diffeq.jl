import MacroTools

"""
    build_ode(rhs::Vector, vs::Vector, ps=[], usym=:u,
                psym=:p, tsym=:t; set_unknowns_zero::Bool=false, check_bounds::Bool=false)

From a set of equations contained in `eqs`, generate a `Meta.Expr` containing the
code for a function which can be directly passed to `OrdinaryDiffEq` in order to solve
it. The variable vector `u` corresponds to the symbols provided in `vs`.

# Arguments
*`eqs::Vector`: The vector containing the right-hand side of equations.
*`vs::Vector`: The vector containing the left-hand side of equations.
*`ps=[]`: List of (symbolic) parameters, which are parsed into parameters
    used in DiffEq functions.
*`usym=:u`: The symbol used for the variable vector.
*`psym=:p`: The symbol used for the parameter vector.
*`tsym=:t`: The symbol used for the time parameter.

# Optional arguments
*`check_bounds::Bool=false`: Choose whether the resulting function should contain
    the `@inbounds` flag, which skips bounds checking for performance.
"""
function build_ode(rhs::Vector, vs::Vector, ps=[], usym=:u, psym=:p, tsym=:t;
                    check_bounds::Bool=false, idx_borders=nothing)
    @assert length(rhs) == length(vs)
    if any(has_indexed.(vs)) || any(has_indexed.(rhs))
        return _build_indexed_ode(rhs, vs, ps, usym, psym, tsym, check_bounds, idx_borders)
    else
        return _build_ode(rhs, vs, ps, usym, psym, tsym, check_bounds)
    end
end

function _build_ode(rhs, vs, ps, usym, psym, tsym, check_bounds)
    vs_adj_ = adjoint.(vs)

    # Check if there are unknown symbols
    missed = find_missing(rhs,vs;vs_adj=vs_adj_,ps=ps)
    isempty(missed) || throw_missing_error(missed)

    dusym = Symbol(:d,usym)
    us = [:($usym[$i]) for i=1:length(vs)]
    dus = [:($dusym[$i]) for i=1:length(vs)]

    vs_ = _to_expression.(vs)
    vs_adj = _to_expression.(vs_adj_)
    rhs_ = _to_expression.(rhs)
    function _pw_func(x)
        if x in vs_
            i = findfirst(isequal(x),vs_)
            return us[i]
        elseif x in vs_adj
            i = findfirst(isequal(x),vs_adj)
            return :( conj($(us[i])) )
        else
            return x
        end
    end
    rhs_ = [MacroTools.postwalk(_pw_func, r) for r in rhs_]

    if !isempty(ps)
        ps_ = _to_expression.(ps)
        psyms = [:($psym[$i]) for i=1:length(ps)]
        rhs_ = [MacroTools.postwalk(x -> (x in ps_) ? psyms[findfirst(isequal(x), ps_)] : x, r) for r in rhs_]
    end

    # From https://github.com/JuliaDiffEq/ModelingToolkit.jl/blob/dca5f38491ae6dea431cb2a7cceb055645086034/src/utils.jl#L44
    line_eqs = [Expr(:(=), dus[i], rhs_[i]) for i=1:length(us)]
    var_eqs = build_expr(:block, line_eqs)

    fargs = :($dusym,$usym,$psym,$tsym)
    if check_bounds
        f_ex = :(
            ($fargs) -> begin
                begin
                    $var_eqs
                end
                return nothing
            end
        )
    else
        f_ex = :(
            ($fargs) -> begin
                @inbounds begin
                    $var_eqs
                end
                return nothing
            end
        )
    end
    return f_ex
end

function _build_indexed_ode(rhs, vs, ps, usym, psym, tsym, check_bounds, idx_borders)
    idx_borders === nothing && error("Need number of elements for indexes as numbers!")
    # Check if there are unknown symbols
    missed = find_missing(rhs,vs;vs_adj=adjoint.(vs))#,ps=ps) TODO check parameters without indices
    isempty(missed) || throw_missing_error(missed)

    _vs, _rhs = expand_indexed(vs, rhs, idx_borders)
    vs_adj_ = adjoint.(_vs)

    dusym = Symbol(:d,usym)
    us = [:($usym[$i]) for i=1:length(_vs)]
    dus = [:($dusym[$i]) for i=1:length(_vs)]

    vs_ = _to_expression.(_vs)
    vs_adj = _to_expression.(vs_adj_)
    rhs_ = _to_expression.(_rhs)
    function _pw_func(x)
        if x in vs_
            i = findfirst(isequal(x),vs_)
            return us[i]
        elseif x in vs_adj
            i = findfirst(isequal(x),vs_adj)
            return :( conj($(us[i])) )
        elseif MacroTools.@capture(x, IndexedParameter(p_, idx_))
            return :(($p)[$(idx...)])
        else
            return x
        end
    end
    rhs_ = [MacroTools.postwalk(_pw_func, r) for r in rhs_]

    if !isempty(ps)
        ps_ = _to_expression.(ps)
        psyms = [:($psym[$i]) for i=1:length(ps)]
        rhs_ = [MacroTools.postwalk(x -> (x in ps_) ? psyms[findfirst(isequal(x), ps_)] : x, r) for r in rhs_]
    end

    # From https://github.com/JuliaDiffEq/ModelingToolkit.jl/blob/dca5f38491ae6dea431cb2a7cceb055645086034/src/utils.jl#L44
    line_eqs = [Expr(:(=), dus[i], rhs_[i]) for i=1:length(us)]
    var_eqs = build_expr(:block, line_eqs)

    fargs = :($dusym,$usym,$psym,$tsym)
    if check_bounds
        f_ex = :(
            ($fargs) -> begin
                begin
                    $var_eqs
                end
                return nothing
            end
        )
    else
        f_ex = :(
            ($fargs) -> begin
                @inbounds begin
                    $var_eqs
                end
                return nothing
            end
        )
    end
    return f_ex
end

expand_indexed(de::DifferentialEquation, idx_borders) = DifferentialEquation(expand_indexed(de.lhs, de.rhs, idx_borders)..., de.hamiltonian, de.jumps, de.rates)
expand_indexed(vs, rhs, idx_borders::Pair) = expand_indexed(vs, rhs, [idx_borders])
function expand_indexed(vs, rhs, idx_borders::Vector)
    rhs_exp = expand_sums.(rhs, idx_borders)
    vs_ = Number[]
    rhs_ = Number[]
    for i=1:length(vs)
        idx = find_index(vs[i])
        if isempty(idx)
            push!(vs_, vs[i])
            push!(rhs_, rhs_exp[i])
        else
            counts = getfield.(idx, :count)
            if idx_borders isa Vector
                borders = [idx_borders[findfirst(x->isequal(x[1],c), idx_borders)][2] for c in counts]
            else
                @assert idx_borders isa Pair
                @assert all((isequal(c,idx_borders[1]) for c in counts))
                borders = [idx_borders[1] for k=1:length(idx)]
            end
            idx_combs = combinations(idx,2)
            skip_equal_inds = []
            for c in idx_combs
                if has_expr(c[1] != c[2], vs[i])
                    push!(skip_equal_inds, c)
                elseif has_expr(c[2] != c[1], vs[i])
                    push!(skip_equal_inds, c)
                end
            end
            # TODO remove i!=j from vs here?
            idx_subs = Iterators.product((1:c for c in borders)...)
            for sub in idx_subs
                check_skip = false
                for (ci,csub) in zip(idx_combs,combinations(sub,2))
                    if csub[1]==csub[2] && ci in skip_equal_inds
                        check_skip = true
                    end
                    check_skip && break
                end
                check_skip && continue
                v_ = swap_index(vs[i], idx[1], sub[1])
                (v_' in vs_) && continue
                r_ = swap_index(rhs_exp[i], idx[1], sub[1])
                for jj=2:length(idx)
                    v_ = swap_index(v_, idx[jj], sub[jj])
                    r_ = swap_index(r_, idx[jj], sub[jj])
                end
                push!(rhs_, r_)
                push!(vs_, v_)
            end
        end
    end
    return vs_, rhs_
end

expand_sums(x, idx_borders) = x
function expand_sums(t::NumberTerm, idx_borders)
    if t.f === Sum
        idx = t.arguments[2:end]
        arg = t.arguments[1]
        # TODO remove i!=j terms when i==j
        if length(idx) > 1
            s = Sum(arg, idx[2:end]...)
            n = idx[1].count
            i = findfirst(x->isequal(x[1],n),idx_borders)
            s_ = expand_sums(s, idx_borders)
            return expand_sums(Sum(s_, idx[1]), idx_borders[i])
        end
        idx_sym = idx[1]
        border = if idx_borders isa Vector
            i = findfirst(x->isequal(x[1],idx_sym.count),idx_borders)
            idx_borders[i]
        else
            idx_borders
        end
        @assert isequal(idx_sym.count, border[1])
        args_ = Number[]
        for i in 1:border[2]
            ex = swap_index(arg, idx_sym, i)
            push!(args_, ex)
        end
        return +(args_...)
    else
        args = [expand_sums(arg, idx_borders) for arg in t.arguments]
        return t.f(args...)
    end
end

"""
    build_ode(eqs::DifferentialEquation, ps=[], usym=:u,
                psym=:p, tsym=:t; set_unknowns_zero::Bool=false, check_bounds::Bool=false)

From a set of differential equations`eqs` of averages, generate a `Meta.Expr`
containing the code for a function which can be directly passed to `OrdinaryDiffEq`
in order to solve it.

# Arguments
*`eqs::DifferentialEquation`: The set of (average) equations.
*`ps=[]`: List of symbolic parameters, which are parsed into parameters
    used in DiffEq functions.
*`usym=:u`: The symbol used for the variable vector.
*`psym=:p`: The symbol used for the parameter vector.
*`tsym=:t`: The symbol used for the time parameter.

# Optional arguments
*`check_bounds::Bool=false`: Choose whether the resulting function should contain
    the `@inbounds` flag, which skips bounds checking for performance.
"""
build_ode(eqs::DifferentialEquation, args...; kwargs...) = build_ode(eqs.rhs,eqs.lhs,args...;kwargs...)

"""
    generate_ode(eqs::DifferentialEquation, ps=[], usym=:u,
                psym=:p, tsym=:t; set_unknowns_zero::Bool=false, check_bounds::Bool=false)

From a set of differential equations `eqs` of averages, generate a `Function`
which can be directly used in `OrdinaryDiffEq`. Essentially, this calls `Meta.eval`
on the output of the `build_ode` function.

# Arguments
*`eqs::DifferentialEquation`: The set of (average) equations.
*`ps=[]`: List of symbolic parameters, which are parsed into parameters
    used in DiffEq functions.
*`usym=:u`: The symbol used for the variable vector.
*`psym=:p`: The symbol used for the parameter vector.
*`tsym=:t`: The symbol used for the time parameter.

# Optional arguments
*`check_bounds::Bool=false`: Choose whether the resulting function should contain
    the `@inbounds` flag, which skips bounds checking for performance.

# Related methods
    generate_ode(eqs::Vector, vs::Vector, ps=[], usym=:u,
                psym=:p, tsym=:t; set_unknowns_zero::Bool=false, check_bounds::Bool=false)
"""
generate_ode(args...;kwargs...) = Meta.eval(build_ode(args...;kwargs...))


# Auxiliary functions
function build_expr(head::Symbol, args)
    ex = Expr(head)
    append!(ex.args, args)
    ex
end

function throw_missing_error(missed)
    error_msg = "The following parameters or averages are missing: "
    for p1=missed
        error_msg *= "$p1 "
    end
    error(error_msg)
end
