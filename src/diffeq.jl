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
    # isempty(missed) || throw_missing_error(missed)

    _vs, _rhs = _expand_indexed(vs, rhs, idx_borders)
    vs_adj_ = adjoint.(_vs)

    dusym = Symbol(:d,usym)
    us = [:($usym[$i]) for i=1:length(_vs)]
    dus = [:($dusym[$i]) for i=1:length(_vs)]

    vs_ = _to_expression.(_vs)
    vs_adj = _to_expression.(vs_adj_)
    rhs_ = [_expand_sums(r, idx_borders, vs_, vs_adj, us) for r in _rhs]
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

function _expand_indexed(vs::Vector, rhs::Vector, idx_borders::Vector)
    vs_ = Number[]
    rhs_ = Number[]
    for i=1:length(vs)
        idx = find_index(vs[i])
        if isempty(idx)
            push!(vs_, vs[i])
        else
            counts = getfield.(idx, :count)
            borders = [idx_borders[findfirst(x->isequal(x[1],c), idx_borders)][2] for c in counts]
            v_, r_ = _expand_indexed(vs[i], rhs[i], idx, borders)
            append!(vs_, v_)
            append!(rhs_, r_)
        end
    end
    vs_unq = unique_ops(vs_)
    rhs_unq = Number[]
    for j=1:length(vs_)
        if vs_[j] in vs_unq
            push!(rhs_unq, rhs_[j])
        end
    end
    return vs_unq, rhs_unq
end
_expand_indexed(vs, rhs, idx_borders::Pair) = _expand_indexed(vs, rhs, [idx_borders])
function _expand_indexed(v, r, idx, borders)
    idx_combs = combinations(idx,2)
    skip_equal_inds = []
    for c in idx_combs
        if has_expr(c[1] != c[2], v)
            push!(skip_equal_inds, c)
        elseif has_expr(c[2] != c[1], v)
            push!(skip_equal_inds, c)
        end
    end

    vs_ = Number[]
    rs_ = Number[]
    for sub in Iterators.product((1:n for n in borders)...)
        check_skip = false
        for (ci,csub) in zip(idx_combs,combinations(sub,2))
            if csub[1]==csub[2]
                for sk in skip_equal_inds
                    check_skip = all(isequal.(ci, sk))
                    check_skip && break
                end
            end
            check_skip && break
        end
        check_skip && continue
        v_ = swap_index(v, idx[1], sub[1])
        r_ = swap_index(r, idx[1], sub[1])
        for j=2:length(idx)
            v_ = swap_index(v_, idx[j], sub[j])
            r_ = swap_index(r_, idx[j], sub[j])
        end
        push!(vs_, v_)
        push!(rs_, r_)
    end
    if isempty(skip_equal_inds)
        return vs_, rs_
    else # Remove things such as 2 != 1 from LHS; TODO: optimize
        return simplify_constants.(vs_), rs_
    end
end

function _expand_sums(r::NumberTerm, idx_borders, vs_, vs_adj, us)
    if r.f === Sum
        return _expand_sum(r.arguments, idx_borders, vs_, vs_adj, us)
    else
        args = [_expand_sums(arg, idx_borders, vs_, vs_adj, us) for arg in r.arguments]
        return Expr(:call, Symbol(r.f), args...)
    end
end
_expand_sums(r, idx_borders, vs_, vs_adj, us) = _to_expression(r)

function _expand_sum(args::Vector, idx_borders, vs_, vs_adj, us)
    s_idx = args[2:end]
    counts = getfield.(s_idx, :count)
    borders = [idx_borders[findfirst(x->isequal(x[1],c), idx_borders)][2] for c in counts]

    args_ = []
    for sub in Iterators.product((1:n for n in borders)...)
        arg_ = swap_index(args[1], s_idx[1], sub[1])
        for j=2:length(s_idx)
            arg_ = swap_index(arg_, s_idx[j], sub[j])
        end
        push!(args_, _to_expression(arg_))
    end

    function _pw_func(x)
        if x in vs_
            i = findfirst(isequal(x),vs_)
            return us[i]
        elseif x in vs_adj
            i = findfirst(isequal(x),vs_adj)
            return :( conj($(us[i])) )
        elseif MacroTools.@capture(x, IndexedParameter(p_, idx_))
            return :(($p)[$(idx...)])
        # elseif MacroTools.@capture(x, SUM_PW(args_))
        #     idx = find_meta_idx.(args)
        #     n = length(idx[1])
        #     @assert all(length.(idx) .== n)
        #     ex_idx = (Symbol.('i'+k for k=0:n-1)...,)
        #     arg = set_meta_idx(args[1], ex_idx, 1)
        #     return :( sum($arg for $ex_idx in zip($idx)) )
        #     # return _expand_sum(arg, idx, vs_, vs_adj, idx_borders)
        else
            return x
        end
    end

    args_pw = [MacroTools.postwalk(_pw_func, arg) for arg in args_]
    idx = find_meta_idx.(args_pw)
    n = sum(length.(idx[1]))
    @assert all([sum(length.(i))==n for i in idx])
    idx_meta = Symbol.(('i':'z')[1:n])
    args_sum = Expr[set_meta_idx(args_pw[1], idx_meta)]
    for arg in args_pw
        arg_ = set_meta_idx(arg, idx_meta)
        arg_ in args_sum || push!(args_sum, arg_)
    end

    ex_args = [:( sum($arg for ($(idx_meta...),) in [$(idx...)]) ) for arg in args_sum]
    if length(ex_args)==1
        return ex_args[1]
    else
        return Expr(:call, :+, ex_args...)
    end

    # return :( sum([$(args_...)]) )
    #
    # itr_idx = Tuple{Vararg{Int,length(s_idx)}}[]
    # idx_expr = [Symbol('i'+j) for j=0:length(arg_idx)-1]
    #
    # function _pw_func(x)
    #     if x in vs_
    #         i = findfirst(isequal(x),vs_)
    #         return usym
    #     elseif x in vs_adj
    #         i = findfirst(isequal(x),vs_adj)
    #         return :( conj($(us[i])) )
    #     elseif MacroTools.@capture(x, IndexedParameter(p_, idx_))
    #         return :(($p)[$(idx...)])
    #     elseif MacroTools.@capture(x, Sum(arg_, idx_))
    #         # return _expand_sum(arg, idx, vs_, vs_adj, idx_borders)
    #     else
    #         return x
    #     end
end

function find_meta_idx(ex)
    idx = Int[]
    idx_ = _find_meta_idx(ex, idx)
    !isnothing(idx_) || error("Could not determine Index in Expr $ex")
    return idx_
end
function _find_meta_idx(ex::Expr, idx)
    if ex.head === :ref || ex.args[1] === :!=
        append!(idx, (ex.args[2:end]...,))
    else
        for arg in ex.args
            _find_meta_idx(arg, idx)
        end
    end
    return idx
end
_find_meta_idx(x,idx) = nothing

function set_meta_idx(ex, idx)
    count = 1
    ex_ = copy(ex)
    count, ex_ = _set_meta_idx!(ex_, idx, count)
    return ex_
end
function _set_meta_idx!(ex::Expr, idx, count)
    if ex.head === :ref || ex.args[1] === :!=
        for i=2:length(ex.args)
            ex.args[i] = idx[count-2+i]
        end
        count += length(ex.args)-1
        return count, ex
    else
        args_ = []
        for arg in ex.args
            count, ex_ = _set_meta_idx!(arg, idx, count)
            push!(args_, ex_)
        end
        return count, Expr(ex.head, args_...)
    end
end
_set_meta_idx!(ex, idx, count) = count, ex

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
