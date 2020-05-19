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
*`set_unknowns_zero::Bool=false`: Choose whether encountered symbols which are not
    contained in either `vs` or `ps` should be neglected (set to 0).
*`check_bounds::Bool=false`: Choose whether the resulting function should contain
    the `@inbounds` flag, which skips bounds checking for performance.
"""
function build_ode(rhs::Vector, vs::Vector, ps=[], usym=:u, psym=:p, tsym=:t;
                    set_unknowns_zero::Bool=false, check_bounds::Bool=false)
    @assert length(rhs) == length(vs)

    vs_adj_ = average_adj.(vs)

    # Check if there are unknown symbols
    missed = find_missing(eqs,vs;vs_adj=vs_adj_,ps=ps)
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

"""
    build_ode(eqs::DifferentialEquationSet, ps=[], usym=:u,
                psym=:p, tsym=:t; set_unknowns_zero::Bool=false, check_bounds::Bool=false)

From a set of differential equations`eqs` of averages, generate a `Meta.Expr`
containing the code for a function which can be directly passed to `OrdinaryDiffEq`
in order to solve it.

# Arguments
*`eqs::DifferentialEquationSet`: The set of (average) equations.
*`ps=[]`: List of parameters (possibly `SymPy.Sym`s), which are parsed into parameters
    used in DiffEq functions.
*`usym=:u`: The symbol used for the variable vector.
*`psym=:p`: The symbol used for the parameter vector.
*`tsym=:t`: The symbol used for the time parameter.

# Optional arguments
*`set_unknowns_zero::Bool=false`: Choose whether encountered symbols which are not
    contained in either `vs` or `ps` should be neglected (set to 0).
*`check_bounds::Bool=false`: Choose whether the resulting function should contain
    the `@inbounds` flag, which skips bounds checking for performance.
"""
build_ode(eqs::DifferentialEquationSet, args...; kwargs...) = build_ode(eqs.rhs,eqs.lhs,args...;kwargs...)

"""
    generate_ode(eqs::DifferentialEquationSet, ps=[], usym=:u,
                psym=:p, tsym=:t; set_unknowns_zero::Bool=false, check_bounds::Bool=false)

From a set of differential equations `eqs` of averages, generate a `Function`
which can be directly used in `OrdinaryDiffEq`. Essentially, this calls `Meta.eval`
on the output of the `build_ode` function.

# Arguments
*`eqs::DifferentialEquationSet`: The set of (average) equations.
*`ps=[]`: List of parameters (possibly `SymPy.Sym`s), which are parsed into parameters
    used in DiffEq functions.
*`usym=:u`: The symbol used for the variable vector.
*`psym=:p`: The symbol used for the parameter vector.
*`tsym=:t`: The symbol used for the time parameter.

# Optional arguments
*`set_unknowns_zero::Bool=false`: Choose whether encountered symbols which are not
    contained in either `vs` or `ps` should be neglected (set to 0).
*`check_bounds::Bool=false`: Choose whether the resulting function should contain
    the `@inbounds` flag, which skips bounds checking for performance.

# Related methods
    generate_ode(eqs::Vector, vs::Vector, ps=[], usym=:u,
                psym=:p, tsym=:t; set_unknowns_zero::Bool=false, check_bounds::Bool=false)
"""
generate_ode(args...;kwargs...) = Meta.eval(build_ode(args...;kwargs...))


"""
    find_missing(rhs::Vector, vs::Vector, vs_adj=adjoint.(vs) ps=[])

For a list of expressions contained in `rhs`, check whether all occurring symbols
are contained either in the variables given in `vs` or `ps`. Returns a list of
missing symbols.
"""
function find_missing(rhs::Vector, vs::Vector; vs_adj::Vector=average_adj.(vs), ps=[])
    missed = Number[]
    for e=rhs
        append!(missed,get_symbolics(e))
    end
    unique!(missed)
    filter!(x->!(x∈vs || x∈ps || x∈vs_adj),missed)
    return missed
end
find_missing(de::DifferentialEquationSet; kwargs...) = find_missing(de.rhs, de.lhs; kwargs...)

# Auxiliary functions
get_symbolics(x::Number) = Number[]
get_symbolics(x::SymbolicNumber) = [x]
function get_symbolics(t::NumberTerm)
    syms = Number[]
    for arg in t.arguments
        append!(syms, get_symbolics(arg))
    end
    return unique(syms)
end

function average_adj(avg::Average)
    # Adjoint reverses order in product; need to simplify to re-order
    # TODO: improve
    op_adj = simplify_operators(adjoint(avg.operator))
    return average(op_adj)
end

function build_expr(head::Symbol, args)
    ex = Expr(head)
    append!(ex.args, args)
    ex
end

function throw_missing_error(missed)
    error_msg = "The following symbols (either parameters or averages) are missing: "
    for p1=missed
        error_msg *= "$p1 "
    end
    error_msg *= "\n"
    error_msg *= "If you want to neglect those, set the `set_unknowns_zero` kwarg to `true`.\n"
    error(error_msg)
end
