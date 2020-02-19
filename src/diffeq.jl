import MacroTools

"""
    build_ode(eqs::Vector{<:SymPy.Sym}, vs::Vector{<:SymPy.Sym}, ps=[], usym=:u,
                psym=:p, tsym=:t; set_unknowns_zero::Bool=false, check_bounds::Bool=false)

From a set of equations contained in `eqs`, generate a `Meta.Expr` containing the
code for a function which can be directly passed to `OrdinaryDiffEq` in order to solve
it. The variable vector `u` corresponds to the symbols provided in `vs`.

# Arguments
*`eqs::Vector{<:SymPy.Sym}`: The vector containing the right-hand side of equations.
*`vs::Vector{<:SymPy.Sym}`: The vector containing the left-hand side of equations.
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
function build_ode(eqs::Vector{<:SymPy.Sym}, vs::Vector{<:SymPy.Sym}, ps=[], usym=:u, psym=:p, tsym=:t;
                    set_unknowns_zero::Bool=false, check_bounds::Bool=false)
    @assert length(eqs) == length(vs)

    # Check if there are unknown symbols
    missed = check_missing(eqs,vs,ps)
    (isempty(missed) || set_unknowns_zero) || throw_missing_error(missed)
    eqs = remove_unknowns(eqs,missed)

    # Substitute using SymPy
    dusym = Symbol(string("d",usym))
    lhs = [:($dusym[$i]) for i=1:length(eqs)]
    u = [SymPy.symbols("$usym[$i]") for i=1:length(vs)]
    p = [SymPy.symbols("$psym[$i]") for i=1:length(ps)]
    subs_u = Dict(vs .=> u)
    subs_p = Dict(ps .=> p)
    subs = merge(subs_u, subs_p)
    rhs = [eq(subs) for eq=eqs]


    # From https://github.com/JuliaDiffEq/ModelingToolkit.jl/blob/dca5f38491ae6dea431cb2a7cceb055645086034/src/utils.jl#L44
    rhs_sym = parse_sympy(rhs)
    line_eqs = [Expr(:(=), lhs[i], rhs_sym[i]) for i=1:length(lhs)]
    var_eqs = build_expr(:block, line_eqs)
    var_eqs = MacroTools.postwalk(ex -> ex == :I ? :im : ex, var_eqs)
    var_eqs = MacroTools.postwalk(ex -> ex == :Dagger ? :adjoint : ex, var_eqs)
    var_eqs = MacroTools.postwalk(ex -> ex == :conjugate ? :conj : ex, var_eqs)

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

function build_indexed_ode(eqs::Vector{<:SymPy.Sym}, vs::Vector{<:SymPy.Sym}, ps=[], usym=:u, psym=:p, tsym=:t;
                    set_unknowns_zero::Bool=false, check_bounds::Bool=false)
    @assert length(eqs) == length(vs)

    # Check if there are unknown symbols
    missed = check_missing(eqs,vs,ps)
    (isempty(missed) || set_unknowns_zero) || throw_missing_error(missed)
    eqs_ = remove_unknowns(eqs,missed)

    # Check for indices
    rhs_inds = Vector{Index}[]
    lhs_inds = Vector{Index}[]
    has_index_lhs = zeros(Bool, length(vs))
    has_index_rhs = zeros(Bool, length(eqs))
    for i=1:length(vs)
        r_inds, _ = find_index(eqs[i])
        if !isempty(r_inds)
            has_index_rhs[i] = true
            push!(rhs_inds, r_inds)
        end

        l_inds, _ = find_index(vs[i])
        if !isempty(l_inds)
            has_index_lhs[i] = true
            push!(lhs_inds, l_inds)
        end
    end
    @assert (any(has_index_lhs) || any(has_index_rhs))

    # Sort equations by number of indices on the lhs (loop depth)
    sp = sortperm(has_index_lhs)
    eqs_ = eqs_[sp]
    vs_ = vs[sp]
    lhs_inds = lhs_ind[sp]
    has_index_lhs = has_index_lhs[sp]
    rhs_inds = rhs_ind[sp]
    has_index_rhs = has_index_rhs[sp]


    # Loop depth
    nloop = length.(lhs_ind)
    # unique!(nloop)
    i_loop = [Symbol("i$k") for k=nloop]
    # (nloop[end] > 1) && error("Nested loops not yet implemented!")

    n_no_index = sum(.!(has_index_lhs))
    u_index = Any[1:length(n_no_index);]
    for i=1:length(vs)-n_no_index
        push!(u_index, i_loop)
    end

    # Substitute using SymPy
    lhs = [Symbol(string("d",usym,"[$i]")) for i=u_index]
    u = [SymPy.symbols("$usym[$i]") for i=u_index]
    p = [SymPy.symbols("$psym[$i]") for i=1:length(ps)]
    subs_u = Dict(vs .=> u)
    subs_p = Dict(ps .=> p)
    subs = merge(subs_u, subs_p)
    rhs = [eq(subs) for eq=eqs]

    # Replacement for arguments of symbolic sums
    for i=1:length(eqs)
        for (keys,vals)=(subs)
            # Check for indexed objects; others have already been replaced
            if classname(keys)=="Indexed"
                k_inds = keys.__pyobject__.indices
                length(k_inds)==1 || continue # TODO: other cases
                k_ = IndexOrder[findfirst(x->sympify(x)==k_inds[1],IndexOrder)]
                # Replace indices by any other known index and try to substitute in expression
                for j=IndexOrder
                    key_ = swap_index(keys, k_, j)
                    val_ = swap_index(vals, k_, j)
                    eqs[i] = eqs[i].__pyobject__.replace(key_,val_)
                end
            end
        end
    end

    # From https://github.com/JuliaDiffEq/ModelingToolkit.jl/blob/dca5f38491ae6dea431cb2a7cceb055645086034/src/utils.jl#L44
    rhs_sym = parse_sympy(rhs)
    line_eqs = [Expr(:(=), lhs[i], rhs_sym[i]) for i=1:length(lhs)]
    var_eqs = build_expr(:block, line_eqs)
    var_eqs = MacroTools.postwalk(ex -> ex == :I ? :im : ex, var_eqs)
    var_eqs = MacroTools.postwalk(ex -> ex == :Dagger ? :adjoint : ex, var_eqs)
    var_eqs = MacroTools.postwalk(ex -> ex == :conjugate ? :conj : ex, var_eqs)
    var_eqs = MacroTools.postwalk(ex -> ex == :Sum ? :sum : ex, var_eqs)

    # for i=IndexOrder
    #     # TODO: replace SymPy sums by actual sums with a for loop
    #     var_eqs = MacroTools.postwalk(x -> @capture(x, Sum_(xs__)) ?

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
function build_ode(eqs::Vector{<:DifferentialEquation}, args...;kwargs...)
    lhs = [e.lhs for e=eqs]
    rhs = [e.rhs for e=eqs]
    return build_ode(rhs,lhs,args...;kwargs...)
end

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
    generate_ode(eqs::Vector{<:SymPy.Sym}, vs::Vector{<:SymPy.Sym}, ps=[], usym=:u,
                psym=:p, tsym=:t; set_unknowns_zero::Bool=false, check_bounds::Bool=false)
"""
generate_ode(args...;kwargs...) = Meta.eval(build_ode(args...;kwargs...))


"""
    check_missing(rhs::Vector, vs::Vector, ps=[])

For a list of expressions contained in `rhs`, check whether all occurring symbols
are contained either in the variables given in `vs` or `ps`. Returns a list of
missing symbols.
"""
function check_missing(eqs::Vector{<:SymPy.Sym}, vs::Vector{<:SymPy.Sym}, ps=[])
    missed = typejoin(SymPy.Sym,eltype(ps))[]
    for e=eqs
        append!(missed,SymPy.free_symbols(e))
    end
    unique!(missed)
    filter!(x->!(x∈vs || x∈ps),missed)
    return missed
end

"""
    check_missing(de::DifferentialEquationSet,ps=[])

For a set of (average) differential equations described by `de`, check whether
all symbols occurring on the right-hand side are either contained in the
left-hand side or in the parameters given in `ps`. In other words, this function
check whether a set of equations is complete.

# Arguments
*`de::DifferentialEquationSet`: The set of differential equations.
*`ps=[]`: The list of parameters which occur in the equations.
"""
check_missing(de::DifferentialEquationSet,ps=[]) = check_missing(de.rhs,de.lhs,ps)

"""
    remove_unknowns(eqs::Vector,unknowns::Vector)

Substitute all `unknowns` in the equations in `eqs` by zero.
"""
function remove_unknowns(eqs::Vector,unknowns::Vector)
    subs = Dict(unknowns .=> 0)
    return [e(subs) for e=eqs]
end

"""
    remove_unknowns(de::DifferentialEquationSet,ps=[])

Substitute all symbols that occur on the right-hand side of the set of equations
in `de`, but are not contained in the left-hand side or the parameters `ps` by zero.
This function uses the `check_missing` function to find any unknown symbols.

# Arguments
*`de::DifferentialEquationSet`: The set of differential equations.
*`ps=[]`: The list of parameters which occur in the equations.
"""
function remove_unknowns(de::DifferentialEquationSet,ps=[])
    missed = check_missing(de,ps)
    rhs = remove_unknowns(de.rhs,missed)
    lhs = remove_unknowns(de.lhs,missed)
    lhs_ = eltype(lhs)[]
    rhs_ = eltype(rhs)[]
    for (l,r)=zip(lhs,rhs)
        iszero(l) || (push!(lhs_,l); push!(rhs_,r))
    end
    return DifferentialEquationSet(lhs_,rhs_)
end

# Auxiliary functions
function build_expr(head::Symbol, args)
    ex = Expr(head)
    append!(ex.args, args)
    ex
end

# TODO: can this string parsing be avoided?
parse_sympy(args) = Meta.parse.(string.(args))

function throw_missing_error(missed)
    error_msg = "The following symbols (either parameters or averages) are missing: "
    for p1=missed
        error_msg *= "$p1 "
    end
    error_msg *= "\n"
    error_msg *= "If you want to neglect those, set the `set_unknowns_zero` kwarg to `true`.\n"
    error(error_msg)
end
