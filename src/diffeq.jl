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
function build_ode(eqs, vs, args...; kwargs...)
    if any(x->(classname(x)=="Indexed"), vs) #TODO: check RHS for indexes
        return build_indexed_ode(eqs, vs, args...; kwargs...)
    else
        _build_ode(eqs, vs, args...; kwargs...)
    end
end
function _build_ode(eqs::Vector{<:SymPy.Sym}, vs::Vector{<:SymPy.Sym}, ps=[], usym=:u, psym=:p, tsym=:t;
                    set_unknowns_zero::Bool=false, check_bounds::Bool=true)
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
                    set_unknowns_zero::Bool=false, check_bounds::Bool=true)
    @assert length(eqs) == length(vs)

    # Check if there are unknown symbols
    missed = check_missing_idx(eqs,vs,ps)
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
    has_index_lhs = has_index_lhs[sp]
    has_index_rhs = has_index_rhs[sp]

    # Loop depth
    nloop = length.(lhs_inds)

    # Loop offset
    offset0 = length(vs_)-length(lhs_inds)

    # Loop boundaries
    lower = [[l_.lower for l_=l] for l=lhs_inds]
    upper = [[u_.upper for u_=l] for l=lhs_inds]

    # Generate linear indices for u
    lower_lin, upper_lin = _gen_lin_inds(lower,upper,offset0)
    offset = cumsum([offset0; upper_lin .- lower_lin .+ 1])[1:end-1]

    idxMap = _gen_idxMap(lower,upper,offset)

    i_loop = [SymPy.symbols.([l.label for l=l1],integer=true) for l1=lhs_inds]

    n_no_index = sum(.!(has_index_lhs))
    u_index = Any[1:length(n_no_index);]
    append!(u_index, [idxMap(k,i_loop[k]) for k=1:length(i_loop)])

    # Substitute using SymPy
    # TODO: substitute ⟨.⟩ with some variables losing the ⟨,⟩, parse, and use MacroTools from here on out
    dusym = Symbol(string("d",usym))
    lhs = [:($dusym[$(i)]) for i=u_index]
    u_sym_base = SymPy.sympy.IndexedBase("$usym")
    u = [u_sym_base[i] for i=u_index]

    p = [SymPy.symbols("$psym[$i]") for i=1:length(ps)]
    subs_u = Dict(vs_ .=> u)
    subs_p = Dict(ps .=> p)
    subs = merge(subs_u, subs_p)
    rhs = [eq(subs) for eq=eqs_]

    # Replacement for arguments of symbolic sums
    for i=1:length(eqs)#n_no_index
        for uu=1:length(i_loop)
            keys = vs_[n_no_index+uu]
            # Check for indexed objects; others have already been replaced
            if classname(keys)=="Indexed"
                k_inds = keys.__pyobject__.indices
                k_ = [IndexOrder[findfirst(x->sympify(x)==kk,IndexOrder)] for kk=k_inds]
                inds_combs = combinations(IndexOrder,length(k_))
                # Replace indices by any other known index and try to substitute in expression
                for j=inds_combs
                    key_ = keys
                    for kk=1:length(k_)
                        key_ = swap_index(key_, k_[kk], j[kk])
                    end
                    # Get value index from idxMap
                    jval = idxMap(uu,sympify.(j))
                    val_ = u_sym_base[jval]
                    rhs[i] = rhs[i].__pyobject__.replace(key_,val_)
                end
            end
        end
    end
    rhs_sym = parse_sympy(rhs)

    for ii=1:length(rhs_sym)
        # Replace SymPy Sums by actual sums
        rhs_sym[ii] = MacroTools.postwalk(x -> MacroTools.@capture(x, Sumsym_(arg_, (i_,l_,u_))) ?
                        :( sum($arg for $(i)=$(l):$(u)) ) : x,
                            rhs_sym[ii])

        # Account for neq indices
        # Replace indexing of u with may have an offset
        rhs_sym[ii] = MacroTools.postwalk(x -> MacroTools.@capture(x, y_[c1_*c2_ + i_≠j_+off_]) ? :(Int($i ≠ $j) * $y[$c1*$c2 + $i+$off]) : x, rhs_sym[ii])
        rhs_sym[ii] = MacroTools.postwalk(x -> MacroTools.@capture(x, y_[i_≠j_+off_]) ? :(Int($i ≠ $j) * $y[$i+$off]) : x, rhs_sym[ii])
        # Replace remaining indices (of parameters) without offset
        rhs_sym[ii] = MacroTools.postwalk(x -> MacroTools.@capture(x, y_[i_≠j_]) ? :(Int($i ≠ $j) * $y[$i]) : x, rhs_sym[ii])
        # Replace != in sum loop iteration
        rhs_sym[ii] = MacroTools.postwalk(x -> MacroTools.@capture(x, i_ ≠ j_ = l_ : u_) ? :($(i)=$(l):$(u)) : x, rhs_sym[ii])
    end
    return rhs_sym

    loop_eqs = [Expr(:(=), lhs[i], rhs_sym[i]) for i=n_no_index+1:length(lhs)]
    loop_block = build_expr(:block, loop_eqs)

    # TODO: nested loops
    loop_ex = :(
        for $(i_loop[1])=$(lower[1]) : $(upper[1])
            $loop_block
        end
    )

    # Non-indexed lines
    line_eqs = [Expr(:(=), lhs[i], rhs_sym[i]) for i=1:n_no_index]
    var_eqs = build_expr(:block, [line_eqs;loop_ex])
    var_eqs = MacroTools.postwalk(ex -> ex == :I ? :im : ex, var_eqs)
    var_eqs = MacroTools.postwalk(ex -> ex == :Dagger ? :adjoint : ex, var_eqs)
    var_eqs = MacroTools.postwalk(ex -> ex == :conjugate ? :conj : ex, var_eqs)

    # TODO: cleaner solution when setting unknowns zero
    # Remove (0)[j]
    var_eqs = MacroTools.postwalk(x -> MacroTools.@capture(x, (0)[j_]) ? 0 : x, var_eqs)

    # TODO: clean up
    # Replace loop borders
    if eltype(lower) <: Symbol
        for l=lower
            # Find substitute
            for s=subs_p
                if l==Symbol(s[1]+1)
                    # TODO: avoid parsing
                    s_ = Meta.parse(string(s[2]))
                    var_eqs = MacroTools.postwalk(x -> x==l ? s_ : x, var_eqs)
                end
            end
        end
    end
    if eltype(upper) <: Symbol
        for l=upper
            # Find substitute
            for s=subs_p
                if l==Symbol(s[1]+1)
                    # TODO: avoid parsing
                    s_ = Meta.parse(string(s[2]))
                    var_eqs = MacroTools.postwalk(x -> x==l ? s_ : x, var_eqs)
                end
            end
        end
    end

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

function _gen_lin_inds(lower,upper,offset)
    n = length.(lower)

    # Convert symbolics to Sympy
    lower_ = [eltype(l)<:Symbol ? SymPy.symbols.(l,integer=true) : l for l=lower]
    upper_ = [eltype(u)<:Symbol ? SymPy.symbols.(u,integer=true) : u for u=upper]

    # Size of each index space
    nsize = [prod([lower_[i][j] - 1 + upper_[i][j] for j=1:n[i]]) for i=1:length(n)]

    # Compute corresponding linear index boundaries
    upper_lin = offset .+ cumsum(nsize)
    lower_lin = upper_lin .- nsize .+ 1

    return lower_lin, upper_lin
end

function _gen_idxMap(lower,upper,offset)
    n = length.(lower)

    # Convert symbolics to Sympy
    lower_ = [eltype(l)<:Symbol ? SymPy.symbols.(l,integer=true) : l for l=lower]
    upper_ = [eltype(u)<:Symbol ? SymPy.symbols.(u,integer=true) : u for u=upper]

    # Size of each index space
    nsize = [[lower_[i][j] - 1 + upper_[i][j] for j=1:n[i]] for i=1:length(n)]

    # For each index space k, return a linear index from a set of indices of length corresponding to k
    _idxMap(k,inds) = (offset[k] .+ _linear_index(inds,nsize[k]))
    return _idxMap
end
function _linear_index(inds,nsize)
    # TODO: test for dims > 3
    strds = Base.size_to_strides(1, nsize...)
    inds_ = inds .- 1
    return sum(inds_ .* strds) + 1
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
function check_missing_idx(eqs::Vector{<:SymPy.Sym}, vs::Vector{<:SymPy.Sym}, ps=[])
    missed = typejoin(SymPy.Sym,eltype(ps))[]
    for e=eqs
        append!(missed,SymPy.free_symbols(e))
    end
    unique!(missed)
    filter!(x->!(classname(x)=="Indexed" || classname(x)=="Idx"),missed)

    vars = eltype(vs)[]
    for v=vs
        append!(vars,SymPy.free_symbols(v))
    end
    unique!(vars)

    pars = if eltype(ps) <: SymPy.Sym
        pars_ = SymPy.Sym[]
        for p=ps
            append!(pars_,SymPy.free_symbols(p))
        end
        pars_
    else
        ps
    end

    filter!(x->!(x∈vars || x∈pars),missed)
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
