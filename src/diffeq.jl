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

    # TODO: cleaner solution to get rid of 0[j]
    ind_len = unique!(length.(lhs_inds))
    inds_combs_ = [combinations(IndexOrder,k) for k=ind_len]
    zr_sym = SymPy.Sym(0)
    for ii=inds_combs_
        for i=ii
            eqs_ = [SymPy.expand(e.__pyobject__.replace(zr_sym[i...],0)) for e=eqs_]
        end
    end

    # Sort equations by number of indices on the lhs (loop depth)
    sp = sortperm(has_index_lhs)
    eqs_ = eqs_[sp]
    vs_ = vs[sp]
    has_index_lhs = has_index_lhs[sp]
    has_index_rhs = has_index_rhs[sp]

    # Loop offset
    offset0 = length(vs_)-length(lhs_inds)

    # Loop boundaries
    lower = [[l_.lower for l_=l] for l=lhs_inds]
    upper_ = [[u_.upper for u_=l] for l=lhs_inds]
    n_neq = [length.([j1.nid for j1=j]) for j=lhs_inds]
    # Need to eliminate one loop step for each neq index
    # TODO: change loop boundaries to skip neq indices
    upper = upper_#[]
    # for i=1:length(upper_)
    #     tmp_ = []
    #     for j=1:length(upper_[i])
    #         if iszero(n_neq[i][j])
    #             push!(tmp_,upper_[i][j])
    #         else
    #             u_ = isa(upper_[i][j],Symbol) ? :($(upper_[i][j]) - $(n_neq[i][j])) : upper_[i][j]-n_neq[i][j]
    #             push!(tmp_,u_)
    #         end
    #     end
    #     push!(upper,tmp_)
    # end

    # Generate linear indices for u
    idxMap = _gen_idxMap(lower,upper_,n_neq,offset0)
    i_loop = [SymPy.symbols.([l.label for l=l1],integer=true) for l1=lhs_inds]
    n_no_index = sum(.!(has_index_lhs))
    # TODO: avoid parsing
    u_index = [1:n_no_index;Meta.parse.(string.([idxMap(k,i_loop[k]) for k=1:length(i_loop)]))]

    # Substitute using SymPy
    # TODO: substitute ⟨.⟩ with some variables losing the ⟨,⟩, parse, and use MacroTools from here on out
    dusym = Symbol(string("d",usym))
    lhs = [:($dusym[$(i)]) for i=u_index]
    # TODO: intrinsically skip indices where indices cannot be equal
    lhs_check = []
    for j=1:length(i_loop)
        nid = [l.nid for l=lhs_inds[j]]
        nid_ind = findall(!isempty,nid)
        if isempty(nid_ind)
            push!(lhs_check, nothing)
        else
            inds = [l.label for l=lhs_inds[j]]
            inds1 = Symbol[]
            inds2 = Symbol[]
            for i=1:length(inds)
                if i∈nid_ind
                    inds_neq = IndexOrder[findall(x->x.id∈nid[i],IndexOrder)]
                    inds_neq_sym = unique([ii.label for ii=inds_neq])
                    append!(inds2,inds_neq_sym)
                    append!(inds1,[inds[i] for ii=1:length(inds_neq_sym)])
                end
            end
            check_ex = :( (( $(inds1[1]) != $(inds2[1]))) )
            for jj=2:length(inds1)
                check_ex = :( $check_ex && ($(inds1[jj]) != $(inds2[jj])) )
            end
            push!(lhs_check, check_ex)
        end
    end
    u = [SymPy.symbols("$usym[$i]") for i=1:n_no_index]

    p = [SymPy.symbols("$psym[$i]") for i=1:length(ps)]
    subs_u = Dict(vs_[1:n_no_index] .=> u)
    subs_p = Dict(ps .=> p)
    subs = merge(subs_u, subs_p)
    rhs = [eq(subs) for eq=eqs_]

    # Symbolic sympy function; place-holder for bool-check of neq indices; can use sympy simplification to reduce number of bools computed
    # TODO: unique symbol, which can still be parsed
    bool_sym = SymPy.symbols("_BOOL_VAR")
    bool_fun = SymPy.Function(bool_sym)
    u_sym_base = SymPy.sympy.IndexedBase("$usym")
    # Replacement for arguments of symbolic sums
    tmp_keys = SymPy.Sym[]
    for i=1:length(eqs)#n_no_index
        for uu=1:length(i_loop)
            keys = vs_[n_no_index+uu]
            # Check for indexed objects; others have already been replaced
            if classname(keys)=="Indexed"
                k_inds = keys.__pyobject__.indices
                k_ = [IndexOrder[findfirst(x->sympify(x)==kk,IndexOrder)] for kk=k_inds]

                # Replace indices by any other known index and try to substitute in expression
                inds_combs = combinations(IndexOrder,length(k_))
                key_base = keys.__pyobject__.base
                for j_=inds_combs
                    for j=permutations(j_)
                        key_ = key_base[sympify.(j)...]
                        jsym = SymPy.symbols.([j1.label for j1=j],integer=true)
                        # Get value index from idxMap
                        jval = idxMap(uu,jsym)
                        _neq = [j1.nid for j1=j]
                        _check = 1
                        for nn=1:length(j)
                            length(_neq[nn]) > 0 || continue
                            _i_neq = sympify.(IndexOrder[findall(x->(x.id∈_neq[nn]&&isempty(x.nid)),IndexOrder)])
                            _j_nn = SymPy.symbols(j[nn].label)
                            for _ii=_i_neq
                                _check *= bool_fun(_j_nn,_ii)
                            end
                        end
                        val_ = _check*u_sym_base[jval]
                        if length(k_) > 1
                            push!(tmp_keys,key_)
                        end
                        rhs[i] = SymPy.expand(rhs[i].__pyobject__.replace(key_,val_))
                    end
                end
            end
        end
    end

    # Replace neq indices in ps; TODO: clean up; use postwalk
    for i=1:length(eqs)
        for pp=1:length(ps)
            keys = p[pp]
            for inds_combs=inds_combs_
                for j_=inds_combs
                    for j=permutations(j_)
                        _neq = [j1.nid for j1=j]
                        all(isempty.(_neq)) && continue
                        _check = 1
                        for nn=1:length(j)
                            length(_neq[nn]) > 0 || continue
                            _i_neq = sympify.(IndexOrder[findall(x->(x.id∈_neq[nn]&&isempty(x.nid)),IndexOrder)])
                            _j_nn = SymPy.symbols(j[nn].label)
                            for _ii=_i_neq
                                _check *= bool_fun(_j_nn,_ii)
                            end
                        end
                        jsym = SymPy.symbols.([j1.label for j1=j], integer=true)
                        key_ = keys[j...]
                        val_ = _check*SymPy.sympy.IndexedBase(keys)[jsym...]
                        rhs[i] = SymPy.expand(rhs[i].__pyobject__.replace(key_,val_))
                    end
                end
            end
        end
    end

    rhs_sym = parse_sympy(rhs)

    # Postwalk function to replace _BOOL_VAR(i,j) by Int(i!=j)
    function _pw_neq_func(x)
        if MacroTools.@capture(x, F_(i_, j_))
            if F == :_BOOL_VAR
                return :( Int($i != $j) )
            else
                return x
            end
        end
        return x
    end

    # Postwalk function to replace SymPy sums by actual sums
    function _pw_sums_func(x)
        if MacroTools.@capture(x, Sumsym_(arg_, (i_,l_,u_)))
            isa(arg,Number) && iszero(arg) && return 0
            if Sumsym==:Sum
                loop_index = MacroTools.@capture(i, j_ ≠ k_) ? j : i
                return :( sum($arg for $(loop_index)=$(l):$(u)) )
            else
                return x
            end
        end
        return x
    end

    for ii=1:length(rhs_sym)
        # Replace SymPy Sums by actual sums
        rhs_sym[ii] = MacroTools.postwalk(_pw_sums_func, rhs_sym[ii])

        # Account for neq indices
        rhs_sym[ii] = MacroTools.postwalk(_pw_neq_func, rhs_sym[ii])
    end

    nloop = Int[]
    loop_block_size = Int[]
    tmp_index = []
    for i=1:length(lhs_inds)
        if !(lhs_inds[i] ∈ tmp_index)
            push!(tmp_index, lhs_inds[i])
            push!(loop_block_size,1)
            push!(nloop,length(lhs_inds[i]))
        else
            loop_block_size[end] += 1
        end
    end

    # TODO: check whether one can partially combine loops
    # TODO: check limits only; replace index label if needed to combine
    can_combine_loops = all([l1∈l for l1=lhs_inds[1], l=lhs_inds])
    loop_eqs_ = [Expr(:(=), lhs[i], rhs_sym[i]) for i=n_no_index+1:length(lhs)]
    # Skip entries of neq indices
    for i=1:length(loop_eqs_)
        if !isa(lhs_check[i], Nothing)
            loop_eqs_[i] = :( ($(lhs_check[i])) && ($(loop_eqs_[i])) )
        end
    end

    if can_combine_loops
        _loop_eqs_ = reverse(loop_eqs_)

        loop_blocks_ = reverse(loop_block_size)
        nloop_ = reverse(nloop)
        lhs_inds_ = reverse(lhs_inds)
        _lower_ = reverse(lower)
        _upper_ = reverse(upper)
        loop_ex_ = Expr[]
        ind0 = 1
        for i=1:length(loop_block_size)
            size = loop_blocks_[i]
            depth = nloop_[i]
            if depth > length(lhs_inds[1])
                inds = lhs_inds_[i][findall(x->!(x∈lhs_inds[1]),lhs_inds_[i])]
            else
                inds = lhs_inds_[i]
            end
            if length(inds) == 1
                l = inds[1]
                iter_ex = :( $(l.label) = $(_lower_[i][1]) : $(_upper_[i][1]) )
            else
                exs = [:( $(inds[j].label) = $(_lower_[i][j]) : $(_upper_[i][j]) ) for j=1:length(inds)]
                iter_ex = build_expr(:block, exs)
            end
            l_eqs_ = build_expr(:block, [reverse(_loop_eqs_[ind0:size+ind0-1]);loop_ex_])
            ind0 += size
            isempty(loop_ex_) && push!(loop_ex_, build_expr(:for, [iter_ex, l_eqs_]))
            loop_ex_[1] = build_expr(:for, [iter_ex, l_eqs_])
        end
        loop_ex = build_expr(:block, loop_ex_)
    else
        loop_eqs = [build_expr(:block,[:( $l )]) for l=loop_eqs_]
        # iter_ex = [[:( $(l[i].label) = $(l[i].lower) : $(l[i].upper) ) for i=1:length(l)] for l=lhs_inds]
        iter_ex = Expr[]
        for i=1:length(lhs_inds)
            if length(lhs_inds[i]) == 1
                l = lhs_inds[i][1]
                ex_ = :( $(l.label) = $(_lower_[i][1]) : $(_upper_[i][1]) )
            else
                exs = [:( $(lhs_inds[i][j].label) = $(_lower_[i][j]) : $(_upper_[i][j]) ) for j=1:length(lhs_inds[i])]
                ex_ = build_expr(:block, exs)
            end
            push!(iter_ex, ex_)
        end
        loop_ex_ = [build_expr(:for, [iter_ex[i], loop_eqs[i]]) for i=1:length(loop_eqs)]
        loop_ex = build_expr(:block, loop_ex_)
    end

    # Non-indexed lines
    line_eqs = [Expr(:(=), lhs[i], rhs_sym[i]) for i=1:n_no_index]
    var_eqs = build_expr(:block, [line_eqs;loop_ex])
    var_eqs = MacroTools.postwalk(ex -> ex == :I ? :im : ex, var_eqs)
    var_eqs = MacroTools.postwalk(ex -> ex == :Dagger ? :adjoint : ex, var_eqs)
    var_eqs = MacroTools.postwalk(ex -> ex == :conjugate ? :conj : ex, var_eqs)

    # TODO: clean up
    # Replace loop borders
    borders = [collect(Iterators.flatten(lower)); collect(Iterators.flatten(upper))]
    for b=borders
        if isa(b, Symbol)
            # Find substitute
            p_ind = findfirst(x->Symbol(x+1)==b, ps)
            isa(p_ind, Nothing) && continue
            p_ = :($psym[$p_ind])
            var_eqs = MacroTools.postwalk(x->x==b ? p_ : x, var_eqs)
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

function _gen_idxMap(lower,upper,n_neq,offset0)
    n = length.(lower)

    # Convert symbolics to Sympy
    lower_ = [eltype(l)<:Symbol ? SymPy.symbols.(l,integer=true) : l for l=lower]
    upper_ = [eltype(u)<:Symbol ? SymPy.symbols.(u,integer=true) : u for u=upper]

    # Size of each index space
    # TODO: reduce size to account for neq indices
    nsize = [[lower_[i][j] - 1 + upper_[i][j] - 0n_neq[i][j] for j=1:n[i]] for i=1:length(n)]
    offset = cumsum([offset0;prod.(nsize)])

    n_neq_ = maximum.(n_neq)
    # For each index space k, return a linear index from a set of indices of length corresponding to k
    _idxMap(k,inds) = (offset[k] .+ _linear_index(inds,nsize[k]))# .- n_neq_[k]*_upper_[k])
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
