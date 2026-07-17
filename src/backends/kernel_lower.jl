# Lowering: completed MomentGraph -> moment-polynomial IR (issue #294, M·v design).
#
# The cumulant guarantee makes this total: every drift is a polynomial in the moments
# (degree <= order + interaction size - 1) with coefficients in parameters only. The IR is:
#
#   du = M * v
#
# where v is the vector of DISTINCT monomials over the states (hash-consed globally, so a
# product shared by many equations is computed once per RHS call), updated incrementally via
# prefix chains: each monomial = (parent monomial) * (one state factor). M is one sparse
# matrix whose values are the numerically evaluated coefficients.
#
# Encoding: a state factor is a signed Int32; j > 0 means u[j], j < 0 means conj(u[-j])
# (conj-folded systems reference folded-out partners as plain averages of the adjoint key).
# Monomial id 1 is the empty product, v[1] == 1.

# Typed error taxonomy (spec finding 4): AutoBackend catches exactly
# NonPolynomialDriftError; everything else stays a hard error.
abstract type KernelLoweringError <: Exception end
struct NonPolynomialDriftError{T} <: KernelLoweringError
    eqindex::Int
    residual::T
end
function Base.showerror(io::IO, e::NonPolynomialDriftError)
    return print(
        io,
        "NonPolynomialDriftError: equation $(e.eqindex) has a non-polynomial part: " *
            "$(e.residual). The moment-kernel path requires drifts polynomial in the " *
            "moments (anything meanfield/complete! produces); use the sharded or " *
            "ModelingToolkit path for rewritten non-polynomial drifts.",
    )
end
struct TimeDependentCoefficientError{T} <: KernelLoweringError
    coeff::T
end
function Base.showerror(io::IO, e::TimeDependentCoefficientError)
    return print(
        io,
        "TimeDependentCoefficientError: coefficient $(e.coeff) depends on the " *
            "independent variable. t-dependent coefficients are not supported by the " *
            "kernel path yet; use `ShardedBackend()` or the ModelingToolkit path.",
    )
end
struct ImParameterCollisionError <: KernelLoweringError end
function Base.showerror(io::IO, ::ImParameterCollisionError)
    return print(
        io,
        "ImParameterCollisionError: a user parameter named `im` collides with the " *
            "algebra's symbolic imaginary unit; rename the parameter.",
    )
end
struct HolomorphicJacobianError <: KernelLoweringError end
function Base.showerror(io::IO, ::HolomorphicJacobianError)
    return print(
        io,
        "HolomorphicJacobianError: the analytic Jacobian is holomorphic-only, and this " *
            "system's drift references conj(state) monomials (a conj-folded closure); its " *
            "true derivative needs the Wirtinger pair, which is not implemented. Close the " *
            "system with `get_adjoints = true` (the default) to unfold the conjugate " *
            "partners, or use an explicit solver without `jac = true`.",
    )
end
struct UnresolvedMomentError{T} <: KernelLoweringError
    moment::T
end
function Base.showerror(io::IO, e::UnresolvedMomentError)
    return print(
        io,
        "UnresolvedMomentError: the right-hand sides reference the average $(e.moment), " *
            "which does not resolve to any state. Call `complete(eqs)` first, or check " *
            "that the system was fully scaled/evaluated.",
    )
end

"""
Resolution of every drift average leaf through the system's recorded treatments, the
`_state_registry` pattern with the state INDEX as the `MomentMap` payload. Returns
`(vars, idx)` in the `_lower_ir` contract: `vars` are the distinct leaf forms as they
appear in the drifts (handed to `polynomial_coeffs`), `idx[leaf]` the signed state index
(negative = conjugate side of the stored representative).
"""
function statevars_resolved(eqs)
    g = eqs.graph
    ctx = build_ctx(eqs)
    treatments = _treatments(eqs, ctx)
    ops = QAdd[(o = undo_average(s); o isa QAdd ? o : o * 1) for s in eqs.states]
    moments = MomentMap(ctx, treatments, ops, collect(Int32, 1:length(ops)))
    idx = Dict{Any, Int32}()
    vars = Any[]
    for nd in values(g.nodes), leaf in eachleaf(Symbolics.unwrap(nd.drift))
        haskey(idx, leaf) && continue
        op = undo_average(leaf)
        r = match_moment(moments, op isa QAdd ? op : op * 1)
        r === nothing && throw(UnresolvedMomentError(leaf))
        i, same = r
        idx[leaf] = same ? i : Int32(-i)
        push!(vars, leaf)
    end
    return vars, idx
end

"""Signed factor list of one monomial, sorted canonically. Errors on unresolvable factors."""
function monomial_factors(mono, idx)
    fs = Int32[]
    addfac(f) = if SymbolicUtils.iscall(f) && SymbolicUtils.operation(f) === (^)
        b, e = SymbolicUtils.arguments(f)
        n = Int(SymbolicUtils.unwrap_const(e))
        j = idx[b]
        for _ in 1:n
            push!(fs, j)
        end
    else
        push!(fs, idx[f])
    end
    if mono isa Number || SymbolicUtils.isconst(mono)
        # empty product (constant term of the drift)
    elseif SymbolicUtils.iscall(mono) && SymbolicUtils.operation(mono) === (*)
        foreach(addfac, SymbolicUtils.arguments(mono))
    else
        addfac(mono)
    end
    return sort!(fs)
end

"""
Moment-polynomial IR. Monomial ids are prefix-closed: `parent[m]` is the id of the monomial
missing the last factor, `leaf[m]` that factor (signed state index). Parents are created
before children, so ids are already a valid update order. `coeffs` is the pooled list of
symbolic coefficient expressions; `coo` holds `(equation, monomial, coeff_id)` triples.
"""
struct MomentIR
    nstates::Int
    parent::Vector{Int32}
    leaf::Vector{Int32}
    coeffs::Vector{Any}
    coo_i::Vector{Int32}
    coo_j::Vector{Int32}
    coo_c::Vector{Int32}
    params::Vector{Any}
end

"""Lower a completed equation set to its moment-polynomial IR (treatments-aware)."""
lower(eqs) = _lower_ir(eqs.graph, statevars_resolved(eqs)..., Symbolics.unwrap(eqs.iv))

"""IR builder over a prepared state resolution (`vars` for `polynomial_coeffs`, `idx`
mapping each average leaf form to its signed state index)."""
function _lower_ir(g, vars, idx, iv_uw)
    mono_ids = Dict{Vector{Int32}, Int32}(Int32[] => Int32(1))
    parent = Int32[0]
    leaf = Int32[0]
    coeff_ids = Dict{Any, Int32}()
    coeffs = Any[]
    coo_i = Int32[]
    coo_j = Int32[]
    coo_c = Int32[]

    function mono_id!(fs::Vector{Int32})
        return get!(mono_ids, fs) do
            p = mono_id!(fs[1:(end - 1)])
            push!(parent, p)
            push!(leaf, fs[end])
            Int32(length(parent))
        end
    end

    # phase 1: polynomial_coeffs per equation, embarrassingly parallel (the dominant
    # lowering cost; measured 2.3x on 12 threads in the prototype)
    drifts = Any[Symbolics.unwrap(nd.drift) for nd in values(g.nodes)]
    neq = length(drifts)
    polys = Vector{Any}(undef, neq)
    Threads.@threads for i in 1:neq
        polys[i] = Symbolics.polynomial_coeffs(drifts[i], vars)
    end

    # phase 2: serial table build in equation order, so `mono_id!` assignment, coefficient
    # pooling, and the residual-check error order are bit-identical to a serial build
    for i in 1:neq
        dict, res = polys[i]
        SymbolicUtils._iszero(res) || throw(NonPolynomialDriftError(i, res))
        for (mono, c) in dict
            j = mono_id!(monomial_factors(mono, idx))
            cid = get!(coeff_ids, c) do
                push!(coeffs, c)
                Int32(length(coeffs))
            end
            push!(coo_i, Int32(i))
            push!(coo_j, j)
            push!(coo_c, cid)
        end
    end
    params = discover_params(coeffs, iv_uw)
    return MomentIR(length(g.nodes), parent, leaf, coeffs, coo_i, coo_j, coo_c, params)
end

"""
Union of the variables of each pooled coefficient. NOT `get_variables(sum(coeffs))`:
summing can cancel a parameter (coefficients `J` and `-J` sum to 0 and lose `J`).
Throws `TimeDependentCoefficientError` if the independent variable appears in a
coefficient, and `ImParameterCollisionError` for a user-created variable named `im`
(the algebra's symbolic imaginary unit has symtype Number; a user `@variables im` has
symtype Real and would silently be bound to `Base.im`).
"""
function discover_params(coeffs, iv = nothing)
    seen = Set{Any}()
    params = Any[]
    for c in coeffs
        (c isa Number || SymbolicUtils.isconst(c)) && continue
        for v in Symbolics.get_variables(c)
            u = Symbolics.unwrap(v)
            iv !== nothing && isequal(u, iv) && throw(TimeDependentCoefficientError(c))
            if SymbolicUtils.issym(u) && Base.nameof(u) === :im
                # the algebra's own imaginary unit (symtype Number) is not a parameter
                SymbolicUtils.symtype(u) === Number || throw(ImParameterCollisionError())
                continue
            end
            u in seen || (push!(seen, u); push!(params, u))
        end
    end
    return params
end

"""
Numeric coefficient values for a parameter assignment (zero codegen, `substitute`-based).
The algebra is SymReal-typed, so the imaginary unit is a symbolic `Sym{Number}(:im)` in the
coefficients (native codegen never sees this: `toexpr` emits the literal symbol `im`, which
resolves to `Base.im` in the generated code); the data path substitutes it explicitly.
"""
function coefficient_values(ir::MomentIR, pdict)
    for k in keys(pdict)
        u = Symbolics.unwrap(k)
        SymbolicUtils.issym(u) && Base.nameof(u) === :im &&
            throw(ImParameterCollisionError())
    end
    pd = Dict{Any, Any}(Symbolics.unwrap(k) => v for (k, v) in pdict)
    for c in ir.coeffs
        c isa Number && continue
        for v in Symbolics.get_variables(c)
            u = SymbolicUtils.unwrap(v)
            SymbolicUtils.issym(u) && Base.nameof(u) === :im && (pd[u] = im)
        end
    end
    return ComplexF64[
        ComplexF64(SymbolicUtils.unwrap_const(Symbolics.substitute(c, pd))) for c in ir.coeffs
    ]
end

"""Materialize the sparse coefficient matrix M (neq × nmonomials) for coefficient values `c`."""
assemble(ir::MomentIR, cvals::Vector{ComplexF64}) = sparse(
    ir.coo_i, ir.coo_j, cvals[ir.coo_c], ir.nstates, length(ir.parent), +,
)

# ---- array-aware parameter values ----------------------------------------------------

_pname(p) = SymbolicUtils.iscall(p) && SymbolicUtils.operation(p) === getindex ?
    Base.nameof(SymbolicUtils.arguments(p)[1]) : _param_name(p)
function _pslots(p)
    if SymbolicUtils.iscall(p) && SymbolicUtils.operation(p) === getindex
        return Int[
            a isa Number ? Int(a) : Int(SymbolicUtils.unwrap_const(a))
                for a in SymbolicUtils.arguments(p)[2:end]
        ]
    end
    return _param_slots(p)
end

"""
Substitution dict for `coefficient_values` from a `parameter_map(eqs, ...)` result:
scalar entries pass through keyed by their unwrapped sym; each discovered kernel
parameter that is an array access (`g[1]`, `Γ[2,1]`, or a callable indexed variable)
is matched by (name, concrete slots) against the array values. A `String` parameter
(cache-loaded kernels store printed names) is matched by printed name instead.

With `strict = false`, parameters that `pmap` does not determine are silently left out
instead of erroring (the partial-update path of `update_parameters!`, where missing
entries keep their stored values).
"""
function kernel_pdict(params::Vector, pmap; strict::Bool = true)
    pd = Dict{Any, Any}()
    arrs = Dict{Symbol, Any}()
    named = Dict{Any, Any}()
    byname = Dict{String, Any}()
    # two passes: `String`-keyed entries (stale values carried by a cache-loaded kernel)
    # first, so a symbolic entry for the same printed name always overwrites them
    for pass in (1, 2), (k, v) in pmap
        (pass == 1) == (k isa String) || continue
        if k isa String
            v isa AbstractArray ? (arrs[Symbol(k)] = v) : (byname[k] = v)
            continue
        end
        ku = Symbolics.unwrap(k)
        if v isa AbstractArray
            arrs[_pname(ku)] = v
        else
            pd[ku] = v
            byname[string(ku)] = v
            n = _pname(ku)
            n === nothing || (named[(n, _pslots(ku))] = v)
        end
    end
    unmatched = Any[]
    for p in params
        haskey(pd, p) && continue
        if p isa String
            m = match(r"^(.+?)\[([0-9,\s]+)\]$", p)
            if haskey(byname, p)
                pd[p] = byname[p]
            elseif m !== nothing && haskey(arrs, Symbol(something(m.captures[1])))
                slots = parse.(Int, split(something(m.captures[2]), ","))
                pd[p] = arrs[Symbol(something(m.captures[1]))][slots...]
            else
                push!(unmatched, p)
            end
            continue
        end
        name = _pname(p)
        slots = _pslots(p)
        if name !== nothing && haskey(arrs, name) && slots !== nothing
            pd[p] = arrs[name][slots...]
        elseif name !== nothing && haskey(named, (name, slots))
            pd[p] = named[(name, slots)]
        elseif name !== nothing && haskey(arrs, name)
            pd[p] = arrs[name]
        else
            push!(unmatched, p)
        end
    end
    strict && !isempty(unmatched) && throw(
        ArgumentError("missing values for kernel parameters: $(unmatched). Pass them in `ps`.")
    )
    return pd
end
