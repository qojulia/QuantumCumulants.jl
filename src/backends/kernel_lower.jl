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

# Typed errors make unsupported direct execution explicit. The direct path never falls back
# to generated code: callers can choose `System(eqs)` and the ModelingToolkit path instead.
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
            "moments (anything meanfield/complete! produces). Use " *
            "`System(eqs)` and the ModelingToolkit path for rewritten non-polynomial drifts.",
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
            "direct path; use `System(eqs)` and the ModelingToolkit path.",
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
            "system with `get_adjoints = true` to unfold the conjugate partners, or use an " *
            "explicit solver without `jac = true`.",
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
`(vars, idx)` in the lowering contract: `vars` are the distinct leaf forms as they appear
in the drifts, and `idx[leaf]` is the signed state index (negative = conjugate side of the
stored representative).
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
    addfac(f) =
    if SymbolicUtils.iscall(f) && SymbolicUtils.operation(f) === (^)
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

"""Lower a completed equation set to its moment-polynomial representation."""
_lower_moment_ir(eqs) =
    _build_moment_ir(eqs.graph, statevars_resolved(eqs)..., Symbolics.unwrap(eqs.iv))

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
        SymbolicUtils.issym(u) &&
            Base.nameof(u) === :im &&
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
    return ComplexF64[_numeric_coefficient(c, pd) for c in ir.coeffs]
end

function _numeric_coefficient(c, pd)
    substituted = Symbolics.substitute(c, pd)
    value = SymbolicUtils.unwrap_const(substituted)
    try
        return ComplexF64(value)
    catch err
        # Array-valued parameter substitution can leave a real/imag wrapper around a concrete
        # complex scalar. Simplify only this cold construction fallback, never on the RHS.
        simplified = Symbolics.simplify(substituted; expand = true)
        value = SymbolicUtils.unwrap_const(simplified)
        try
            return ComplexF64(value)
        catch
            throw(err)
        end
    end
end

"""Materialize the sparsity pattern of `Mᵀ` (nmonomials × neq).

The numeric coefficient values live in `KernelParameters`, which is the `prob.p` payload;
the `MomentKernel` retains this pattern only. Storing the transpose (M in CSR) makes the RHS
a row-parallel gather; see `MomentKernel`.
"""
assemble_pattern(ir::MomentIR) = sparse(
    ir.coo_j,
    ir.coo_i,
    fill(Int32(1), length(ir.coo_i)),
    length(ir.parent),
    ir.nstates,
    +,
)

# ---- array-aware parameter values ----------------------------------------------------

_pname(p) =
    SymbolicUtils.iscall(p) && SymbolicUtils.operation(p) === getindex ?
    Base.nameof(SymbolicUtils.arguments(p)[1]) : _param_name(p)
function _pslots(p)
    if SymbolicUtils.iscall(p) && SymbolicUtils.operation(p) === getindex
        return Int[
            a isa Number ? Int(a) : Int(SymbolicUtils.unwrap_const(a)) for
                a in SymbolicUtils.arguments(p)[2:end]
        ]
    end
    return _param_slots(p)
end

"""
Substitution dict for `coefficient_values` from a `parameter_map(eqs, ...)` result:
scalar entries pass through keyed by their unwrapped symbolic identity; each discovered
kernel parameter that is an array access (`g[1]`, `Γ[2,1]`, or a callable indexed variable)
is matched by its structural name and concrete slots against the array value.

With `strict = false`, parameters that `pmap` does not determine are silently left out
instead of erroring (the partial-update path of `update_parameters!`, where missing
entries keep their stored values).
"""
function kernel_pdict(params::Vector, pmap; strict::Bool = true)
    pd = Dict{Any, Any}()
    arrs = Dict{Symbol, Any}()
    named = Dict{Any, Any}()
    for (k, v) in pmap
        ku = Symbolics.unwrap(k)
        name = _pname(ku)
        if v isa AbstractArray
            name === nothing ? (pd[ku] = v) : (arrs[name] = v)
        else
            pd[ku] = v
            name === nothing || (named[(name, _pslots(ku))] = v)
        end
    end
    unmatched = Any[]
    for p in params
        haskey(pd, p) && continue
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
    strict &&
        !isempty(unmatched) &&
        throw(
        ArgumentError(
            "missing values for kernel parameters: $(unmatched). Pass them in `ps`.",
        ),
    )
    return pd
end
