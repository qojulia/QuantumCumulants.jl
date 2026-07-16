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
# (conj-folded systems reference folded-out partners as plain averages of the adjoint key,
# see the Kerr probe in the design doc). Monomial id 1 is the empty product, v[1] == 1.

using QuantumCumulants
using QuantumCumulants: canon_key
using Symbolics
using SymbolicUtils
using SymbolicUtils: iscall, operation, arguments
using SparseArrays

# Typed error taxonomy (spec finding 4): AutoBackend catches exactly
# NonPolynomialDriftError; everything else stays a hard error.
abstract type KernelLoweringError <: Exception end
struct NonPolynomialDriftError <: KernelLoweringError
    eqindex::Int
    residual::Any
end
function Base.showerror(io::IO, e::NonPolynomialDriftError)
    print(
        io,
        "NonPolynomialDriftError: equation $(e.eqindex) has a non-polynomial part: " *
            "$(e.residual). The moment-kernel path requires drifts polynomial in the " *
            "moments (anything meanfield/complete! produces); use the sharded or " *
            "ModelingToolkit path for rewritten non-polynomial drifts.",
    )
end
struct TimeDependentCoefficientError <: KernelLoweringError
    coeff::Any
end
function Base.showerror(io::IO, e::TimeDependentCoefficientError)
    print(
        io,
        "TimeDependentCoefficientError: coefficient $(e.coeff) depends on the " *
            "independent variable. t-dependent coefficients are not supported by the " *
            "kernel path yet; use the ModelingToolkit path.",
    )
end
struct ImParameterCollisionError <: KernelLoweringError end
function Base.showerror(io::IO, ::ImParameterCollisionError)
    print(
        io,
        "ImParameterCollisionError: a user parameter named `im` collides with the " *
            "algebra's symbolic imaginary unit; rename the parameter.",
    )
end
struct UnresolvedMomentError <: KernelLoweringError
    moment::Any
end
function Base.showerror(io::IO, e::UnresolvedMomentError)
    print(
        io,
        "UnresolvedMomentError: the right-hand sides reference the average $(e.moment), " *
            "which does not resolve to any state. Call `complete(eqs)` first, or check " *
            "that the system was fully scaled/evaluated.",
    )
end

"""
States plus folded conjugate partners of a completed graph, as unwrapped average terms.
Returns `(vars, idx)` where `vars` is the list to hand `polynomial_coeffs` and
`idx[var] = signed state index` (negative = conjugate of the folded-out partner).
State order is the `keys(g.nodes)` order, the same order `initial_values` uses.
"""
function statevars(g)
    ks = collect(keys(g.nodes))
    idx = Dict{Any, Int32}()
    vars = Any[]
    for (i, k) in enumerate(ks)
        v = Symbolics.unwrap(average(k))
        idx[v] = Int32(i)
        push!(vars, v)
        kc = canon_key(adjoint(k), g.ctx)
        if !isequal(kc, k) && !haskey(g.nodes, kc)
            vc = Symbolics.unwrap(average(kc))
            idx[vc] = Int32(-i)
            push!(vars, vc)
        end
    end
    return vars, idx
end

"""Signed factor list of one monomial, sorted canonically. Errors on unresolvable factors."""
function monomial_factors(mono, idx)
    fs = Int32[]
    addfac(f) = if iscall(f) && operation(f) === (^)
        b, e = arguments(f)
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
    elseif iscall(mono) && operation(mono) === (*)
        foreach(addfac, arguments(mono))
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

function lower(eqs)
    g = eqs.graph
    vars, idx = statevars(g)
    return _lower_ir(g, vars, idx, Symbolics.unwrap(eqs.iv))
end

"""IR builder over a prepared state resolution (`vars` for `polynomial_coeffs`, `idx`
mapping each average leaf form to its signed state index)."""
function _lower_ir(g, vars, idx, iv_uw)
    mono_ids = Dict{Vector{Int32}, Int32}(Int32[] => Int32(1))
    parent = Int32[0]
    leaf = Int32[0]
    coeff_ids = Dict{Any, Int32}()
    coeffs = Any[]
    coo_i = Int32[]; coo_j = Int32[]; coo_c = Int32[]

    function mono_id!(fs::Vector{Int32})
        get!(mono_ids, fs) do
            p = mono_id!(fs[1:(end - 1)])
            push!(parent, p)
            push!(leaf, fs[end])
            Int32(length(parent))
        end
    end

    for (i, nd) in enumerate(values(g.nodes))
        dict, res = Symbolics.polynomial_coeffs(Symbolics.unwrap(nd.drift), vars)
        SymbolicUtils._iszero(res) || throw(NonPolynomialDriftError(i, res))
        for (mono, c) in dict
            j = mono_id!(monomial_factors(mono, idx))
            cid = get!(coeff_ids, c) do
                push!(coeffs, c)
                Int32(length(coeffs))
            end
            push!(coo_i, Int32(i)); push!(coo_j, j); push!(coo_c, cid)
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
            if SymbolicUtils.issym(u) && SymbolicUtils.nameof(u) === :im
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
        SymbolicUtils.issym(u) && SymbolicUtils.nameof(u) === :im &&
            throw(ImParameterCollisionError())
    end
    pd = Dict{Any, Any}(Symbolics.unwrap(k) => v for (k, v) in pdict)
    for c in ir.coeffs
        c isa Number && continue
        for v in Symbolics.get_variables(c)
            u = SymbolicUtils.unwrap(v)
            SymbolicUtils.issym(u) && SymbolicUtils.nameof(u) === :im && (pd[u] = im)
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
