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
        SymbolicUtils._iszero(res) || error(
            "equation $i has a non-polynomial part: $res. The moment-kernel path requires " *
                "drifts polynomial in the moments (anything meanfield/complete! produces); " *
                "use the ModelingToolkit path for rewritten non-polynomial drifts.",
        )
        for (mono, c) in dict
            j = mono_id!(monomial_factors(mono, idx))
            cid = get!(coeff_ids, c) do
                push!(coeffs, c)
                Int32(length(coeffs))
            end
            push!(coo_i, Int32(i)); push!(coo_j, j); push!(coo_c, cid)
        end
    end
    params = Symbolics.unwrap.(Symbolics.get_variables(sum(coeffs)))
    return MomentIR(length(g.nodes), parent, leaf, coeffs, coo_i, coo_j, coo_c, params)
end

"""
Numeric coefficient values for a parameter assignment (zero codegen, `substitute`-based).
The algebra is SymReal-typed, so the imaginary unit is a symbolic `Sym{Number}(:im)` in the
coefficients (native codegen never sees this: `toexpr` emits the literal symbol `im`, which
resolves to `Base.im` in the generated code); the data path substitutes it explicitly.
"""
function coefficient_values(ir::MomentIR, pdict)
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
