# Structured lowering of completed deterministic cumulant equations.
#
# A completed MeanfieldEquations is represented as
#
#     duᵢ = Σₖ cₖ(p) mₖ(u)
#
# without committing to a numerical execution strategy. Shared monomials form a
# prefix-closed DAG. Equation rows are stored directly as CSR-like term ranges; this is
# the mathematical incidence structure, not a SparseArrays storage trick.

abstract type MomentLoweringError <: Exception end

struct NonPolynomialDriftError{T} <: MomentLoweringError
    eqindex::Int
    expression::T
end
function Base.showerror(io::IO, e::NonPolynomialDriftError)
    return print(
        io,
        "NonPolynomialDriftError: equation $(e.eqindex) contains unsupported " *
            "state-dependent structure: $(e.expression). The direct moment representation " *
            "accepts polynomial drift expressions only; use `System(eqs)` for a general " *
            "symbolic system.",
    )
end

struct TimeDependentCoefficientError{T} <: MomentLoweringError
    coeff::T
end
function Base.showerror(io::IO, e::TimeDependentCoefficientError)
    return print(
        io,
        "TimeDependentCoefficientError: coefficient $(e.coeff) depends on the independent " *
            "variable. Time-dependent coefficients are outside the direct moment " *
            "representation; use `System(eqs)` for that model.",
    )
end

struct UnresolvedMomentError{T} <: MomentLoweringError
    moment::T
end
function Base.showerror(io::IO, e::UnresolvedMomentError)
    return print(
        io,
        "UnresolvedMomentError: the right-hand sides reference the average $(e.moment), " *
            "which does not resolve to any stored state. Call `complete(eqs)` first, or " *
            "check that scaling/evaluation produced a closed system.",
    )
end

"""
    MomentIR

Evaluator-independent representation of a completed deterministic moment hierarchy.

`states[i]` identifies row/state `i`. Monomial `1` is the empty product. For `m > 1`,
`parent[m]` and `leaf[m]` encode `monomial[m] = monomial[parent[m]] * statefactor(leaf[m])`.
A positive leaf `j` denotes state `j`; a negative leaf `-j` denotes `conj(state[j])`.

Equation `i` owns the term range `rowptr[i]:(rowptr[i + 1] - 1)`. Each term references one
shared monomial and one pooled symbolic coefficient. Symbolic metadata is intentionally cold:
numerical evaluators must compile this representation into concrete runtime storage rather
than inspect `states`, `coeffs`, or `params` on the hot RHS path.
"""
struct MomentIR
    states::Vector{Any}
    parent::Vector{Int32}
    leaf::Vector{Int32}
    rowptr::Vector{Int32}
    monomial::Vector{Int32}
    coeff_id::Vector{Int32}
    coeffs::Vector{Any}
    params::Vector{Any}
end

Base.length(ir::MomentIR) = length(ir.states)

"""Map every average leaf occurring in the equations to its signed stored-state index."""
function _moment_state_indices(eqs::MeanfieldEquations)
    ctx = build_ctx(eqs)
    treatments = _treatments(eqs, ctx)
    ops = QAdd[(o = undo_average(s); o isa QAdd ? o : o * 1) for s in eqs.states]
    moments = MomentMap(ctx, treatments, ops, collect(Int32, 1:length(ops)))

    idx = Dict{Any, Int32}()
    for eq in eqs.equations, leaf in eachleaf(Symbolics.unwrap(eq.rhs))
        haskey(idx, leaf) && continue
        op = undo_average(leaf)
        r = match_moment(moments, op isa QAdd ? op : op * 1)
        r === nothing && throw(UnresolvedMomentError(leaf))
        i, same = r
        idx[leaf] = same ? Int32(i) : Int32(-i)
    end
    return idx
end

"""Discover scalar symbolic coefficient dependencies in deterministic order."""
function _moment_params(coeffs, iv)
    iv_uw = Symbolics.unwrap(iv)
    seen = Set{Any}()
    params = Any[]
    for coeff in coeffs
        coeff isa Number && continue
        for var in Symbolics.get_variables(coeff)
            u = Symbolics.unwrap(var)
            isequal(u, iv_uw) && throw(TimeDependentCoefficientError(coeff))
            # SQA represents the algebraic imaginary unit symbolically. It is a constant,
            # not a model parameter. A user variable also named `im` remains a parameter
            # because its symbolic type is not Number.
            if SymbolicUtils.issym(u) &&
                    Base.nameof(u) === :im &&
                    SymbolicUtils.symtype(u) === Number
                continue
            end
            if !(u in seen)
                push!(seen, u)
                push!(params, u)
            end
        end
    end
    sort!(params; by = string)
    return params
end

"""
    _lower_moment_ir(eqs::MeanfieldEquations) -> MomentIR

Lower a completed deterministic hierarchy to its structured polynomial representation.
The output is independent of parameter values and numerical evaluator choice.
"""
function _lower_moment_ir(eqs::MeanfieldEquations)
    idx = _moment_state_indices(eqs)
    return _build_moment_ir(eqs, idx)
end

function _build_moment_ir(eqs::MeanfieldEquations, idx::Dict{Any, Int32})
    edges = Dict{Tuple{Int32, Int32}, Int32}()
    parent = Int32[0]
    leaf = Int32[0]

    coeff_ids = Dict{Any, Int32}()
    coeffs = Any[]
    rowptr = Vector{Int32}(undef, length(eqs.equations) + 1)
    rowptr[1] = Int32(1)
    monomial = Int32[]
    coeff_id = Int32[]

    function mono_id!(factors::Tuple)
        p = Int32(1)
        @inbounds for factor in factors
            f = Int32(factor)
            key = (p, f)
            p = get!(edges, key) do
                push!(parent, p)
                push!(leaf, f)
                Int32(length(parent))
            end
        end
        return p
    end

    state_cache = IdDict{Any, Bool}()
    for (i, eq) in enumerate(eqs.equations)
        drift = Symbolics.unwrap(eq.rhs)
        poly = _compile_moment_polynomial(drift, idx, state_cache)
        poly === nothing && throw(NonPolynomialDriftError(i, drift))

        terms = collect(poly)
        sort!(terms; by = first)
        for (factors, coeff) in terms
            _moment_coeff_iszero(coeff) && continue
            mid = mono_id!(factors)
            cid = get!(coeff_ids, coeff) do
                push!(coeffs, coeff)
                Int32(length(coeffs))
            end
            push!(monomial, mid)
            push!(coeff_id, cid)
        end
        rowptr[i + 1] = Int32(length(monomial) + 1)
    end

    states = Any[Symbolics.unwrap(s) for s in eqs.states]
    params = _moment_params(coeffs, eqs.iv)
    return MomentIR(states, parent, leaf, rowptr, monomial, coeff_id, coeffs, params)
end
