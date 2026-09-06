# Fast polynomial lowering for the direct moment kernel.
#
# Completed cumulant drifts are polynomials in moment states. Preserve that structure while
# lowering instead of expanding every equation through Symbolics.polynomial_coeffs. The
# generic Symbolics path remains an equation-local compatibility fallback for rewritten
# expressions outside the small polynomial grammar handled here.

_unwrap_poly(x::Symbolics.Num) = SymbolicUtils.unwrap(x)
_unwrap_poly(x) = x

function _moment_state_factor(x, idx)
    x = _unwrap_poly(x)
    haskey(idx, x) && return idx[x]
    if x isa SymbolicUtils.BasicSymbolic && SymbolicUtils.iscall(x)
        op = SymbolicUtils.operation(x)
        args = SymbolicUtils.arguments(x)
        if op === conj && length(args) == 1
            y = _unwrap_poly(args[1])
            haskey(idx, y) && return Int32(-idx[y])
        end
    end
    return nothing
end

function _contains_moment_state(x, idx, cache::IdDict{Any, Bool})
    x = _unwrap_poly(x)
    haskey(cache, x) && return cache[x]
    result = if _moment_state_factor(x, idx) !== nothing
        true
    elseif !(x isa SymbolicUtils.BasicSymbolic) || !SymbolicUtils.iscall(x)
        false
    else
        any(a -> _contains_moment_state(a, idx, cache), SymbolicUtils.arguments(x))
    end
    cache[x] = result
    return result
end

function _nonnegative_integer_exponent(e)
    e = _unwrap_poly(e)
    value = if e isa Number
        e
    else
        try
            SymbolicUtils.unwrap_const(e)
        catch
            return nothing
        end
    end
    n = if value isa Integer
        Int(value)
    elseif value isa Rational && denominator(value) == 1
        Int(value)
    elseif value isa Real && isinteger(value)
        Int(value)
    else
        return nothing
    end
    return n < 0 ? nothing : n
end

_kernel_coeff_iszero(c) = try
    _iszero_part(c)
catch
    false
end

function _merge_factor_tuples(a::Tuple, b::Tuple)
    isempty(a) && return b
    isempty(b) && return a
    out = Vector{Int32}(undef, length(a) + length(b))
    ia = 1
    ib = 1
    io = 1
    while ia <= length(a) && ib <= length(b)
        if a[ia] <= b[ib]
            out[io] = a[ia]
            ia += 1
        else
            out[io] = b[ib]
            ib += 1
        end
        io += 1
    end
    while ia <= length(a)
        out[io] = a[ia]
        ia += 1
        io += 1
    end
    while ib <= length(b)
        out[io] = b[ib]
        ib += 1
        io += 1
    end
    return Tuple(out)
end

function _poly_add_term!(out::Dict{Tuple, Any}, mono::Tuple, coeff)
    _kernel_coeff_iszero(coeff) && return out
    if haskey(out, mono)
        value = out[mono] + coeff
        if _kernel_coeff_iszero(value)
            delete!(out, mono)
        else
            out[mono] = value
        end
    else
        out[mono] = coeff
    end
    return out
end

function _poly_add!(out::Dict{Tuple, Any}, p::Dict{Tuple, Any})
    for (mono, coeff) in p
        _poly_add_term!(out, mono, coeff)
    end
    return out
end

function _poly_scale(p::Dict{Tuple, Any}, coeff)
    _kernel_coeff_iszero(coeff) && return Dict{Tuple, Any}()
    isequal(coeff, 1) && return p
    out = Dict{Tuple, Any}()
    sizehint!(out, length(p))
    for (mono, value) in p
        _poly_add_term!(out, mono, coeff * value)
    end
    return out
end

function _constant_poly_coeff(p::Dict{Tuple, Any})
    length(p) == 1 || return nothing
    return get(p, (), nothing)
end

function _poly_mul(a::Dict{Tuple, Any}, b::Dict{Tuple, Any})
    isempty(a) && return Dict{Tuple, Any}()
    isempty(b) && return Dict{Tuple, Any}()
    ca = _constant_poly_coeff(a)
    ca === nothing || return _poly_scale(b, ca)
    cb = _constant_poly_coeff(b)
    cb === nothing || return _poly_scale(a, cb)

    out = Dict{Tuple, Any}()
    sizehint!(out, length(a) * length(b))
    for (ma, caa) in a, (mb, cbb) in b
        _poly_add_term!(out, _merge_factor_tuples(ma, mb), caa * cbb)
    end
    return out
end

function _poly_pow(base::Dict{Tuple, Any}, n::Int)
    n == 0 && return Dict{Tuple, Any}(() => 1)
    n == 1 && return base
    result = Dict{Tuple, Any}(() => 1)
    power = base
    k = n
    while k > 0
        isodd(k) && (result = _poly_mul(result, power))
        k >>= 1
        k == 0 || (power = _poly_mul(power, power))
    end
    return result
end

"""
    _compile_moment_polynomial(x, idx, state_cache) -> Union{Dict, Nothing}

Recursively lower a symbolic expression to a sparse polynomial whose keys are canonical
signed state-factor tuples. Subtrees without moment states are preserved intact as
coefficients. Returns `nothing` when a state-dependent operation lies outside the supported
polynomial grammar; the caller then uses the generic Symbolics fallback for that equation.
"""
function _compile_moment_polynomial(x, idx, state_cache::IdDict{Any, Bool})
    x = _unwrap_poly(x)

    j = _moment_state_factor(x, idx)
    j === nothing || return Dict{Tuple, Any}((Int32(j),) => 1)

    if !_contains_moment_state(x, idx, state_cache)
        return Dict{Tuple, Any}(() => x)
    end

    if !(x isa SymbolicUtils.BasicSymbolic) || !SymbolicUtils.iscall(x)
        return nothing
    end

    op = SymbolicUtils.operation(x)
    args = SymbolicUtils.arguments(x)

    if op === (+)
        out = Dict{Tuple, Any}()
        for arg in args
            p = _compile_moment_polynomial(arg, idx, state_cache)
            p === nothing && return nothing
            _poly_add!(out, p)
        end
        return out
    elseif op === (*)
        out = Dict{Tuple, Any}(() => 1)
        for arg in args
            p = _compile_moment_polynomial(arg, idx, state_cache)
            p === nothing && return nothing
            out = _poly_mul(out, p)
        end
        return out
    elseif op === (^) && length(args) == 2
        _contains_moment_state(args[2], idx, state_cache) && return nothing
        n = _nonnegative_integer_exponent(args[2])
        n === nothing && return nothing
        base = _compile_moment_polynomial(args[1], idx, state_cache)
        base === nothing && return nothing
        return _poly_pow(base, n)
    elseif op === (/) && length(args) == 2
        _contains_moment_state(args[2], idx, state_cache) && return nothing
        numerator = _compile_moment_polynomial(args[1], idx, state_cache)
        numerator === nothing && return nothing
        return _poly_scale(numerator, inv(_unwrap_poly(args[2])))
    elseif op === (-) && length(args) == 1
        p = _compile_moment_polynomial(args[1], idx, state_cache)
        p === nothing && return nothing
        return _poly_scale(p, -1)
    end

    return nothing
end

function _generic_moment_terms(drift, vars, idx, eqindex)
    dict, residual = Symbolics.polynomial_coeffs(drift, vars)
    _iszero_part(residual) || throw(NonPolynomialDriftError(eqindex, residual))
    terms = Tuple{Tuple, Any}[]
    sizehint!(terms, length(dict))
    for (mono, coeff) in dict
        push!(terms, (Tuple(monomial_factors(mono, idx)), coeff))
    end
    sort!(terms; by = first)
    return terms
end

function _moment_equation_terms(drift, vars, idx, eqindex, state_cache)
    poly = _compile_moment_polynomial(drift, idx, state_cache)
    poly === nothing && return _generic_moment_terms(drift, vars, idx, eqindex)
    terms = Tuple{Tuple, Any}[(mono, coeff) for (mono, coeff) in poly]
    sort!(terms; by = first)
    return terms
end

# More-specific production builder. The generic builder in kernel_lower.jl remains available
# as the compatibility/reference implementation, while normal statevars_resolved output
# dispatches here without the expensive whole-equation polynomial conversion.
function _build_moment_ir(
        g::MomentGraph,
        vars::Vector{Any},
        idx::Dict{Any, Int32},
        iv_uw,
    )
    edges = Dict{Tuple{Int32, Int32}, Int32}()
    parent = Int32[0]
    leaf = Int32[0]
    coeff_ids = Dict{Any, Int32}()
    coeffs = Any[]
    coo_i = Int32[]
    coo_j = Int32[]
    coo_c = Int32[]

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

    drifts = Any[Symbolics.unwrap(nd.drift) for nd in values(g.nodes)]
    state_cache = IdDict{Any, Bool}()

    # Equation order and monomial order are deterministic. Unsupported rewritten expressions
    # fall back only for that equation, preserving the existing error taxonomy.
    for i in eachindex(drifts)
        terms = _moment_equation_terms(drifts[i], vars, idx, i, state_cache)
        for (factors, coeff) in terms
            _kernel_coeff_iszero(coeff) && continue
            j = mono_id!(factors)
            cid = get!(coeff_ids, coeff) do
                push!(coeffs, coeff)
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
