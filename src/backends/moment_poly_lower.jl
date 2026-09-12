# Recursive sparse-polynomial compiler for MomentIR.
#
# The supported grammar is deliberately explicit. State-free subtrees are coefficients;
# state-dependent sums, products, nonnegative integer powers, unary/binary subtraction, and
# division by a state-free denominator are polynomial structure. Any other state-dependent
# operation declines lowering and is reported by the MomentIR builder. There is intentionally
# no Symbolics.polynomial_coeffs fallback in production.

_unwrap_moment_poly(x::Symbolics.Num) = SymbolicUtils.unwrap(x)
_unwrap_moment_poly(x) = x

function _moment_state_factor(x, idx)
    x = _unwrap_moment_poly(x)
    haskey(idx, x) && return idx[x]
    if x isa SymbolicUtils.BasicSymbolic && SymbolicUtils.iscall(x)
        op = SymbolicUtils.operation(x)
        args = SymbolicUtils.arguments(x)
        if op === conj && length(args) == 1
            y = _unwrap_moment_poly(args[1])
            haskey(idx, y) && return Int32(-idx[y])
        end
    end
    return nothing
end

function _contains_moment_state(x, idx, cache::IdDict{Any, Bool})
    x = _unwrap_moment_poly(x)
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
    e = _unwrap_moment_poly(e)
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

"""Conservative zero test local to cold coefficient construction."""
function _moment_coeff_iszero(coeff)
    u = _unwrap_moment_poly(coeff)
    u isa Number && return iszero(u)
    if u isa SymbolicUtils.BasicSymbolic && SymbolicUtils.isconst(u)
        return try
            iszero(SymbolicUtils.unwrap_const(u))
        catch
            false
        end
    end
    if u isa SymbolicUtils.BasicSymbolic &&
            SymbolicUtils.iscall(u) &&
            SymbolicUtils.operation(u) === (+)
        expanded = try
            _unwrap_moment_poly(Symbolics.expand(coeff))
        catch
            u
        end
        expanded isa Number && return iszero(expanded)
        if expanded isa SymbolicUtils.BasicSymbolic && SymbolicUtils.isconst(expanded)
            return try
                iszero(SymbolicUtils.unwrap_const(expanded))
            catch
                false
            end
        end
    end
    return false
end

function _merge_moment_factors(a::Tuple, b::Tuple)
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

function _moment_poly_add_term!(out::Dict{Tuple, Any}, mono::Tuple, coeff)
    _moment_coeff_iszero(coeff) && return out
    if haskey(out, mono)
        value = out[mono] + coeff
        if _moment_coeff_iszero(value)
            delete!(out, mono)
        else
            out[mono] = value
        end
    else
        out[mono] = coeff
    end
    return out
end

function _moment_poly_add!(out::Dict{Tuple, Any}, p::Dict{Tuple, Any})
    for (mono, coeff) in p
        _moment_poly_add_term!(out, mono, coeff)
    end
    return out
end

function _moment_poly_scale(p::Dict{Tuple, Any}, coeff)
    _moment_coeff_iszero(coeff) && return Dict{Tuple, Any}()
    isequal(coeff, 1) && return p
    out = Dict{Tuple, Any}()
    sizehint!(out, length(p))
    for (mono, value) in p
        _moment_poly_add_term!(out, mono, coeff * value)
    end
    return out
end

function _constant_moment_poly_coeff(p::Dict{Tuple, Any})
    length(p) == 1 || return nothing
    return get(p, (), nothing)
end

function _moment_poly_mul(a::Dict{Tuple, Any}, b::Dict{Tuple, Any})
    isempty(a) && return Dict{Tuple, Any}()
    isempty(b) && return Dict{Tuple, Any}()
    ca = _constant_moment_poly_coeff(a)
    ca === nothing || return _moment_poly_scale(b, ca)
    cb = _constant_moment_poly_coeff(b)
    cb === nothing || return _moment_poly_scale(a, cb)

    out = Dict{Tuple, Any}()
    sizehint!(out, length(a) * length(b))
    for (ma, caa) in a, (mb, cbb) in b
        _moment_poly_add_term!(out, _merge_moment_factors(ma, mb), caa * cbb)
    end
    return out
end

function _moment_poly_pow(base::Dict{Tuple, Any}, n::Int)
    n == 0 && return Dict{Tuple, Any}(() => 1)
    n == 1 && return base
    result = Dict{Tuple, Any}(() => 1)
    power = base
    k = n
    while k > 0
        isodd(k) && (result = _moment_poly_mul(result, power))
        k >>= 1
        k == 0 || (power = _moment_poly_mul(power, power))
    end
    return result
end

"""
    _compile_moment_polynomial(x, idx, state_cache)

Compile `x` to `Dict{Tuple,Any}` mapping canonical signed state-factor tuples to symbolic
coefficients. Returns `nothing` precisely when a state-dependent subtree is outside the
supported polynomial grammar.
"""
function _compile_moment_polynomial(x, idx, state_cache::IdDict{Any, Bool})
    x = _unwrap_moment_poly(x)

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
            _moment_poly_add!(out, p)
        end
        return out
    elseif op === (*)
        out = Dict{Tuple, Any}(() => 1)
        for arg in args
            p = _compile_moment_polynomial(arg, idx, state_cache)
            p === nothing && return nothing
            out = _moment_poly_mul(out, p)
        end
        return out
    elseif op === (-)
        if length(args) == 1
            p = _compile_moment_polynomial(args[1], idx, state_cache)
            p === nothing && return nothing
            return _moment_poly_scale(p, -1)
        elseif length(args) == 2
            a = _compile_moment_polynomial(args[1], idx, state_cache)
            a === nothing && return nothing
            b = _compile_moment_polynomial(args[2], idx, state_cache)
            b === nothing && return nothing
            return _moment_poly_add!(a, _moment_poly_scale(b, -1))
        end
        return nothing
    elseif op === (^) && length(args) == 2
        _contains_moment_state(args[2], idx, state_cache) && return nothing
        n = _nonnegative_integer_exponent(args[2])
        n === nothing && return nothing
        base = _compile_moment_polynomial(args[1], idx, state_cache)
        base === nothing && return nothing
        return _moment_poly_pow(base, n)
    elseif op === (/) && length(args) == 2
        _contains_moment_state(args[2], idx, state_cache) && return nothing
        numerator = _compile_moment_polynomial(args[1], idx, state_cache)
        numerator === nothing && return nothing
        return _moment_poly_scale(numerator, inv(_unwrap_moment_poly(args[2])))
    end

    return nothing
end
