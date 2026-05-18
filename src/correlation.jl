"""
    CorrelationFunction(op1, op2, eqs::MeanFieldEquations;
                        steady_state=true, simplify=true)

Two-time correlation `⟨op1(0) · op2(τ)⟩` in a new independent variable τ.
Reuses the underlying meanfield/cumulant/complete pipeline on a new
equation set with `iv = τ` and the original parameters held fixed.
"""
struct CorrelationFunction{T<:MeanFieldEquations, O1<:QField, O2<:QField}
    op1::O1
    op2::O2
    eqs::T
    τ::Symbolics.Num
    steady_state::Bool
end

function CorrelationFunction(op1::QField, op2::QField,
                             eqs0::MeanFieldEquations;
                             steady_state::Bool = true,
                             simplify::Bool = true)
    τ = first(MTK.@independent_variables τ)
    new_op = op1 * op2
    eqs_c = meanfield([new_op], eqs0.hamiltonian, eqs0.jumps;
                      Jdagger = eqs0.jumps_dagger,
                      rates   = eqs0.rates,
                      order   = eqs0.order,
                      simplify, iv = τ)
    complete!(eqs_c)
    return CorrelationFunction(op1, op2, eqs_c, τ, steady_state)
end

"""
    Spectrum(c::CorrelationFunction, ps; simplify=true)

The Laplace-transformed correlation function as a callable. Calling
`S(ω, u_end, p0)` returns the (real, normalized) spectrum at frequency `ω`,
given the steady-state vector `u_end` of the original system and the numeric
parameter values `p0` for `ps`.

The current implementation solves the (Laplace-domain) linear system
`A(ω) x(ω) = b(ω)` derived from the τ-evolution equations of the correlation
function. `S(ω) = real(x_1(ω))` is returned and normalized at the end.
"""
struct Spectrum{T<:MeanFieldEquations,P}
    c::CorrelationFunction{T,<:QField,<:QField}
    ω::Symbolics.Num
    ps::P
end

function Spectrum(c::CorrelationFunction, ps; simplify::Bool = true)
    ω = first(@variables ω_spectrum)
    return Spectrum{typeof(c.eqs), typeof(ps)}(c, ω, ps)
end

Spectrum(c::CorrelationFunction; simplify::Bool = true) = Spectrum(c, (); simplify)

function (S::Spectrum)(ω_vals::AbstractVector, u_end, p0)
    return [S(ω_val, u_end, p0) for ω_val in ω_vals]
end

function (S::Spectrum)(ω_val::Real, u_end, p0)
    c = S.c
    eqs = c.eqs
    n = length(eqs.equations)
    A = zeros(ComplexF64, n, n)
    b = zeros(ComplexF64, n)
    rhss = [eq.rhs for eq in eqs.equations]
    p_sub = _build_p_sub(S.ps, p0, c.τ, u_end, eqs)
    @inbounds for (i, rhs) in enumerate(rhss)
        # Linear decomposition: rhs = sum_j A[i,j] * state_j + b[i]
        for (j, dv) in enumerate(eqs.states)
            coeff = Symbolics.derivative(rhs, dv)
            coeff_n = SymbolicUtils.substitute(coeff, p_sub)
            A[i, j] = _scalarize(coeff_n)
        end
        # Constant term: substitute states with 0 to get the affine part
        zero_sub = merge(p_sub,
            Dict{Any, Any}(SymbolicUtils.unwrap(dv) => 0.0 for dv in eqs.states))
        const_n = SymbolicUtils.substitute(rhs, zero_sub)
        b[i] = _scalarize(const_n)
    end
    M = (im * ω_val) .* Matrix{ComplexF64}(I, n, n) .- A
    x = M \ b
    return real(x[1])
end

function _build_p_sub(ps, p0, τ, u_end, eqs)
    sub = Dict{Any, Any}()
    if ps !== nothing && !isempty(ps)
        for (p, v) in zip(ps, p0)
            sub[SymbolicUtils.unwrap(p)] = ComplexF64(v)
        end
    end
    if u_end !== nothing
        # Map steady-state ⟨op⟩ to numeric values
        if u_end isa AbstractDict
            for (k, v) in u_end
                sub[SymbolicUtils.unwrap(k)] = ComplexF64(v)
            end
        elseif u_end isa AbstractVector
            for (avg, v) in zip(eqs.states, u_end)
                sub[avg] = ComplexF64(v)
            end
        end
    end
    sub[SymbolicUtils.unwrap(τ)] = 0.0
    return sub
end

function _scalarize(x::Symbolics.Num)
    return _scalarize(SymbolicUtils.unwrap(x))
end
function _scalarize(x::SymbolicUtils.BasicSymbolic)
    if SymbolicUtils.isconst(x)
        v = x.val
        v isa Number && return ComplexF64(v)
        return ComplexF64(0)
    end
    return ComplexF64(0)
end
_scalarize(x::Number) = ComplexF64(x)
_scalarize(x) = ComplexF64(0)

# Coefficient of `var` in linear `expr` (ignoring nonlinear terms).
function _coeff_of(expr, var)
    expr_u = SymbolicUtils.unwrap(expr)
    var_u = SymbolicUtils.unwrap(var)
    derivative = Symbolics.derivative(expr_u, var_u)
    return derivative
end

"""
    correlation_u0(c, u_end)

Build u0 for the τ-evolution from the steady-state of the original system.
Currently returns the supplied `u_end` mapping verbatim, padded with zeros
for averages introduced by the τ-system.
"""
function correlation_u0(c::CorrelationFunction, u_end::AbstractDict)
    u0 = Dict{Symbolics.Num, ComplexF64}()
    dict, _ = _avg_to_var_dict(c.eqs)
    for (avg, v) in dict
        u0[v] = ComplexF64(get(u_end, avg, 0))
    end
    return u0
end

"""
    correlation_p0(c, p, u_end)

Pass-through of `p` for now. Frozen-state injection of `op1`-touching
averages can be added on top of this base.
"""
function correlation_p0(c::CorrelationFunction, p::AbstractDict,
                        u_end::AbstractDict)
    return Dict{Symbolics.Num, ComplexF64}(k => ComplexF64(v) for (k, v) in p)
end
