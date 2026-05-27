"""
    CorrelationFunction(op1, op2, eqs::MeanFieldEquations;
                        steady_state=true, filter_func=nothing)

Two-time correlation `⟨op1(τ) · op2(0)⟩` in a new independent variable τ.

Implementation uses the Quantum Regression Theorem via an *ancilla* subspace:
`op2` is reconstructed on a fresh subspace `aon_anc = max_aon + 1` so that
``[H, \\mathrm{op1\\,op2\\_anc}] = [H, \\mathrm{op1}] \\cdot \\mathrm{op2\\_anc}``,
and the τ-equation of `⟨op1 · op2_anc⟩` reproduces the regression dynamics.

When `steady_state=true`, the completion stage closes the system only on
averages that *touch the ancilla* — every other average appearing on the RHS
is interpreted as a constant steady-state coefficient.
"""
struct CorrelationFunction{T <: MeanFieldEquations, O1 <: QField, O2 <: QField, O2A <: QField}
    op1::O1
    op2::O2          # original op2 (acts on the original subspaces)
    op2_anc::O2A     # op2 lifted to the ancilla subspace
    aon_anc::Int
    eqs::T           # τ-equations for the ancilla-touching averages
    eqs0::MeanFieldEquations   # original-system equations (kept for u0/p0)
    τ::Symbolics.Num
    steady_state::Bool
end

# ---- ancilla embedding helpers ----------------------------------------------

# Recreate `op2` on subspace `aon`. Operators in SQA only carry an integer
# `space_index`, so embedding into an extended product space simply means
# constructing the same operator type with a different space index.
_embed_on(op::SQA.Destroy, aon::Int) = SQA.Destroy(op.name, aon, op.index)
_embed_on(op::SQA.Create, aon::Int) = SQA.Create(op.name, aon, op.index)
function _embed_on(op::SQA.Transition, aon::Int)
    return SQA.Transition(
        op.name, op.i, op.j, aon, op.index,
        op.ground_state, op.n_levels
    )
end

# Pick a fresh subspace index strictly greater than anything in `eqs0`.
function _ancilla_aon(eqs0::MeanFieldEquations, op1::QField, op2::QField)
    aons = Int[]
    append!(aons, SQA.acts_on(eqs0.hamiltonian))
    for j in eqs0.jumps
        append!(aons, SQA.acts_on(j))
    end
    for j in eqs0.jumps_dagger
        append!(aons, SQA.acts_on(j))
    end
    for op in eqs0.operators
        append!(aons, SQA.acts_on(op))
    end
    append!(aons, SQA.acts_on(op1))
    append!(aons, SQA.acts_on(op2))
    return isempty(aons) ? 2 : (maximum(aons) + 1)
end

# Filter accepting only LEAF averages whose operator product touches `aon_anc`
# AND at least one other subspace. Used as the `find_missing` predicate so we
# only derive equations for proper correlation states; pure ancilla averages
# (⟨op2_anc⟩ alone) are constants equal to ⟨op2⟩ at steady state, and
# non-ancilla averages are steady-state coefficients on the RHS.
function _ancilla_filter(aon_anc::Int)
    return function (x)
        SQA.is_average(x) || return true
        SymbolicUtils.iscall(x) || return true
        SymbolicUtils.operation(x) === SQA.sym_average || return true
        aons = SQA.acts_on(x)
        return (aon_anc in aons) && length(aons) > 1
    end
end

# Iterative completion that closes only on ancilla-touching averages and
# leaves the rest of the RHS untouched (those are steady-state coefficients).
function _complete_ancilla!(
        eqs::MeanFieldEquations, aon_anc::Int, steady_state::Bool,
        user_filter; max_iter::Int = 200, mix_choice = maximum,
        parent_canon = nothing,
    )
    anc_filter = _ancilla_filter(aon_anc)
    closure_filter = if user_filter === nothing
        anc_filter
    else
        x -> anc_filter(x) && user_filter(x)
    end
    for _ in 1:max_iter
        missing_states = find_missing(
            eqs; filter_func = closure_filter,
            get_adjoints = false, canon = parent_canon,
        )
        # When NOT in steady state we also need equations for the non-ancilla
        # averages (they evolve in τ too). Add those here.
        if !steady_state
            non_anc_missing = find_missing(
                eqs;
                filter_func = x -> !anc_filter(x),
                get_adjoints = false, canon = parent_canon,
            )
            append!(missing_states, non_anc_missing)
        end
        if isempty(missing_states)
            # If the user supplied a filter, zero out any leaf averages on the
            # RHS that the user rejected (and still touch the ancilla). Do not
            # zero ambient steady-state averages: those are valid coefficients.
            user_filter === nothing || _filter_anc_leaves!(eqs, anc_filter, user_filter)
            return eqs
        end
        new_ops = QField[_undo_for_derivation(m) for m in missing_states]
        # Same alpha-rename dance as `_derive_for(::MeanFieldEquations, ...)`:
        # if a LHS free index collides with an H/J bound name, the inner
        # meanfield commutator silently drops cross-atom terms. Rename to a
        # fresh successor, derive, then undo so the appended equations land
        # in the canonical (user-visible) names.
        bound = _bound_indices(eqs)
        fresh_ops, undo = _alpha_rename_away(new_ops, bound)
        new_eqs = _meanfield_deterministic(
            eqs.direction, fresh_ops,
            eqs.hamiltonian, eqs.jumps, eqs.jumps_dagger, eqs.rates,
            eqs.order, mix_choice, eqs.iv
        )
        if !isempty(undo)
            new_eqs = _apply_undo(new_eqs, undo)
        end
        _append!(eqs, new_eqs)
    end
    error("CorrelationFunction completion did not close within $max_iter iterations")
end

# Walk RHS and replace each LEAF average ⟨X⟩ that touches the ancilla but
# fails user_filter with 0. Ambient steady-state averages are untouched.
function _filter_anc_leaves!(eqs, anc_filter, user_filter)
    for (i, eq) in enumerate(eqs.equations)
        new_rhs = _filter_anc_expr(eq.rhs, anc_filter, user_filter)
        eqs.equations[i] = eq.lhs ~ new_rhs
    end
    return eqs
end

function _filter_anc_expr(x, anc_filter, user_filter)
    x isa SymbolicUtils.BasicSymbolic || return x
    if SQA.is_average(x) && SymbolicUtils.iscall(x) &&
            SymbolicUtils.operation(x) === SQA.sym_average
        # Leaf average. Only touch ancilla-touching leaves.
        return (anc_filter(x) && !user_filter(x)) ? 0 : x
    end
    SymbolicUtils.iscall(x) || return x
    op = SymbolicUtils.operation(x)
    args = SymbolicUtils.arguments(x)
    new_args = Any[_filter_anc_expr(a, anc_filter, user_filter) for a in args]
    op === complex && length(new_args) == 2 &&
        return new_args[1] + new_args[2] * Symbolics.IM
    return op(new_args...)
end

# ---- constructor ------------------------------------------------------------

function CorrelationFunction(
        op1::QField, op2::QField,
        eqs0::MeanFieldEquations;
        steady_state::Bool = true,
        filter_func = nothing,
    )
    τ = first(MTK.@independent_variables τ)
    aon_anc = _ancilla_aon(eqs0, op1, op2)
    op2_anc = _embed_on(op2, aon_anc)
    new_op = op1 * op2_anc

    # Extend the per-subspace order vector to cover the ancilla. The ancilla
    # is a single one-body operator, so its order doesn't matter; use the max
    # of the existing orders so it never restricts factorization.
    order_ext = if eqs0.order === nothing
        nothing
    else
        ord = eqs0.order
        ord_pad = vcat(ord, fill(maximum(ord), aon_anc - length(ord)))
        ord_pad
    end

    eqs_c = meanfield(
        [new_op], eqs0.hamiltonian, eqs0.jumps;
        Jdagger = eqs0.jumps_dagger,
        rates = eqs0.rates,
        order = order_ext,
        iv = τ,
    )

    # Iteratively derive equations only for averages that touch the ancilla,
    # leaving non-ancilla averages on the RHS as steady-state coefficients.
    # The user's filter_func (if given) further restricts which ancilla-touching
    # averages get their own equation. We thread the parent eqs0's canonical
    # index pool so the derived correlation states use the same atom-index
    # names as eqs0 (otherwise `_ambient_avgs` lookup against eqs0.states
    # fails and the spectrum collapses to a trivial constant).
    parent_canon = _build_canonical_indices(eqs0)
    _complete_ancilla!(eqs_c, aon_anc, steady_state, filter_func; parent_canon)

    return CorrelationFunction(op1, op2, op2_anc, aon_anc, eqs_c, eqs0, τ, steady_state)
end

# ---- spectrum ---------------------------------------------------------------

"""
    Spectrum(c::CorrelationFunction, ps)

The Laplace-transformed correlation function as a callable. Calling
`S(ω, u_end, p0)` returns the (real, normalized) spectrum at frequency `ω`,
given the steady-state mapping `u_end` (`avg => value`) and the numeric
parameter values `p0` for `ps`.

Solves the (Laplace-domain) linear system `(iω·I - A) X̃(ω) = x̃(0)`, where
`x(τ)` are the ancilla-touching averages, `A` is the Jacobian of the τ-equation
RHS w.r.t. those states, `x̃(τ) = x(τ) - x_∞` is the centred process (so the
correlation decays to zero, matching the Wiener-Khinchin assumption), and
`x_∞ = -A⁻¹ b_const` is the asymptote driven by the constant source from
ambient steady-state averages. Returns `2·Re(X̃(ω)[1])`, the symmetric (two-sided)
power spectrum `S(ω) = 2 Re{∫₀^∞ g(τ) e^{-iωτ} dτ}` consistent with master and
with `QuantumOptics.timecorrelations.correlation2spectrum`.
"""
struct Spectrum{C <: CorrelationFunction, P}
    c::C
    ω::Symbolics.Num
    ps::P
end

function Spectrum(c::CorrelationFunction, ps)
    ω = first(@variables ω_spectrum)
    return Spectrum{typeof(c), typeof(ps)}(c, ω, ps)
end

Spectrum(c::CorrelationFunction) = Spectrum(c, ())

function (S::Spectrum)(ω_vals::AbstractVector, u_end, p0)
    A, rhs_b, n = _spectrum_kernel(S, u_end, p0)
    I_n = Matrix{ComplexF64}(I, n, n)
    out = Vector{Float64}(undef, length(ω_vals))
    @inbounds for (k, ω_val) in enumerate(ω_vals)
        M = (im * ω_val) .* I_n .- A
        x = M \ rhs_b
        out[k] = 2 * real(x[1])
    end
    return out
end

(S::Spectrum)(ω_val::Real, u_end, p0) = first(S([ω_val], u_end, p0))

# Build (A, rhs_b, n) once per call by symbolic finite differences on the
# τ-equations (linear in the states, so this is exact). All ω-dependence
# lives in the per-ω solve below. We deliberately avoid
# `Symbolics.build_function`: its JIT-compile of expressions containing
# averaged-operator atoms takes seconds at higher orders, far more than
# substitution-based extraction for a one-shot Spectrum invocation.
# `rhs_b = u_τ - x_∞` centres x̃ on the asymptote x_∞ = -A⁻¹ b_const so the
# one-sided FT converges. For phase-invariant systems b_const = 0.
function _spectrum_kernel(S::Spectrum, u_end, p0)
    c = S.c
    eqs = c.eqs
    n = length(eqs.equations)

    rhss_u = [SymbolicUtils.unwrap(eq.rhs) for eq in eqs.equations]
    p_sub = _build_p_sub(S.ps, p0, c.τ, u_end, eqs)
    u_end_dict = _as_avg_dict(c, u_end)
    for avg in _ambient_avgs(c)
        aons = SQA.acts_on(avg)
        lookup_avg = (c.aon_anc in aons && length(aons) == 1) ?
            _undo_ancilla(c, avg) : avg
        p_sub[avg] = ComplexF64(_lookup_avg(u_end_dict, lookup_avg))
    end

    state_us = [SymbolicUtils.unwrap(dv) for dv in eqs.states]
    zero_sub = merge(p_sub, Dict{Any, Any}(s => 0.0 + 0.0im for s in state_us))

    b_const = Vector{ComplexF64}(undef, n)
    @inbounds for (i, rhs) in enumerate(rhss_u)
        b_const[i] = _scalarize(SymbolicUtils.substitute(rhs, zero_sub))
    end
    A = Matrix{ComplexF64}(undef, n, n)
    @inbounds for j in 1:n
        sub_j = merge(zero_sub, Dict{Any, Any}(state_us[j] => 1.0 + 0.0im))
        for (i, rhs) in enumerate(rhss_u)
            f_ej = _scalarize(SymbolicUtils.substitute(rhs, sub_j))
            A[i, j] = f_ej - b_const[i]
        end
    end

    u_τ = zeros(ComplexF64, n)
    @inbounds for (i, s) in enumerate(eqs.states)
        avg0 = _undo_ancilla(c, SymbolicUtils.unwrap(s))
        u_τ[i] = ComplexF64(_lookup_avg(u_end_dict, avg0))
    end

    rhs_b = any(!iszero, b_const) ? u_τ .+ (A \ b_const) : u_τ
    return A, rhs_b, n
end

function _build_p_sub(ps, p0, τ, u_end, eqs)
    sub = Dict{Any, Any}()
    if ps !== nothing && !isempty(ps)
        for (p, v) in zip(ps, p0)
            sub[SymbolicUtils.unwrap(p)] = ComplexF64(v)
        end
    end
    if u_end !== nothing
        if u_end isa AbstractDict
            for (k, v) in u_end
                ku = SymbolicUtils.unwrap(k)
                sub[ku] = ComplexF64(v)
                # also map ⟨k'⟩ ↦ conj(v) so RHS coefficients with the conjugate
                # average are resolved without requiring the user to spell them out.
                if SQA.is_average(ku)
                    sub[SymbolicUtils.unwrap(_avg_conj_of(ku))] = ComplexF64(conj(v))
                end
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

function _avg_conj_of(x)
    SQA.is_average(x) || return x
    op = SQA.undo_average(x)
    return average(adjoint(op))
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
    # The substituted RHS lives in SymReal land where the imaginary unit is
    # the symbolic constant `Symbolics.IM`. Replace it with the numeric `im`
    # and let Symbolics fold the constant expression.
    y = SymbolicUtils.substitute(x, Dict(SymbolicUtils.unwrap(Symbolics.IM) => im))
    v = Symbolics.value(y)
    v isa Number && return ComplexF64(v)
    return ComplexF64(0)
end
_scalarize(x::Number) = ComplexF64(x)
_scalarize(x) = ComplexF64(0)

# ---- ODE / MTK bridge -------------------------------------------------------

# Collect "ambient" averages that appear in the RHS of `c.eqs` but are NOT
# τ-evolving. Two cases:
#  - averages not touching the ancilla at all (pure original-system
#    coefficients);
#  - averages touching ONLY the ancilla (⟨op2_anc⟩-type), which are constants
#    equal to the original-system ⟨op2⟩ at steady state.
# Both are mapped to MTK parameters in `System(c)` and given numeric values
# in `correlation_p0`.
function _ambient_avgs(c::CorrelationFunction)
    found = OrderedCollections.OrderedSet{SymbolicUtils.BasicSymbolic}()
    states_set = Set(SymbolicUtils.unwrap.(c.eqs.states))
    for eq in c.eqs.equations
        _collect_ambient!(found, eq.rhs, c.aon_anc, states_set)
    end
    return collect(found)
end

function _collect_ambient!(found, x, aon_anc, states_set)
    x isa SymbolicUtils.BasicSymbolic || return
    if SQA.is_average(x) && SymbolicUtils.operation(x) === SQA.sym_average
        x in states_set && return
        aons = SQA.acts_on(x)
        if !(aon_anc in aons) || length(aons) == 1
            push!(found, x)
        end
        return
    end
    SymbolicUtils.iscall(x) || return
    for a in SymbolicUtils.arguments(x)
        _collect_ambient!(found, a, aon_anc, states_set)
    end
    return
end

function _ambient_param(avg::SymbolicUtils.BasicSymbolic)
    # Prefix with `ss_` so ambient steady-state parameters never collide with
    # state variable names. Different subspaces aren't reflected in
    # `_stable_avg_name`, so `⟨a · σ⟩` (original) and `⟨a_anc · σ⟩` (ancilla
    # state) would otherwise hash to the same Symbol.
    name = Symbol("ss_", _stable_avg_name(avg))
    v = first(@variables $name)
    return v
end

"""
    correlation_u0(c, u_end)

Build u0 for the τ-evolution. Each ancilla-touching state `⟨X(τ)·op2_anc⟩`
gets its initial value from the steady-state `⟨X·op2⟩` (where `op2` is on
the *original* subspace), looked up in `u_end`. Missing averages default to
0. `u_end` may be either an `avg => value` `Dict` or a `Vector` aligned with
`c.eqs0.states`.
"""
function correlation_u0(c::CorrelationFunction, u_end)
    u_end_dict = _as_avg_dict(c, u_end)
    u0 = Dict{Symbolics.Num, ComplexF64}()
    dict, _ = _avg_to_var_dict(c.eqs)
    for s in c.eqs.states
        avg0 = _undo_ancilla(c, s)
        v = _lookup_avg(u_end_dict, avg0)
        u0[dict[s]] = ComplexF64(v)
    end
    return u0
end

# Replace op2_anc inside the operator product of `avg` with op2 (on the
# original subspace), producing the τ=0 representative. Rebuilds the product
# through SQA arithmetic so the canonical pipeline fires and the resulting
# average matches whatever `eqs0` would have stored.
function _undo_ancilla(c::CorrelationFunction, avg::SymbolicUtils.BasicSymbolic)
    SQA.is_average(avg) || return avg
    op = SQA.undo_average(avg)
    op isa SQA.QAdd || return avg
    aon0 = c.op2.space_index
    aon_anc = c.aon_anc
    result = zero(op)
    for (term, coeff) in op.arguments
        prod = one(SQA.QAdd) * coeff
        for o in term.ops
            o2 = (o.space_index == aon_anc) ? _embed_on(o, aon0) : o
            prod = prod * o2
        end
        result = result + prod
    end
    return average(result)
end

function _as_avg_dict(c::CorrelationFunction, u_end)
    if u_end isa AbstractDict
        return Dict{SymbolicUtils.BasicSymbolic, Any}(
            SymbolicUtils.unwrap(k) => v for (k, v) in u_end
        )
    elseif u_end isa AbstractVector
        return Dict{SymbolicUtils.BasicSymbolic, Any}(
            SymbolicUtils.unwrap(avg) => v
                for (avg, v) in zip(c.eqs0.states, u_end)
        )
    else
        throw(ArgumentError("u_end must be Dict or Vector, got $(typeof(u_end))"))
    end
end

function _lookup_avg(dict, avg)
    avg_u = SymbolicUtils.unwrap(avg)
    # Fast path: leaf average present as-is (or via its conjugate) in `dict`.
    haskey(dict, avg_u) && return dict[avg_u]
    avg_u isa Number && return avg_u
    if avg_u isa SymbolicUtils.BasicSymbolic
        if _is_leaf_average(avg_u)
            conj_u = SymbolicUtils.unwrap(_avg_conj_of(avg_u))
            haskey(dict, conj_u) && return conj(dict[conj_u])
            return 0
        end
        # Compound expression: `_undo_ancilla` rebuilds via SQA arithmetic,
        # which normal-orders products (e.g. `a·a' → a'·a + 1`). Substitute
        # every known leaf average (and its conjugate) with its numeric value
        # and evaluate. Unknown leaves default to 0 via `_scalarize`.
        sub = Dict{Any, Any}()
        for (k, v) in dict
            sub[k] = v
            k isa SymbolicUtils.BasicSymbolic && _is_leaf_average(k) || continue
            ck = SymbolicUtils.unwrap(_avg_conj_of(k))
            haskey(sub, ck) || (sub[ck] = conj(v))
        end
        return _scalarize(SymbolicUtils.substitute(avg_u, sub))
    end
    return 0
end

"""
    correlation_p0(c, u_end, ps_p0)

Build the parameter dict for the τ-evolution ODE. Combines the user-provided
`ps_p0` (parameter symbol => value pairs) with the numeric values of the
ambient steady-state averages on `c.eqs` RHS (looked up in `u_end`).
"""
function correlation_p0(c::CorrelationFunction, u_end, ps_p0)
    u_end_dict = _as_avg_dict(c, u_end)
    out = Dict{Any, ComplexF64}()
    for (k, v) in ps_p0
        out[k] = ComplexF64(v)
    end
    for avg in _ambient_avgs(c)
        p = _ambient_param(avg)
        aons = SQA.acts_on(avg)
        lookup_avg = (c.aon_anc in aons && length(aons) == 1) ?
            _undo_ancilla(c, avg) : avg
        out[p] = ComplexF64(_lookup_avg(u_end_dict, lookup_avg))
    end
    return out
end

"""
    ModelingToolkitBase.System(c::CorrelationFunction; name::Symbol)

Build an MTK `System` for the τ-evolution of the correlation function.
Ambient steady-state averages on the RHS are substituted with MTK parameters
whose values are supplied via [`correlation_p0`](@ref).
"""
function MTK.System(c::CorrelationFunction; name::Symbol)
    eqs = c.eqs
    iv = eqs.iv
    iv_uw = SymbolicUtils.unwrap(iv)
    D = Symbolics.Differential(iv)
    dict, dvs = _avg_to_var_dict(eqs)
    conj_dict = _conj_substitution_dict(eqs, dict)

    # Ambient steady-state averages → real parameters.
    ambient = _ambient_avgs(c)
    ambient_subs = Dict{Any, Any}()
    ambient_params = Symbolics.Num[]
    for avg in ambient
        p = _ambient_param(avg)
        ambient_subs[avg] = p
        # Conjugate average maps to conj(p) (only used if it sneaks into the RHS).
        conj_avg = SymbolicUtils.unwrap(_avg_conj_of(avg))
        if conj_avg !== avg
            ambient_subs[conj_avg] = conj(p)
        end
        push!(ambient_params, p)
    end

    new_eqs = Vector{Symbolics.Equation}(undef, length(eqs.equations))
    ps_set = Set{SymbolicUtils.BasicSymbolic}()
    @inbounds for (i, eq) in enumerate(eqs.equations)
        # Custom leaf substitution: SymbolicUtils.substitute reconstructs
        # interior nodes via their operation, which drops the AvgSym symtype
        # on avg nodes that are *not* themselves substituted (they become
        # symtype=Any and stop being recognised by is_average / substitute).
        # So we collapse conj, ambient, and dv substitution into one tree walk.
        merged = Dict{Any, Any}()
        for (k, v) in conj_dict
            merged[k] = v
        end
        for (k, v) in ambient_subs
            merged[k] = v
        end
        for (k, v) in dict
            merged[k] = v
        end
        rhs = _safe_substitute(eq.rhs, merged)
        new_eqs[i] = D(dict[eq.lhs]) ~ rhs
        _collect_params!(ps_set, rhs, dict, iv_uw)
    end
    ps_user = [MTK.toparam(p) for p in collect(ps_set)]
    ps_ambient = [MTK.toparam(SymbolicUtils.unwrap(p)) for p in ambient_params]
    ps = vcat(ps_user, ps_ambient)
    return MTK.System(new_eqs, iv, dvs, ps; name = name)
end

"""
    scale(c::CorrelationFunction)

Scale the ancilla τ-equations of a `CorrelationFunction`. The underlying
`eqs0` (original-system equations) is also scaled so the spectrum and
steady-state lookup at `Spectrum(c, ps)(ω, u_end, p0)` resolves against the
same canonical state names as the τ system.
"""
function scale(c::CorrelationFunction)
    return CorrelationFunction(
        c.op1, c.op2, c.op2_anc, c.aon_anc,
        scale(c.eqs), scale(c.eqs0), c.τ, c.steady_state,
    )
end
