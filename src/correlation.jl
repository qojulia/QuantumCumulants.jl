"""
    CorrelationFunction

Two-time correlation function ⟨op1(t+τ) op2(t)⟩ of a solved mean-field system,
constructed via the quantum regression theorem. Build one with the
[`CorrelationFunction(op1, op2, eqs0)`](@ref) constructor, then feed it to
`ModelingToolkitBase.System` together with [`correlation_u0`](@ref) and
[`correlation_p0`](@ref) to solve the τ-evolution, or to [`Spectrum`](@ref) for the
power spectrum.

# Fields
- `op1`, `op2`: the operators whose correlation is computed.
- `op2_ancilla`: `op2` re-embedded on the ancilla subspace `aon_ancilla` (internal).
- `aon_ancilla`: Hilbert-subspace index of the ancilla carrying `op2` (internal).
- `eqs`: the closed equations of motion in the delay `τ`.
- `eqs0`: the steady-state system the correlation is built on.
- `τ`: the delay independent variable.
- `steady_state`: if `true`, ambient averages are held constant; if `false`, they are
  evolved alongside the τ-states.
"""
struct CorrelationFunction{T <: MeanfieldEquations, O1 <: QField, O2 <: QField, O2A <: QField}
    op1::O1
    op2::O2
    op2_ancilla::O2A
    aon_ancilla::Int
    eqs::T
    eqs0::MeanfieldEquations
    τ::Symbolics.Num
    steady_state::Bool
end

# ---- ancilla embedding -------------------------------------------------------

_embed_on(op::SQA.Destroy, aon::Int) = SQA.Destroy(op.name, aon, op.index)
_embed_on(op::SQA.Create, aon::Int) = SQA.Create(op.name, aon, op.index)
_embed_on(op::SQA.Transition, aon::Int) =
    SQA.Transition(op.name, op.i, op.j, aon, op.index, op.ground_state, op.n_levels)

function _ancilla_aon(eqs0::MeanfieldEquations, op1::QField, op2::QField)
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
    append!(aons, SQA.acts_on(op1)); append!(aons, SQA.acts_on(op2))
    return isempty(aons) ? 2 : (maximum(aons) + 1)
end

"""
Closure filter accepting only leaf averages whose product touches `aon_ancilla` and at
least one other subspace: the proper correlation states. Pure-ancilla averages
(⟨op2_ancilla⟩) and non-ancilla averages are steady-state coefficients, not τ-states.
"""
function _ancilla_filter(aon_ancilla::Int)
    return function (x)
        SQA.is_average(x) || return true
        SymbolicUtils.iscall(x) || return true
        SymbolicUtils.operation(x) === SQA.sym_average || return true
        aons = SQA.acts_on(x)
        return (aon_ancilla in aons) && length(aons) > 1
    end
end

# ---- constructor -------------------------------------------------------------

"""
    CorrelationFunction(op1, op2, eqs0; steady_state=true, filter_func=nothing, max_iter=100_000)

Build the two-time correlation ⟨op1(τ)·op2(0)⟩ from a solved mean-field system
`eqs0` via the quantum regression theorem. `op2` is lifted onto a fresh ancilla
subspace and the resulting system is closed only on averages that touch the
ancilla; the remaining averages are steady-state coefficients. With
`steady_state=false` the ambient averages are evolved too. `filter_func` further
restricts the closure, and `max_iter` bounds the closure iterations.

# Examples
```jldoctest
julia> h = FockSpace(:cavity);

julia> @qnumbers a::Destroy(h);

julia> eqs = meanfield([a' * a], a' * a, [a]; rates = [1.0], order = 2);

julia> CorrelationFunction(a', a, eqs)
⟨a' * a⟩
```
"""
function CorrelationFunction(
        op1::QField, op2::QField, eqs0::MeanfieldEquations;
        steady_state::Bool = true, filter_func = nothing, max_iter::Int = 100_000,
    )
    τ = first(MTK.@independent_variables τ)
    aon_ancilla = _ancilla_aon(eqs0, op1, op2)
    op2_ancilla = _embed_on(op2, aon_ancilla)
    new_op = op1 * op2_ancilla

    order_ext = eqs0.order === nothing ? nothing :
        vcat(eqs0.order, fill(maximum(eqs0.order), aon_ancilla - length(eqs0.order)))

    eqs_c = meanfield(
        [new_op], eqs0.hamiltonian, eqs0.jumps;
        Jdagger = eqs0.jumps_dagger, rates = eqs0.rates, order = order_ext, iv = τ,
    )

    ancilla_filter = _ancilla_filter(aon_ancilla)
    closure_filter = filter_func === nothing ? ancilla_filter : x -> ancilla_filter(x) && filter_func(x)
    g = _graph_from_eqs(eqs_c)
    closure!(g; filter = closure_filter, get_adjoints = false, max_iter)
    if !steady_state
        closure!(g; filter = x -> !ancilla_filter(x), get_adjoints = false, max_iter)
    end
    # Inherit eqs0's per-subspace treatments so the resolver matches the
    # correlation's leaves in the same treatment as the parent system.
    if !isempty(eqs0.treatments)
        merged = merge(g.treatments, eqs0.treatments)
        g = MomentGraph(g.nodes, g.sys, g.ctx, merged)
    end
    eqs_c = assemble_equations(g)
    filter_func === nothing || _filter_ancilla_leaves!(eqs_c, ancilla_filter, filter_func)

    return CorrelationFunction(op1, op2, op2_ancilla, aon_ancilla, eqs_c, eqs0, τ, steady_state)
end

"""
Zero the ancilla-touching leaves the user's filter rejects; ambient leaves stay.
"""
function _filter_ancilla_leaves!(eqs, ancilla_filter, user_filter)
    for (i, eq) in enumerate(eqs.equations)
        eqs.equations[i] = eq.lhs ~ _filter_ancilla_expr(eq.rhs, ancilla_filter, user_filter)
    end
    return eqs
end
function _filter_ancilla_expr(x, ancilla_filter, user_filter)
    x isa SymbolicUtils.BasicSymbolic || return x
    if _is_avg_leaf(x)
        return (ancilla_filter(x) && !user_filter(x)) ? 0 : x
    end
    SymbolicUtils.iscall(x) || return x
    op = SymbolicUtils.operation(x)
    new_args = Any[_filter_ancilla_expr(a, ancilla_filter, user_filter) for a in SymbolicUtils.arguments(x)]
    op === complex && length(new_args) == 2 && return new_args[1] + new_args[2] * Symbolics.IM
    return op(new_args...)
end

# ---- ambient / lookup helpers ------------------------------------------------

"""
Averages appearing on the RHS that are not τ-states: non-ancilla coefficients or
pure-ancilla constants. These become MTK parameters in `System(c)`.
"""
function _ambient_avgs(c::CorrelationFunction)
    found = OrderedCollections.OrderedSet{SymbolicUtils.BasicSymbolic}()
    states_set = Set(SymbolicUtils.unwrap.(c.eqs.states))
    for eq in c.eqs.equations
        _collect_ambient!(found, eq.rhs, c.aon_ancilla, states_set)
    end
    return collect(found)
end
function _collect_ambient!(found, x, aon_ancilla, states_set)
    x isa SymbolicUtils.BasicSymbolic || return
    if SQA.is_average(x) && SymbolicUtils.operation(x) === SQA.sym_average
        x in states_set && return
        aons = SQA.acts_on(x)
        (!(aon_ancilla in aons) || length(aons) == 1) && push!(found, x)
        return
    end
    SymbolicUtils.iscall(x) || return
    for a in SymbolicUtils.arguments(x)
        _collect_ambient!(found, a, aon_ancilla, states_set)
    end
    return
end

_ambient_param(avg::SymbolicUtils.BasicSymbolic) =
    first(@variables $(Symbol("ss_", serialize(undo_average(avg)))))

_avg_conj_of(x) = SQA.is_average(x) ? average(adjoint(undo_average(x))) : x

"""
Replace `op2_ancilla` with `op2` (the original subspace) inside `avg`, recovering the
τ=0 representative.
"""
function _undo_ancilla(c::CorrelationFunction, avg::SymbolicUtils.BasicSymbolic)
    SQA.is_average(avg) || return avg
    op = undo_average(avg)
    op isa QAdd || return avg
    aon0 = c.op2.space_index; aon_ancilla = c.aon_ancilla
    result = zero(op)
    for (term, coeff) in op.arguments
        prod = one(QAdd) * coeff
        for o in term.ops
            prod = prod * ((o.space_index == aon_ancilla) ? _embed_on(o, aon0) : o)
        end
        result = result + prod
    end
    return average(result)
end

function _as_avg_dict(c::CorrelationFunction, u_end)
    if u_end isa AbstractDict
        return Dict{SymbolicUtils.BasicSymbolic, Any}(SymbolicUtils.unwrap(k) => v for (k, v) in u_end)
    elseif u_end isa AbstractVector
        return Dict{SymbolicUtils.BasicSymbolic, Any}(
            SymbolicUtils.unwrap(avg) => v for (avg, v) in zip(c.eqs0.states, u_end)
        )
    end
    throw(ArgumentError("u_end must be Dict or Vector, got $(typeof(u_end))"))
end

"""
Steady-state value resolver for the correlation `c`, built once from `u_end`. Returns a
closure mapping an average (a leaf, or a normal-ordered expression of averages) to its
numeric steady-state value, `0` for an absent leaf. Matching is by the parent system's
Hermitian-conjugate representative (`canonical_rep` in `c.eqs0`'s context), not by raw
symbol identity, so a correlation average resolves to its steady-state value even when the
correlation and the parent mint different base index names for the same representative atom
(⟨σ_iₓ₂₂⟩ vs ⟨σ_jₓ₂₂⟩, as when the Hamiltonian's sum index differs from the user's LHS
index). The conjugate of ⟨X⟩ resolves to ⟨X⟩* via the representative's side bit.
"""
function _ss_resolver(c::CorrelationFunction, u_end)
    u_end_dict = _as_avg_dict(c, u_end)
    ctx = build_ctx(c.eqs0)
    treatments = _treatments(c.eqs0, ctx)
    ks = [k for k in keys(u_end_dict) if undo_average(k) isa QAdd]
    m = MomentMap(
        ctx, treatments,
        QAdd[undo_average(k) for k in ks],
        ComplexF64[ComplexF64(u_end_dict[k]) for k in ks],
    )
    resolve_leaf = function (leaf)
        r = match_moment(m, undo_average(leaf))
        r === nothing && return nothing
        v, same = r
        return same ? v : conj(v)
    end
    return function (avg)
        avg_u = SymbolicUtils.unwrap(avg)
        avg_u isa Number && return ComplexF64(avg_u)
        avg_u isa SymbolicUtils.BasicSymbolic || return ComplexF64(0)
        if _is_avg_leaf(avg_u)
            r = resolve_leaf(avg_u)
            return r === nothing ? ComplexF64(0) : r
        end
        sub = Dict{Any, Any}()
        for l in eachleaf(avg_u)
            r = resolve_leaf(l)
            sub[l] = r === nothing ? 0.0 + 0.0im : r
        end
        return _scalarize(SymbolicUtils.substitute(avg_u, sub))
    end
end

_scalarize(x::Symbolics.Num) = _scalarize(SymbolicUtils.unwrap(x))
function _scalarize(x::SymbolicUtils.BasicSymbolic)
    SymbolicUtils.isconst(x) && (v = x.val; return v isa Number ? ComplexF64(v) : ComplexF64(0))
    y = SymbolicUtils.substitute(x, Dict(SymbolicUtils.unwrap(Symbolics.IM) => im))
    v = Symbolics.value(y)
    return v isa Number ? ComplexF64(v) : ComplexF64(0)
end
_scalarize(x::Number) = ComplexF64(x)
_scalarize(x) = ComplexF64(0)

# ---- u0 / p0 -----------------------------------------------------------------

"""
    correlation_u0(c, u_end)

Initial values for the τ-evolution of the [`CorrelationFunction`](@ref) `c`. Each
ancilla state ⟨X(τ) op2_ancilla⟩ starts from the steady-state value ⟨X op2⟩ (with `op2` on
its original subspace), looked up in `u_end`.

See also: [`CorrelationFunction`](@ref), [`correlation_p0`](@ref).
"""
function correlation_u0(c::CorrelationFunction, u_end)
    resolve = _ss_resolver(c, u_end)
    reg = _state_registry(c.eqs)
    u0 = Dict{Symbolics.Num, ComplexF64}()
    for (i, s) in enumerate(c.eqs.states)
        u0[reg.vars[i]] = resolve(_undo_ancilla(c, SymbolicUtils.unwrap(s)))
    end
    return u0
end

"""
    correlation_p0(c, u_end, ps_p0)

Parameter dictionary for the τ-ODE of the [`CorrelationFunction`](@ref) `c`: the user
parameters `ps_p0` together with the numeric values of the ambient steady-state averages,
looked up in `u_end`.

See also: [`CorrelationFunction`](@ref), [`correlation_u0`](@ref).
"""
function correlation_p0(c::CorrelationFunction, u_end, ps_p0)
    resolve = _ss_resolver(c, u_end)
    out = Dict{Any, ComplexF64}()
    for (k, v) in ps_p0
        out[k] = ComplexF64(v)
    end
    for avg in _ambient_avgs(c)
        aons = SQA.acts_on(avg)
        lookup = (c.aon_ancilla in aons && length(aons) == 1) ? _undo_ancilla(c, avg) : avg
        out[_ambient_param(avg)] = resolve(lookup)
    end
    return out
end

# ---- System(c) ---------------------------------------------------------------

"""
    ModelingToolkitBase.System(c::CorrelationFunction; name)

MTK `System` for the τ-evolution. Ambient steady-state averages become parameters
(values supplied via [`correlation_p0`](@ref)).
"""
function MTK.System(c::CorrelationFunction; name::Symbol)
    eqs = c.eqs
    iv = eqs.iv; iv_uw = SymbolicUtils.unwrap(iv); D = Symbolics.Differential(iv)
    reg = _state_registry(eqs)
    # Ambient steady-state averages, keyed by the same Hermitian-conjugate representative as
    # the state registry (so a folded conjugate ambient resolves via the side bit).
    amb_avgs = _ambient_avgs(c)
    amb = MomentMap(
        reg.ctx, reg.treatments,
        [undo_average(avg) for avg in amb_avgs],
        Symbolics.Num[_ambient_param(avg) for avg in amb_avgs],
    )
    # State variable first, then ambient parameter; both via the same conjugation-side rule.
    resolve = function (leaf)
        op = undo_average(leaf)
        op isa QAdd || return leaf
        r = resolve_moment_sym(reg.moments, op)
        r === nothing || return r
        a = resolve_moment_sym(amb, op)
        return a === nothing ? leaf : a
    end
    new_eqs = Vector{Symbolics.Equation}(undef, length(eqs.equations))
    ps_set = Set{SymbolicUtils.BasicSymbolic}()
    @inbounds for (i, eq) in enumerate(eqs.equations)
        rhs = mapleaves(resolve, eq.rhs)
        new_eqs[i] = D(reg.vars[i]) ~ rhs
        _collect_params!(ps_set, SymbolicUtils.unwrap(rhs), iv_uw)
    end
    ps = [MTK.toparam(p) for p in ps_set]
    return MTK.System(new_eqs, iv, reg.vars, ps; name = name)
end

"""
    scale(c::CorrelationFunction)

Scale both the τ-equations and the underlying `eqs0` so the spectrum's
steady-state lookup resolves against the same scaled state names.
"""
scale(c::CorrelationFunction) = CorrelationFunction(
    c.op1, c.op2, c.op2_ancilla, c.aon_ancilla, scale(c.eqs), scale(c.eqs0), c.τ, c.steady_state,
)

# ---- Spectrum ----------------------------------------------------------------

"""
    Spectrum(c::CorrelationFunction, ps)

The Laplace-transformed correlation as a callable: `S(ω, u_end, p0)` returns the
symmetric power spectrum `2·Re{∫₀^∞ g(τ)e^{-iωτ}dτ}`. Solves `(iω·I - A)X̃ = x̃(0)`
where `A` is the τ-system Jacobian at steady state and `x̃ = x - x_∞` is centred so
the one-sided transform converges. The linear/anti-linear split (the conjugate
columns from `get_adjoints=false`) is handled by the 2n augmentation.

# Examples
```jldoctest
julia> h = FockSpace(:cavity);

julia> @qnumbers a::Destroy(h);

julia> eqs = meanfield([a' * a], a' * a, [a]; rates = [1.0], order = 2);

julia> c = CorrelationFunction(a', a, eqs);

julia> Spectrum(c)
ℱ(⟨a' * a⟩)(ω)
```
"""
struct Spectrum{C <: CorrelationFunction, P}
    c::C
    ω::Symbolics.Num
    ps::P
end
function Spectrum(c::CorrelationFunction, ps)
    c.steady_state || throw(
        ArgumentError(
            "Spectrum requires `CorrelationFunction(...; steady_state = true)`. " *
                "Got `steady_state = false`."
        ),
    )
    return Spectrum{typeof(c), typeof(ps)}(c, first(@variables ω_spectrum), ps)
end
Spectrum(c::CorrelationFunction) = Spectrum(c, ())

function (S::Spectrum)(ω_vals::AbstractVector, u_end, p0)
    A, rhs_b, n = _spectrum_kernel(S, u_end, p0)
    I_n = Matrix{ComplexF64}(I, n, n)
    out = Vector{Float64}(undef, length(ω_vals))
    @inbounds for (k, ω_val) in enumerate(ω_vals)
        x = ((im * ω_val) .* I_n .- A) \ rhs_b
        out[k] = 2 * real(x[1])
    end
    return out
end
(S::Spectrum)(ω_val::Real, u_end, p0) = first(S([ω_val], u_end, p0))

"""
Assemble the linear system `(iω·I - A)X̃ = x̃(0)` for the spectrum. The τ-system is
linear, so its Jacobian `A` is computed exactly by symbolic finite differences:
each RHS leaf is classified once as linear (matches a state), anti-linear
(conjugate of a state), or ambient (a steady-state coefficient). Probing the
linear leaves of state `j` gives column `j` of `A_lin`, the anti-linear leaves
give `A_alin`, and the 2n augmentation closes the conjugate response.
"""
function _spectrum_kernel(S::Spectrum, u_end, p0)
    c = S.c
    eqs = c.eqs
    n = length(eqs.equations)
    rhss_u = [SymbolicUtils.unwrap(eq.rhs) for eq in eqs.equations]
    p_sub = _build_p_sub(S.ps, p0, c.τ, u_end, eqs)
    resolve = _ss_resolver(c, u_end)

    # Classify each leaf by Hermitian conjugation, recording its side: same side as the
    # matched state is linear (A_lin), opposite side is anti-linear (A_alin). The
    # conjugate representative stays robust for scaled states, where adjoint does not round-trip.
    spec_ctx = build_ctx(eqs)
    state_map = MomentMap(
        spec_ctx, _treatments(eqs, spec_ctx),
        [undo_average(s) for s in eqs.states],
        collect(1:n),
    )
    rhs_leaves = SymbolicUtils.BasicSymbolic[]
    seen = Set{SymbolicUtils.BasicSymbolic}()
    for rhs in rhss_u, l in eachleaf(rhs)
        l in seen || (push!(seen, l); push!(rhs_leaves, l))
    end
    lin_by_state = Dict{Int, Vector{SymbolicUtils.BasicSymbolic}}()
    alin_by_state = Dict{Int, Vector{SymbolicUtils.BasicSymbolic}}()
    classified = Set{SymbolicUtils.BasicSymbolic}()
    for leaf in rhs_leaves
        r = match_moment(state_map, undo_average(leaf))
        r === nothing && continue
        i, same = r
        bucket = same ? lin_by_state : alin_by_state
        push!(get!(bucket, i, SymbolicUtils.BasicSymbolic[]), leaf)
        push!(classified, leaf)
    end
    for avg in _ambient_avgs(c)
        avg in classified && continue
        aons = SQA.acts_on(avg)
        lookup = (c.aon_ancilla in aons && length(aons) == 1) ? _undo_ancilla(c, avg) : avg
        p_sub[avg] = resolve(lookup)
    end

    zero_sub = copy(p_sub)
    for (_, leaves) in lin_by_state, l in leaves
        zero_sub[l] = 0.0 + 0.0im
    end
    for (_, leaves) in alin_by_state, l in leaves
        zero_sub[l] = 0.0 + 0.0im
    end

    b_const = ComplexF64[_scalarize(SymbolicUtils.substitute(rhs, zero_sub)) for rhs in rhss_u]

    A_lin = zeros(ComplexF64, n, n)
    for (j, leaves) in lin_by_state
        sub_j = copy(zero_sub)
        for l in leaves
            sub_j[l] = 1.0 + 0.0im
        end
        for (i, rhs) in enumerate(rhss_u)
            A_lin[i, j] = _scalarize(SymbolicUtils.substitute(rhs, sub_j)) - b_const[i]
        end
    end
    u_τ = ComplexF64[
        resolve(_undo_ancilla(c, SymbolicUtils.unwrap(s))) for s in eqs.states
    ]
    if isempty(alin_by_state)
        rhs_b = any(!iszero, b_const) ? u_τ .+ (A_lin \ b_const) : u_τ
        return A_lin, rhs_b, n
    end
    A_alin = zeros(ComplexF64, n, n)
    for (j, leaves) in alin_by_state
        sub_j = copy(zero_sub)
        for l in leaves
            sub_j[l] = 1.0 + 0.0im
        end
        for (i, rhs) in enumerate(rhss_u)
            A_alin[i, j] = _scalarize(SymbolicUtils.substitute(rhs, sub_j)) - b_const[i]
        end
    end
    M = [A_lin A_alin; conj.(A_alin) conj.(A_lin)]
    b_aug = vcat(b_const, conj.(b_const))
    u_aug = vcat(u_τ, conj.(u_τ))
    rhs_b = any(!iszero, b_aug) ? u_aug .+ (M \ b_aug) : u_aug
    return M, rhs_b, 2n
end

function _build_p_sub(ps, p0, τ, u_end, eqs)
    sub = Dict{Any, Any}()
    if ps !== nothing && !isempty(ps)
        for (p, v) in zip(ps, p0)
            sub[SymbolicUtils.unwrap(p)] = ComplexF64(v)
        end
    end
    if u_end isa AbstractDict
        for (k, v) in u_end
            ku = SymbolicUtils.unwrap(k)
            sub[ku] = ComplexF64(v)
            SQA.is_average(ku) && (sub[SymbolicUtils.unwrap(_avg_conj_of(ku))] = ComplexF64(conj(v)))
        end
    elseif u_end isa AbstractVector
        for (avg, v) in zip(eqs.states, u_end)
            sub[SymbolicUtils.unwrap(avg)] = ComplexF64(v)
        end
    end
    sub[SymbolicUtils.unwrap(τ)] = 0.0
    return sub
end
