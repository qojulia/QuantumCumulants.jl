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
- `ambient`: the τ-system's ambient averages (steady-state coefficients, not τ-states),
  computed once from the closed `eqs` so the `u0`/`System`/spectrum paths need not re-walk it.
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
    ambient::Vector{SymbolicUtils.BasicSymbolic}
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
`eqs0` via the quantum regression theorem. `op2` is re-embedded onto a fresh ancilla
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

    ord = eqs0.order
    order_ext = ord === nothing ? nothing :
        vcat(ord, fill(maximum(ord), aon_ancilla - length(ord)))

    eqs_c = meanfield(
        [new_op], eqs0.hamiltonian, eqs0.jumps;
        Jdagger = eqs0.jumps_dagger, rates = eqs0.rates, order = order_ext, iv = τ,
    )

    ancilla_filter = _ancilla_filter(aon_ancilla)
    closure_filter = filter_func === nothing ? ancilla_filter : x -> ancilla_filter(x) && filter_func(x)
    # `meanfield` seeds (does not close) the ancilla system, so `eqs_c.graph` is the seeded
    # graph; close it here, purely.
    g = closure(eqs_c.graph; filter = closure_filter, get_adjoints = false, max_iter)
    if !steady_state
        g = closure(g; filter = x -> !ancilla_filter(x), get_adjoints = false, max_iter)
    end
    if !isempty(eqs0.graph.treatments)
        merged = merge(g.treatments, eqs0.graph.treatments)
        g = MomentGraph(g.nodes, g.sys, g.ctx, merged)
    end
    if filter_func !== nothing
        g = map_drifts(g, (_, d) -> _filter_ancilla_expr(d, ancilla_filter, filter_func))
    end
    eqs_tau = MeanfieldEquations(g)
    ambient = _ambient_avgs(eqs_tau, aon_ancilla)

    return CorrelationFunction(
        op1, op2, op2_ancilla, aon_ancilla, eqs_tau, eqs0, τ, steady_state, ambient,
    )
end

function _filter_ancilla_expr(x, ancilla_filter, user_filter)
    return rewrite(x) do y
        _is_avg_leaf(y) ? ((ancilla_filter(y) && !user_filter(y)) ? 0 : y) : nothing
    end
end

# ---- ambient / lookup helpers ------------------------------------------------

"""
Averages on the τ-system's RHS that are not τ-states: non-ancilla coefficients or
pure-ancilla constants. These become MTK parameters in `System(c)`. Computed once at
construction and stored in `c.ambient`; the `u0`/`System`/spectrum paths read that field.
"""
function _ambient_avgs(eqs::MeanfieldEquations, aon_ancilla::Int)
    found = OrderedCollections.OrderedSet{SymbolicUtils.BasicSymbolic}()
    states_set = Set(SymbolicUtils.unwrap.(eqs.states))
    for eq in eqs.equations
        _collect_ambient!(found, eq.rhs, aon_ancilla, states_set)
    end
    return collect(found)
end
function _collect_ambient!(found, x, aon_ancilla, states_set)
    walk(x) do n
        if SQA.is_average(n) && SymbolicUtils.operation(n) === SQA.sym_average
            if !(n in states_set)
                aons = SQA.acts_on(n)
                (!(aon_ancilla in aons) || length(aons) == 1) && push!(found, n)
            end
            return false
        end
        return true
    end
    return
end

_ambient_param(avg::SymbolicUtils.BasicSymbolic) =
    first(@variables $(Symbol("ss_", avg_name(undo_average(avg)))))

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
    # Keys are unwrapped `Number`-symtype time-dependent averages, not `Num`-wrappable.
    u0 = Dict{Any, ComplexF64}()
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
    for avg in c.ambient
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
    amb_avgs = c.ambient
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
    new_eqs = Vector{Symbolics.Equation}(undef, length(eqs))
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
function scale(c::CorrelationFunction)
    eqs = scale(c.eqs)
    # `ambient` is moment-set dependent, so recompute it for the scaled τ-system rather than
    # carrying the unscaled set.
    return CorrelationFunction(
        c.op1, c.op2, c.op2_ancilla, c.aon_ancilla, eqs, scale(c.eqs0), c.τ, c.steady_state,
        _ambient_avgs(eqs, c.aon_ancilla),
    )
end
