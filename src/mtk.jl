# Pretty display label (⟨a'*a⟩) for a moment's unknown; not the identity, dedup is structural.
avg_name(op::QAdd) = Symbol(string(SQA.average(op)))

# The moment's time-dependent `Number` unknown, lifted by SQA's `make_time_dependent`.
# `iv` is MTK's `t_nounits` (see `_make_iv`), so the result is a first-class MTK unknown.
time_dependent_var(op::QAdd, iv) = SymbolicUtils.unwrap(SQA.make_time_dependent(SQA.average(op), iv))

_state_vars(ops::AbstractVector, iv) = SymbolicUtils.BasicSymbolic[time_dependent_var(op, iv) for op in ops]

# ---- state variable registry -------------------------------------------------

"""
Build the per-state `u(t)` variables and the `moments` lookup shared by `System`,
`initial_values` and `get_solution`, so variable naming agrees across them. `moments`
keys each state var by its Hermitian-conjugate representative; resolution (leaves,
`get_solution`) goes through that one structural matcher.
"""
function _state_registry(eqs::AbstractMeanfieldEquations)
    ctx = build_ctx(eqs)
    # Key/match in the system's recorded treatment, not a hardcoded `concrete_rep`
    # (an empty map, scalar systems, reads as all-Free).
    treatments = _treatments(eqs, ctx)
    ops = QAdd[(o = undo_average(s); o isa QAdd ? o : o * 1) for s in eqs.states]
    # Lifted averages are `Number`-symtype, which `Symbolics.Num` (requires `<:Real`)
    # cannot wrap, so the state variables are carried as unwrapped `BasicSymbolic`.
    vars = _state_vars(ops, eqs.iv)
    moments = MomentMap(ctx, treatments, ops, vars)
    return (; ctx, treatments, vars, moments)
end

"""
    moment_variable_map(eqs::AbstractMeanfieldEquations)

Ordered map from each moment average ``⟨op⟩`` to the time-dependent MTK unknown
`u(t)` that [`System`](@ref) assigns it. The bridge `initial_values`, `get_solution`
and `parameter_map` resolve moments against.
"""
function moment_variable_map(eqs::AbstractMeanfieldEquations)
    reg = _state_registry(eqs)
    return OrderedCollections.OrderedDict(
        eqs.states[i] => reg.vars[i] for i in eachindex(eqs.states)
    )
end

"""
Build a closure that resolves a RHS average leaf to its state variable (or that
variable's conjugate). A closed system has every RHS moment among its states, so a
non-matching leaf means an unclosed system; it is left untouched (the caller is
expected to `complete` first). Correlation/spectrum handle their ambient
steady-state moments separately.
"""
function _leaf_resolver(reg)
    return function (leaf)
        op = undo_average(leaf)
        op isa QAdd || return leaf
        r = resolve_moment_sym(reg.moments, op)
        return r === nothing ? leaf : r
    end
end

"""
Collect the bare `Real`-symtype parameter symbols from a substituted RHS into `set`,
skipping constants, the independent variable, time-dependent state variables and
average leaves.
"""
function _collect_params!(set, x, iv_uw)
    walk(x) do n
        SymbolicUtils.isconst(n) && return false
        if !SymbolicUtils.iscall(n)
            (n !== iv_uw && SymbolicUtils.symtype(n) <: Real) && push!(set, n)
            return false
        end
        _is_avg_leaf(n) && return false
        op = SymbolicUtils.operation(n)
        args = SymbolicUtils.arguments(n)
        # Array access `δ[k]`: collect the array BASE symbol as a parameter (the scalar
        # branch above skips it, symtype not <: Real); MTK scalarizes it at compile.
        if op === getindex && !isempty(args) && args[1] isa SymbolicUtils.BasicSymbolic &&
                SymbolicUtils.symtype(args[1]) <: AbstractArray
            push!(set, args[1])
            return false
        end
        (length(args) == 1 && args[1] === iv_uw && op isa SymbolicUtils.BasicSymbolic) && return false
        return true
    end
    return
end

function _sorted_params(ps_set)
    ps = collect(ps_set)
    return ps[sortperm(string.(ps))]
end

"""
    noise_channels(eqs::AbstractMeanfieldEquations)

The monitored measurement channels: one entry per jump with nonzero efficiency, each an
independent Brownian in the SDE [`System`](@ref) builds. Each entry is
`(; index, jump, rate, efficiency)`; a deterministic system has none.
"""
noise_channels(::MeanfieldEquations) = NamedTuple[]
function noise_channels(eqs::NoiseMeanfieldEquations)
    return [
        (; index = k, jump = eqs.jumps[k], rate = eqs.rates[k], efficiency = eqs.efficiencies[k])
            for k in eachindex(eqs.efficiencies) if !iszero(eqs.efficiencies[k])
    ]
end

function _noise_channel_rhss(eqs::NoiseMeanfieldEquations)
    active = Int[ch.index for ch in noise_channels(eqs)]
    channels = Vector{Vector{Any}}(undef, length(active))
    build_noise = _noise_builder(eqs.direction)
    for (j, k) in enumerate(active)
        _, noise_eqs = build_noise(
            eqs.operators,
            [eqs.jumps[k]],
            [eqs.jumps_dagger[k]],
            [eqs.rates[k]],
            [eqs.efficiencies[k]],
        )
        channels[j] = Any[
            SymbolicUtils.unwrap(
                    _reduce_ground_in_drift(
                        eqs.order === nothing ? Symbolics.Num(eq.rhs) :
                        Symbolics.Num(cumulant_expansion(eq.rhs, eqs.order))
                    ),
                ) for eq in noise_eqs
        ]
    end
    return active, channels
end

# Build a Brownian variable programmatically (mirrors the `@brownians` macro
# expansion) so System construction avoids a runtime `eval` per call.
function _make_brownian(name::Symbol)
    sym = SymbolicUtils.Sym{Symbolics.SymReal}(name; type = Real, shape = UnitRange{Int64}[])
    tagged = SymbolicUtils.setmetadata(sym, Symbolics.VariableSource, (:brownian, name))
    return MTK.tobrownian(Symbolics.wrap(tagged))
end
_make_brownians(names::Vector{Symbol}) = [_make_brownian(name) for name in names]

# ---- System ------------------------------------------------------------------

"""
    ModelingToolkitBase.System(eqs::MeanfieldEquations; name)
    ModelingToolkitBase.System(eqs::NoiseMeanfieldEquations; name)

Build an MTK `System` from the equation set: one `u(t)` variable per state, the
drift resolved to those variables. A `NoiseMeanfieldEquations` becomes an SDE
with one independent Brownian column per monitored channel, each carrying that
channel's noise drift (sign +1 Forward, -1 Backward). RHS leaves not matching a
state become parameters.
"""
function MTK.System(eqs::MeanfieldEquations; name::Symbol)
    iv = eqs.iv
    iv_uw = SymbolicUtils.unwrap(iv)
    D = Symbolics.Differential(iv)
    reg = _state_registry(eqs)
    resolve = _leaf_resolver(reg)
    new_eqs = Vector{Symbolics.Equation}(undef, length(eqs.equations))
    ps_set = Set{SymbolicUtils.BasicSymbolic}()
    @inbounds for (i, eq) in enumerate(eqs.equations)
        rhs = mapleaves(resolve, SymbolicUtils.unwrap(eq.rhs))
        new_eqs[i] = D(reg.vars[i]) ~ rhs
        _collect_params!(ps_set, rhs, iv_uw)
    end
    ps = [MTK.toparam(p) for p in _sorted_params(ps_set)]
    return MTK.System(new_eqs, iv, reg.vars, ps; name = name)
end

MTK.System(eqs::NoiseMeanfieldEquations{O, H, Op, Jt, Jdt, R, E, S, Forward}; name::Symbol) where {O, H, Op, Jt, Jdt, R, E, S} =
    _to_system_sde(eqs, name, +1)
MTK.System(eqs::NoiseMeanfieldEquations{O, H, Op, Jt, Jdt, R, E, S, Backward}; name::Symbol) where {O, H, Op, Jt, Jdt, R, E, S} =
    _to_system_sde(eqs, name, -1)

function _to_system_sde(eqs::NoiseMeanfieldEquations, name::Symbol, sign::Int)
    iv = eqs.iv
    iv_uw = SymbolicUtils.unwrap(iv)
    D = Symbolics.Differential(iv)
    reg = _state_registry(eqs)
    resolve = _leaf_resolver(reg)
    active, channel_noise = _noise_channel_rhss(eqs)
    # One independent Brownian per monitored channel. With no active channel keep a
    # single zero-noise column so the result is still a valid SDE.
    nch = max(length(active), 1)
    ws = _make_brownians([Symbol("_qc_dW_", j) for j in 1:nch])
    ws_uw = SymbolicUtils.unwrap.(ws)
    T = typeof(first(ws_uw))
    new_eqs = Vector{Symbolics.Equation}(undef, length(eqs.equations))
    ps_set = Set{SymbolicUtils.BasicSymbolic}()
    @inbounds for (i, eq) in enumerate(eqs.equations)
        rhs = mapleaves(resolve, SymbolicUtils.unwrap(eq.rhs))
        signed_rhs = sign == 1 ? rhs : TermInterface.maketerm(T, *, Any[sign, rhs], nothing)
        total_rhs = signed_rhs
        for (j, w_uw) in enumerate(ws_uw)
            noise_rhs = isempty(active) ? 0 :
                SymbolicUtils.unwrap(mapleaves(resolve, channel_noise[j][i]))
            noise_term = TermInterface.maketerm(T, *, Any[noise_rhs, w_uw], nothing)
            total_rhs = TermInterface.maketerm(T, +, Any[total_rhs, noise_term], nothing)
            _collect_params!(ps_set, noise_rhs, iv_uw)
        end
        new_eqs[i] = D(reg.vars[i]) ~ total_rhs
        _collect_params!(ps_set, rhs, iv_uw)
    end
    ps = [MTK.toparam(p) for p in _sorted_params(ps_set)]
    return MTK.System(new_eqs, iv, reg.vars, ps, ws; name = name)
end

# ---- initial values / parameter map / solution ------------------------------

"""
    initial_values(eqs::AbstractMeanfieldEquations; defaults=Dict())

Map every state's `u(t)` variable to its initial value, defaulting to
`zero(ComplexF64)` for averages absent from `defaults`.
"""
function initial_values(eqs::AbstractMeanfieldEquations; defaults::AbstractDict = Dict())
    reg = _state_registry(eqs)
    # Keys are unwrapped `Number`-symtype average variables, not `Num`-wrappable.
    out = Dict{Any, ComplexF64}()
    for (k, avg) in enumerate(eqs.states)
        out[reg.vars[k]] = ComplexF64(get(defaults, avg, 0))
    end
    return out
end

"""
    initial_values(eqs, state)            # from a numeric quantum state
    initial_values(eqs, u0::AbstractVector)  # from a value vector aligned with eqs.states
"""
function initial_values(eqs::AbstractMeanfieldEquations, state)
    return ComplexF64[ComplexF64(SQA.numeric_average(v, state)) for v in eqs.states]
end

function initial_values(eqs::AbstractMeanfieldEquations, u0::AbstractVector{<:Number}; kwargs...)
    length(u0) == length(eqs.states) || throw(
        DimensionMismatch(
            "initial value vector length $(length(u0)) does not match number of states $(length(eqs.states))",
        ),
    )
    reg = _state_registry(eqs)
    return Dict{Any, ComplexF64}(reg.vars[k] => ComplexF64(u0[k]) for k in eachindex(u0))
end

"""
    parameter_map(sys::MTK.System, pairs)

Keep only the entries of `pairs` whose key is a live unknown or parameter of the
compiled system `sys` (MTK rejects superfluous entries).
"""
function parameter_map(sys::MTK.System, pairs)
    live = Set{SymbolicUtils.BasicSymbolic}()
    for p in MTK.parameters(sys)
        push!(live, SymbolicUtils.unwrap(p))
    end
    for u in MTK.unknowns(sys)
        push!(live, SymbolicUtils.unwrap(u))
    end
    out = Dict{Any, Any}()
    for (k, v) in pairs
        SymbolicUtils.unwrap(k) in live && (out[k] = v)
    end
    return out
end

"""
    parameter_map(eqs::AbstractMeanfieldEquations, pairs) -> Dict

Translate a user-facing `Pair`/`Dict` of parameter assignments into an
MTK-compatible map. An `IndexedVariable` callable key (e.g. `g` from
`g(i) = IndexedVariable(:g, i)`) is matched by name to the Symbolics array
parameter that `evaluate` synthesised from the per-atom references; a scalar
value is broadcast to fill that array, an `AbstractArray` value passes through
as the explicit per-atom values. Scalar (non-array) parameters pass through.
"""
function parameter_map(eqs::AbstractMeanfieldEquations, pairs)
    arr_by_name = Dict{Symbol, SymbolicUtils.BasicSymbolic}()
    scalar_by_name = Dict{Symbol, SymbolicUtils.BasicSymbolic}()
    for eq in eqs.equations
        _collect_named_params!(arr_by_name, scalar_by_name, SymbolicUtils.unwrap(eq.rhs))
        _collect_named_params!(arr_by_name, scalar_by_name, SymbolicUtils.unwrap(eq.lhs))
    end
    pmap = Dict{Any, Any}()
    # Per-slot accumulator: keys like `δ(i_2_1)=>v1` carry a concrete slot, so each
    # fills a distinct array slot (a bare scalar broadcasts via `fill` instead).
    slot_acc = Dict{Symbol, Dict{Vector{Int}, Any}}()
    for (k, v) in pairs
        ku = SymbolicUtils.unwrap(k)
        name = _param_name(ku)
        if name !== nothing && haskey(arr_by_name, name)
            arr = arr_by_name[name]
            slots = _param_slots(ku)
            if v isa AbstractArray
                pmap[arr] = v
            elseif slots !== nothing
                get!(slot_acc, name, Dict{Vector{Int}, Any}())[slots] = v
            else
                dims = _array_dims(arr)
                pmap[arr] = dims === nothing ? v : fill(v, dims[1])
            end
        elseif name !== nothing && haskey(scalar_by_name, name)
            pmap[scalar_by_name[name]] = v
        else
            pmap[k] = v
        end
    end
    for (name, slots_to_v) in slot_acc
        pmap[arr_by_name[name]] = _build_slot_array(arr_by_name[name], slots_to_v)
    end
    return pmap
end

"""
Concrete-index slots of an indexed-variable key, e.g. `δ(i(1))` becomes `[1]`. Reads each
argument's slot from its `SQA.index_slot` metadata; returns `nothing` when any argument
carries no slot, in which case the value is broadcast over the array rather than targeting
a single slot.
"""
function _param_slots(x::SymbolicUtils.BasicSymbolic)
    SymbolicUtils.iscall(x) || return nothing
    slots = Int[]
    for a in SymbolicUtils.arguments(x)
        (a isa SymbolicUtils.BasicSymbolic && !SymbolicUtils.iscall(a)) || return nothing
        s = SQA.index_slot(a)
        s === nothing && return nothing
        push!(slots, s)
    end
    return isempty(slots) ? nothing : slots
end
_param_slots(_) = nothing

"""
Concrete per-axis lengths of a shaped array-variable sym, or `nothing` for scalars and
symbolically-sized arrays. The single place that inspects `SymbolicUtils.shape`'s layout
(scalars yield an empty shape, so emptiness distinguishes them from arrays).
"""
function _array_dims(x)
    sh = SymbolicUtils.shape(x)
    (sh isa SymbolicUtils.SmallVec{UnitRange{Int}} && !isempty(sh)) || return nothing
    return ntuple(d -> length(sh[d]), length(sh))
end

"""
Assemble a concrete array parameter from per-slot scalar assignments. The size comes
from the synthesised array's declared shape when available, otherwise from the maximum
slot per dimension; unassigned slots stay zero.
"""
function _build_slot_array(arr, slots_to_v::Dict{Vector{Int}, Any})
    ndim = length(first(keys(slots_to_v)))
    dims = _array_dims(arr)
    if dims === nothing || length(dims) != ndim
        dims = ntuple(d -> maximum(s[d] for s in keys(slots_to_v)), ndim)
    end
    V = promote_type((typeof(v) for v in values(slots_to_v))...)
    out = zeros(V, dims...)
    for (slots, v) in slots_to_v
        out[slots...] = v
    end
    return out
end

"""
User-visible name of a parameter expression: a bare `Sym` `κ` gives `:κ`, a callable
`g(i)` (the `IndexedVariable` shape) gives `:g`, and anything else gives `nothing`.
"""
_param_name(x) = nothing
function _param_name(x::SymbolicUtils.BasicSymbolic)
    SymbolicUtils.iscall(x) || return Base.nameof(x)
    op = SymbolicUtils.operation(x)
    (op isa SymbolicUtils.BasicSymbolic && !SymbolicUtils.iscall(op)) && return Base.nameof(op)
    return nothing
end

"""
Populate the name → parameter dictionaries from an expression tree: array-typed `Sym`s
go into `arr_dict`, bare `Real` scalar `Sym`s into `scalar_dict`.
"""
function _collect_named_params!(arr_dict, scalar_dict, x)
    x isa SymbolicUtils.BasicSymbolic || return
    if !SymbolicUtils.iscall(x)
        SymbolicUtils.issym(x) || return
        if SymbolicUtils.symtype(x) <: AbstractArray
            arr_dict[Base.nameof(x)] = x
        elseif SymbolicUtils.symtype(x) <: Real
            scalar_dict[Base.nameof(x)] = x
        end
        return
    end
    for a in SymbolicUtils.arguments(x)
        _collect_named_params!(arr_dict, scalar_dict, a)
    end
    return
end

"""
    get_solution(sol, avg_or_op, eqs)

Query an ODE/SDE solution for the trajectory of an average (or a `QField`,
averaged internally). Resolves the user's operator to a tracked state through the
single structural matcher `match_moment`, which folds Hermitian conjugation (the
side bit recovers a stored conjugate) and the system's index/symmetry treatment, so
a query posed in any equivalent symbolic form hits the same state.
"""
function get_solution(sol, avg::SymbolicUtils.BasicSymbolic, eqs::AbstractMeanfieldEquations)
    reg = _state_registry(eqs)
    op = undo_average(avg)
    if op isa QAdd
        r = match_moment(reg.moments, op)
        if r !== nothing
            var, same = r
            return same ? (τ -> _eval_at(sol, var, τ)) :
                (τ -> conj.(_eval_at(sol, var, τ)))
        end
    end
    throw(KeyError(avg))
end
get_solution(sol, op::QField, eqs::AbstractMeanfieldEquations) = get_solution(sol, average(op), eqs)

_eval_at(sol, var, τ::AbstractVector) = Array(sol(τ; idxs = var))
_eval_at(sol, var, τ) = sol(τ; idxs = var)
