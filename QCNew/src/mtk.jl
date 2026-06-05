# MTK bridge (Layer 6). Lower a (completed/scaled/evaluated) struct to a
# ModelingToolkit `System`. Because the struct is self-consistent (drift leaves
# carry the same naming regime as the states), codegen is one `mapleaves` pass
# that resolves each leaf to a state variable by `literal_key` (commuting-order
# and NE normalised, but NOT alpha-renamed, so concrete per-site atoms stay
# distinct), with a conjugate tier for the partner reps that `get_adjoints=false`
# leaves on the RHS. No per-level key juggling and no re-canonicalisation pass.

# ---- stable variable naming (replaces _stable_avg_name) ----------------------

serialize(op::QAdd) = Symbol("avg_" * _op_name_chunk(op))

function _op_name_chunk(op::QAdd)
    isempty(op.arguments) && return "zero"
    chunks = String[]
    for (term, _) in op.arguments
        base = join((_op_name_chunk(o) for o in term.ops), "_")
        if !isempty(term.ne)
            base *= "_" * join((string(a.name) * "neq" * string(b.name) for (a, b) in term.ne), "_")
        end
        push!(chunks, base)
    end
    return join(chunks, "_plus_")
end
_op_name_chunk(op::SQA.QSym) = string(op.name) * _op_name_extra(op)

_op_name_extra(op::SQA.QSym) = _op_index_suffix(op)
function _op_name_extra(op::SQA.Transition)
    i = op.i isa Symbol ? string(op.i) : string(Int(op.i))
    j = op.j isa Symbol ? string(op.j) : string(Int(op.j))
    return "_" * i * j * _op_index_suffix(op)
end
_op_name_extra(op::SQA.Destroy) = _op_index_suffix(op)
_op_name_extra(op::SQA.Create) = "_dag" * _op_index_suffix(op)
_op_name_extra(op::SQA.Pauli) = "_" * string(Int(op.axis)) * _op_index_suffix(op)
_op_name_extra(op::SQA.Spin) = "_" * string(Int(op.axis)) * _op_index_suffix(op)

function _op_index_suffix(op::SQA.QSym)
    isdefined(op, :index) || return ""
    idx = op.index
    idx === SQA.NO_INDEX && return ""
    return "_" * string(idx.name)
end

_make_var(name::Symbol, iv::Symbolics.Num) = first(@variables $name(iv))

# ---- state variable registry -------------------------------------------------

# Build the per-state variables and the literal/canon lookup maps shared by
# System / initial_values / get_solution, so naming agrees across them.
function _state_registry(eqs::AbstractMeanFieldEquations)
    ctx = build_ctx(eqs.operators, eqs.hamiltonian, eqs.jumps, eqs.jumps_dagger)
    # The system's recorded per-subspace coordinate (spec "The single source of
    # truth"): the resolver MUST key/match in this coordinate, not a hardcoded
    # `concrete_rep`. An empty map (scalar systems) reads as all-Free.
    coords = isempty(eqs.coords) ? all_free_coords(ctx) :
        Dict{Int, Coordinate}(sp => Coordinate(c) for (sp, c) in eqs.coords)
    ops = QAdd[undo_average(s) isa QAdd ? undo_average(s) : undo_average(s) * 1 for s in eqs.states]
    vars = Symbolics.Num[_make_var(serialize(op), eqs.iv) for op in ops]
    # ONE map keyed by the conjugation orbit rep in the system coordinate, valued
    # by `(state var, side-of-rep)`. A leaf on the SAME side as the stored rep
    # resolves to the var; the opposite side to its conjugate.
    by_rep = Dict{QAdd, Tuple{Symbolics.Num, Bool}}()
    by_canon = Dict{QAdd, Symbolics.Num}()
    for (op, v) in zip(ops, vars)
        rep, side = canonical_rep(op, ctx; coords)
        get!(by_rep, rep, (v, side))
        get!(by_canon, canon_key(op, ctx), v)
    end
    return (; ctx, coords, ops, vars, by_rep, by_canon)
end

# Resolve a RHS leaf to a state variable (or its conjugate), else leave it (an
# ambient parameter / steady-state coefficient). One orbit lookup plus a sign:
# the leaf and the stored rep agree iff they sit on the same conjugation side.
#
# The conjugate partner is emitted as a raw `term(conj, v; type=Number)` rather
# than `Base.conj(v)`. The state vars carry Real symtype (so only genuine
# parameters survive `_collect_params!`), and `Base.conj` on a `SymReal` folds
# to identity, which would silently zero every `⟨X⟩ - ⟨X†⟩` driving term. The
# raw term skips that simplifier and the `conj` node survives mtkcompile.
function _leaf_resolver(reg)
    return function (leaf)
        op = undo_average(leaf)
        op isa QAdd || return leaf
        rep, side = canonical_rep(op, reg.ctx; coords = reg.coords)
        haskey(reg.by_rep, rep) || return leaf
        var, rep_side = reg.by_rep[rep]
        return side == rep_side ? SymbolicUtils.unwrap(var) :
            SymbolicUtils.term(conj, SymbolicUtils.unwrap(var); type = Number)
    end
end

# Collect bare Real-symtype parameter symbols from a substituted RHS. Skips
# constants, the iv, time-dependent state vars (`u(iv)` calls), and average
# leaves (matched leaves are already vars; unmatched ambient averages are not
# parameters and are handled by the correlation layer).
function _collect_params!(set, x, iv_uw)
    x isa SymbolicUtils.BasicSymbolic || return
    SymbolicUtils.isconst(x) && return
    if !SymbolicUtils.iscall(x)
        (x !== iv_uw && SymbolicUtils.symtype(x) <: Real) && push!(set, x)
        return
    end
    _is_avg_leaf(x) && return
    op = SymbolicUtils.operation(x)
    args = SymbolicUtils.arguments(x)
    # Array-element access `δ[k]`: collect the array BASE symbol (symtype
    # `Vector{Real}`/`Matrix{Real}`) as a parameter, so the array-typed coupling
    # is bound (spec Task 7b; the scalar-only branch above skips it because its
    # symtype is not <: Real). MTK scalarizes the array param at compile.
    if op === getindex && !isempty(args) && args[1] isa SymbolicUtils.BasicSymbolic &&
            SymbolicUtils.symtype(args[1]) <: AbstractArray
        push!(set, args[1])
        return
    end
    (length(args) == 1 && args[1] === iv_uw && op isa SymbolicUtils.BasicSymbolic) && return
    for a in args
        _collect_params!(set, a, iv_uw)
    end
    return
end

# ---- System ------------------------------------------------------------------

"""
    ModelingToolkitBase.System(eqs::MeanFieldEquations; name)
    ModelingToolkitBase.System(eqs::NoiseMeanFieldEquations; name)

Lower the equation set to an MTK `System`: one `u(t)` variable per state, the
drift resolved to those variables. A `NoiseMeanFieldEquations` lowers to an SDE
whose single Brownian column is the aggregated noise drift (sign +1 Forward, -1
Backward). RHS leaves not matching a state become parameters.
"""
function MTK.System(eqs::MeanFieldEquations; name::Symbol)
    iv = eqs.iv
    iv_uw = SymbolicUtils.unwrap(iv)
    D = Symbolics.Differential(iv)
    reg = _state_registry(eqs)
    resolve = _leaf_resolver(reg)
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

MTK.System(eqs::NoiseMeanFieldEquations{O, H, Op, Jt, Jdt, R, E, S, Forward}; name::Symbol) where {O, H, Op, Jt, Jdt, R, E, S} =
    _to_system_sde(eqs, name, +1)
MTK.System(eqs::NoiseMeanFieldEquations{O, H, Op, Jt, Jdt, R, E, S, Backward}; name::Symbol) where {O, H, Op, Jt, Jdt, R, E, S} =
    _to_system_sde(eqs, name, -1)

function _to_system_sde(eqs::NoiseMeanFieldEquations, name::Symbol, sign::Int)
    iv = eqs.iv
    iv_uw = SymbolicUtils.unwrap(iv)
    D = Symbolics.Differential(iv)
    reg = _state_registry(eqs)
    resolve = _leaf_resolver(reg)
    w = first(MTK.@brownians _qc_dW)
    w_uw = SymbolicUtils.unwrap(w)
    T = typeof(w_uw)
    new_eqs = Vector{Symbolics.Equation}(undef, length(eqs.equations))
    ps_set = Set{SymbolicUtils.BasicSymbolic}()
    @inbounds for (i, eq) in enumerate(eqs.equations)
        rhs = SymbolicUtils.unwrap(mapleaves(resolve, eq.rhs))
        noise_rhs = SymbolicUtils.unwrap(mapleaves(resolve, eqs.noise_equations[i].rhs))
        signed_rhs = sign == 1 ? rhs : TermInterface.maketerm(T, *, Any[sign, rhs], nothing)
        noise_term = TermInterface.maketerm(T, *, Any[noise_rhs, w_uw], nothing)
        total_rhs = TermInterface.maketerm(T, +, Any[signed_rhs, noise_term], nothing)
        new_eqs[i] = D(reg.vars[i]) ~ total_rhs
        _collect_params!(ps_set, rhs, iv_uw)
        _collect_params!(ps_set, noise_rhs, iv_uw)
    end
    ps = [MTK.toparam(p) for p in ps_set]
    return MTK.System(new_eqs, iv, reg.vars, ps, [w]; name = name)
end

# ---- initial values / parameter map / solution ------------------------------

"""
    initial_values(eqs::AbstractMeanFieldEquations; defaults=Dict())

Map every state's `u(t)` variable to its initial value, defaulting to
`zero(ComplexF64)` for averages absent from `defaults`.
"""
function initial_values(eqs::AbstractMeanFieldEquations; defaults::AbstractDict = Dict())
    reg = _state_registry(eqs)
    out = Dict{Symbolics.Num, ComplexF64}()
    for (k, avg) in enumerate(eqs.states)
        out[reg.vars[k]] = ComplexF64(get(defaults, avg, 0))
    end
    return out
end

"""
    initial_values(eqs, state)            # from a numeric quantum state
    initial_values(eqs, u0::AbstractVector)  # from a value vector aligned with eqs.states
"""
function initial_values(eqs::AbstractMeanFieldEquations, state)
    return ComplexF64[ComplexF64(SQA.numeric_average(v, state)) for v in eqs.states]
end

function initial_values(eqs::AbstractMeanFieldEquations, u0::AbstractVector{<:Number}; kwargs...)
    length(u0) == length(eqs.states) || throw(
        DimensionMismatch(
            "initial value vector length $(length(u0)) does not match number of states $(length(eqs.states))",
        ),
    )
    reg = _state_registry(eqs)
    return Dict{Symbolics.Num, ComplexF64}(reg.vars[k] => ComplexF64(u0[k]) for k in eachindex(u0))
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
    parameter_map(eqs::AbstractMeanFieldEquations, pairs) -> Dict

Translate a user-facing `Pair`/`Dict` of parameter assignments into an
MTK-compatible map. An `IndexedVariable` callable key (e.g. `g` from
`g(i) = IndexedVariable(:g, i)`) is matched by name to the Symbolics array
parameter that `evaluate` synthesised from the per-atom references; a scalar
value is broadcast to fill that array, an `AbstractArray` value passes through
as the explicit per-atom values. Scalar (non-array) parameters pass through.
"""
function parameter_map(eqs::AbstractMeanFieldEquations, pairs)
    arr_by_name = Dict{Symbol, SymbolicUtils.BasicSymbolic}()
    scalar_by_name = Dict{Symbol, SymbolicUtils.BasicSymbolic}()
    for eq in eqs.equations
        _collect_named_params!(arr_by_name, scalar_by_name, SymbolicUtils.unwrap(eq.rhs))
        _collect_named_params!(arr_by_name, scalar_by_name, SymbolicUtils.unwrap(eq.lhs))
    end
    pmap = Dict{Any, Any}()
    # Per-slot accumulator: separate keys `δ(i_2_1)=>v1, δ(i_2_2)=>v2, …` each
    # carry a CONCRETE index (parseable slot), so they fill distinct slots of the
    # `δ` array. Without this each scalar `v` was broadcast over the whole array
    # via `fill`, so only the LAST assignment survived (every mode got the same
    # detuning). Bucket them by array name and assemble the array once at the end.
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
                shp = SymbolicUtils.shape(arr)
                n = (shp isa SymbolicUtils.SmallVec{UnitRange{Int}} && !isempty(shp)) ?
                    length(shp[1]) : nothing
                pmap[arr] = n === nothing ? v : fill(v, n)
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

# Concrete-index slots of an indexed-variable key, e.g. `δ(i_2_1)` -> `[1]`,
# `Γ(i_2_1, i_2_2)` -> `[1, 2]`. Returns `nothing` when any argument is not a
# concrete `name_…_<int>` index (a bare symbolic index `δ(i)` has no slot, so it
# broadcasts instead of targeting a single position).
function _param_slots(x::SymbolicUtils.BasicSymbolic)
    SymbolicUtils.iscall(x) || return nothing
    slots = Int[]
    for a in SymbolicUtils.arguments(x)
        (a isa SymbolicUtils.BasicSymbolic && !SymbolicUtils.iscall(a)) || return nothing
        s = _parse_slot(Base.nameof(a))
        s === nothing && return nothing
        push!(slots, s)
    end
    return isempty(slots) ? nothing : slots
end
_param_slots(_) = nothing

# Assemble a concrete array parameter from per-slot scalar assignments. Sizes from
# the synthesised array's declared shape when available, else from the maximum slot
# per dimension. Unassigned slots stay zero (the caller is expected to specify
# every active site, as the evaluate-time array shape implies).
function _build_slot_array(arr, slots_to_v::Dict{Vector{Int}, Any})
    ndim = length(first(keys(slots_to_v)))
    shp = SymbolicUtils.shape(arr)
    dims = (shp isa SymbolicUtils.SmallVec{UnitRange{Int}} && length(shp) == ndim) ?
        Tuple(length(r) for r in shp) :
        Tuple(maximum(s[d] for s in keys(slots_to_v)) for d in 1:ndim)
    V = promote_type((typeof(v) for v in values(slots_to_v))...)
    out = zeros(V, dims...)
    for (slots, v) in slots_to_v
        out[slots...] = v
    end
    return out
end

# User-visible name of a parameter expression: bare Sym `:κ` -> :κ; callable-Sym
# Term `g(i.sym)` (IndexedVariable shape) -> :g; nothing otherwise.
_param_name(x) = nothing
function _param_name(x::SymbolicUtils.BasicSymbolic)
    SymbolicUtils.iscall(x) || return Base.nameof(x)
    op = SymbolicUtils.operation(x)
    (op isa SymbolicUtils.BasicSymbolic && !SymbolicUtils.iscall(op)) && return Base.nameof(op)
    return nothing
end

# Populate name -> MTK-parameter dicts from a tree: array-typed Syms into
# `arr_dict`, bare Real scalar Syms into `scalar_dict`.
function _collect_named_params!(arr_dict, scalar_dict, x)
    x isa SymbolicUtils.BasicSymbolic || return
    if !SymbolicUtils.iscall(x)
        SymbolicUtils.issym(x) || return
        shp = SymbolicUtils.shape(x)
        if shp isa SymbolicUtils.SmallVec{UnitRange{Int}} && !isempty(shp)
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
averaged internally). Resolves the user's operator to a tracked state by the
conjugation orbit (`concrete_rep`, recovering a folded conjugate via the side
bit), then falls back to `canon_key` so a query posed in a different but
equivalent symbolic form still hits.
"""
function get_solution(sol, avg::SymbolicUtils.BasicSymbolic, eqs::AbstractMeanFieldEquations)
    reg = _state_registry(eqs)
    op = undo_average(avg)
    if op isa QAdd
        rep, side = canonical_rep(op, reg.ctx; coords = reg.coords)
        if haskey(reg.by_rep, rep)
            var, rep_side = reg.by_rep[rep]
            return side == rep_side ? (τ -> _eval_at(sol, var, τ)) :
                (τ -> conj.(_eval_at(sol, var, τ)))
        end
        haskey(reg.by_canon, canon_key(op, reg.ctx)) && return τ -> _eval_at(sol, reg.by_canon[canon_key(op, reg.ctx)], τ)
        kcc = canon_key(adjoint(op), reg.ctx)
        haskey(reg.by_canon, kcc) && return τ -> conj.(_eval_at(sol, reg.by_canon[kcc], τ))
    end
    throw(KeyError(avg))
end
get_solution(sol, op::QField, eqs::AbstractMeanFieldEquations) = get_solution(sol, average(op), eqs)

_eval_at(sol, var, τ::AbstractVector) = Array(sol(τ; idxs = var))
_eval_at(sol, var, τ) = sol(τ; idxs = var)
