serialize(op::QAdd) = Symbol("avg_" * _op_name_part(op))

function _op_name_part(op::QAdd)
    isempty(op.arguments) && return "zero"
    chunks = String[]
    for (term, _) in op.arguments
        base = join((_op_name_part(o) for o in term.ops), "_")
        if !isempty(term.ne)
            base *= "_" * join((string(a.name) * "neq" * string(b.name) for (a, b) in term.ne), "_")
        end
        push!(chunks, base)
    end
    return join(chunks, "_plus_")
end
_op_name_part(op::SQA.QSym) = string(op.name) * _op_name_suffix(op)

_op_name_suffix(op::SQA.QSym) = _op_index_suffix(op)
function _op_name_suffix(op::SQA.Transition)
    i = op.i isa Symbol ? string(op.i) : string(Int(op.i))
    j = op.j isa Symbol ? string(op.j) : string(Int(op.j))
    return "_" * i * j * _op_index_suffix(op)
end
_op_name_suffix(op::SQA.Destroy) = _op_index_suffix(op)
_op_name_suffix(op::SQA.Create) = "_dag" * _op_index_suffix(op)
_op_name_suffix(op::SQA.Pauli) = "_" * string(Int(op.axis)) * _op_index_suffix(op)
_op_name_suffix(op::SQA.Spin) = "_" * string(Int(op.axis)) * _op_index_suffix(op)

function _op_index_suffix(op::SQA.QSym)
    isdefined(op, :index) || return ""
    idx = op.index
    idx === SQA.NO_INDEX && return ""
    return "_" * string(idx.name)
end

_make_var(name::Symbol, iv::Symbolics.Num) = first(@variables $name(iv))

# ---- state variable registry -------------------------------------------------

"""
Build the per-state `u(t)` variables and the scaled/`canon_key` lookup maps shared by
`System`, `initial_values` and `get_solution`, so variable naming agrees across them.
"""
function _state_registry(eqs::AbstractMeanFieldEquations)
    ctx = build_ctx(eqs)
    # Key/match in the system's recorded treatment, not a hardcoded `concrete_rep`
    # (an empty map, scalar systems, reads as all-Free).
    treatments = _treatments(eqs, ctx)
    ops = QAdd[undo_average(s) isa QAdd ? undo_average(s) : undo_average(s) * 1 for s in eqs.states]
    vars = Symbolics.Num[_make_var(serialize(op), eqs.iv) for op in ops]
    # Keyed by the Hermitian-conjugate representative, valued by `(state var, side-of-rep)`: a leaf
    # on the same side as the rep resolves to the var, the opposite side to its conjugate.
    by_rep = Dict{QAdd, Tuple{Symbolics.Num, Bool}}()
    by_canon = Dict{QAdd, Symbolics.Num}()
    for (op, v) in zip(ops, vars)
        rep, side = canonical_rep(op, ctx; treatments)
        get!(by_rep, rep, (v, side))
        get!(by_canon, canon_key(op, ctx), v)
    end
    return (; ctx, treatments, ops, vars, by_rep, by_canon)
end

"""
Resolve a `(rep, side)` against a `rep → (symbol, rep_side)` table: the symbol on the
matching conjugation side, its conjugate on the opposite side, or `nothing` when `rep`
is absent. Raw `term(conj, …)` not `Base.conj`, which folds to identity on Real symtype.
"""
function _resolve_side(table, rep, side)
    haskey(table, rep) || return nothing
    sym, rep_side = table[rep]
    u = SymbolicUtils.unwrap(sym)
    return side == rep_side ? u : SymbolicUtils.term(conj, u; type = Number)
end

"""
Build a closure that resolves a RHS average leaf to its state variable (or that
variable's conjugate), leaving non-state leaves untouched as ambient parameters. The
leaf and the stored representative agree iff they sit on the same conjugation side.
"""
function _leaf_resolver(reg)
    return function (leaf)
        op = undo_average(leaf)
        op isa QAdd || return leaf
        rep, side = canonical_rep(op, reg.ctx; treatments = reg.treatments)
        r = _resolve_side(reg.by_rep, rep, side)
        return r === nothing ? leaf : r
    end
end

"""
Collect the bare `Real`-symtype parameter symbols from a substituted RHS into `set`,
skipping constants, the independent variable, time-dependent state variables and
average leaves.
"""
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
    # Array access `δ[k]`: collect the array BASE symbol as a parameter (the scalar
    # branch above skips it, symtype not <: Real); MTK scalarizes it at compile.
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

Build an MTK `System` from the equation set: one `u(t)` variable per state, the
drift resolved to those variables. A `NoiseMeanFieldEquations` becomes an SDE
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
                shape = SymbolicUtils.shape(arr)
                n = (shape isa SymbolicUtils.SmallVec{UnitRange{Int}} && !isempty(shape)) ?
                    length(shape[1]) : nothing
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

"""
Concrete-index slots of an indexed-variable key, e.g. `δ(i_2_1)` becomes `[1]`. Returns
`nothing` when any argument is not a concrete `name_…_<int>` index, in which case the
value is broadcast over the array rather than targeting a single slot.
"""
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

"""
Assemble a concrete array parameter from per-slot scalar assignments. The size comes
from the synthesised array's declared shape when available, otherwise from the maximum
slot per dimension; unassigned slots stay zero.
"""
function _build_slot_array(arr, slots_to_v::Dict{Vector{Int}, Any})
    ndim = length(first(keys(slots_to_v)))
    shape = SymbolicUtils.shape(arr)
    dims = (shape isa SymbolicUtils.SmallVec{UnitRange{Int}} && length(shape) == ndim) ?
        Tuple(length(r) for r in shape) :
        Tuple(maximum(s[d] for s in keys(slots_to_v)) for d in 1:ndim)
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
        shape = SymbolicUtils.shape(x)
        if shape isa SymbolicUtils.SmallVec{UnitRange{Int}} && !isempty(shape)
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
Hermitian conjugation (`concrete_rep`, recovering a folded conjugate via the side
bit), then falls back to `canon_key` so a query posed in a different but
equivalent symbolic form still hits.
"""
function get_solution(sol, avg::SymbolicUtils.BasicSymbolic, eqs::AbstractMeanFieldEquations)
    reg = _state_registry(eqs)
    op = undo_average(avg)
    if op isa QAdd
        rep, side = canonical_rep(op, reg.ctx; treatments = reg.treatments)
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
