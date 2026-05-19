function _stable_avg_name(avg::SymbolicUtils.BasicSymbolic)
    @assert SQA.is_average(avg) "expected an Average BasicSymbolic"
    op = SQA.undo_average(avg)
    s = "avg_" * _op_name_chunk(op)
    return Symbol(s)
end

function _op_name_chunk(op::QAdd)
    isempty(op.arguments) && return "zero"
    chunks = String[]
    for (term, _) in op.arguments
        base = join((_op_name_chunk(o) for o in term.ops), "_")
        # `ne` (non-equal-index constraints) is semantic state and must
        # appear in the MTK identifier — otherwise two terms that differ
        # only in their constraint set generate clashing variable names
        # (see TODO.md unique_squeezing dedup note).
        if !isempty(term.ne)
            ne_chunks = String[]
            for (a, b) in term.ne
                push!(ne_chunks, string(a.name) * "neq" * string(b.name))
            end
            base *= "_" * join(ne_chunks, "_")
        end
        push!(chunks, base)
    end
    return join(chunks, "_plus_")
end
function _op_name_chunk(op::SQA.QSym)
    base = string(op.name)
    extra = _op_name_extra(op)
    return isempty(extra) ? base : base * extra
end

function _op_name_extra(op::SQA.QSym)
    return _op_index_suffix(op)
end
function _op_name_extra(op::SQA.Transition)
    i = op.i isa Symbol ? string(op.i) : string(Int(op.i))
    j = op.j isa Symbol ? string(op.j) : string(Int(op.j))
    return "_" * i * j * _op_index_suffix(op)
end
_op_name_extra(op::SQA.Destroy) = _op_index_suffix(op)
_op_name_extra(op::SQA.Create) = "_dag" * _op_index_suffix(op)
function _op_name_extra(op::SQA.Pauli)
    return "_" * string(Int(op.axis)) * _op_index_suffix(op)
end
function _op_name_extra(op::SQA.Spin)
    return "_" * string(Int(op.axis)) * _op_index_suffix(op)
end

function _op_index_suffix(op::SQA.QSym)
    isdefined(op, :index) || return ""
    idx = op.index
    idx === SQA.NO_INDEX && return ""
    return "_" * string(idx.name)
end

function _avg_to_var_dict(eqs::AbstractMeanFieldEquations)
    iv = eqs.iv
    dict = Dict{SymbolicUtils.BasicSymbolic, Symbolics.Num}()
    dvs = Symbolics.Num[]
    for avg in eqs.states
        v = _make_time_dependent_var(_stable_avg_name(avg), iv)
        dict[avg] = v
        push!(dvs, v)
    end
    return dict, dvs
end

function _make_time_dependent_var(name::Symbol, iv::Symbolics.Num)
    v = first(@variables $name(iv))
    return v
end

function _collect_params!(set, x, dict, iv_uw)
    if x isa SymbolicUtils.BasicSymbolic
        SymbolicUtils.isconst(x) && return
        if !SymbolicUtils.iscall(x)
            if !haskey(dict, x) && x !== iv_uw && SymbolicUtils.symtype(x) <: Real
                push!(set, x)
            end
            return
        end
        op = SymbolicUtils.operation(x)
        args = SymbolicUtils.arguments(x)
        if length(args) == 1 && args[1] === iv_uw && op isa SymbolicUtils.BasicSymbolic
            return
        end
        for a in args
            _collect_params!(set, a, dict, iv_uw)
        end
    end
    return
end

"""
    to_system(eqs::MeanFieldEquations; name::Symbol)

Build a `ModelingToolkitBase.System` from the QC equation set. Substitutes
Averages with real-typed `u(t)` Num variables and passes `dvs`/`ps` explicitly.
"""
function to_system(
        eqs::NoiseMeanFieldEquations{O, H, Op, Jt, Jdt, R, E, S, Forward};
        name::Symbol
    ) where {O, H, Op, Jt, Jdt, R, E, S}
    return _to_system_sde(eqs, name, +1)
end

function to_system(
        eqs::NoiseMeanFieldEquations{O, H, Op, Jt, Jdt, R, E, S, Backward};
        name::Symbol
    ) where {O, H, Op, Jt, Jdt, R, E, S}
    return _to_system_sde(eqs, name, -1)
end

function _to_system_sde(eqs::NoiseMeanFieldEquations, name::Symbol, sign::Int)
    iv = eqs.iv
    iv_uw = SymbolicUtils.unwrap(iv)
    D = Symbolics.Differential(iv)
    dict, dvs = _avg_to_var_dict(eqs)
    conj_dict = _conj_substitution_dict(eqs, dict)
    # Single Brownian per system: `eqs.noise_equations` already aggregates the
    # per-jump noise drifts into one column. Multiple independent measurement
    # processes would need one Brownian per jump (TODO).
    w = first(MTK.@brownians _qc_dW)
    new_eqs = Vector{Symbolics.Equation}(undef, length(eqs.equations))
    ps_set = Set{SymbolicUtils.BasicSymbolic}()
    @inbounds for (i, eq) in enumerate(eqs.equations)
        rhs_conj = _substitute_conj_avgs(eq.rhs, conj_dict)
        rhs = _safe_substitute(rhs_conj, dict)
        noise_eq = eqs.noise_equations[i]
        noise_rhs_conj = _substitute_conj_avgs(noise_eq.rhs, conj_dict)
        noise_rhs = _safe_substitute(noise_rhs_conj, dict)
        # Substituted noise / drift trees can carry symtype `Any` (mixed
        # average products). SymbolicUtils refuses arithmetic between
        # mismatched symtypes, so build the product/sum nodes via `maketerm`,
        # which preserves structure without dispatching the worker buffer.
        w_uw = SymbolicUtils.unwrap(w)
        T = typeof(w_uw)
        signed_rhs = sign == 1 ? rhs :
            TermInterface.maketerm(T, *, Any[sign, rhs], nothing)
        noise_term = TermInterface.maketerm(T, *, Any[noise_rhs, w_uw], nothing)
        total_rhs = TermInterface.maketerm(T, +, Any[signed_rhs, noise_term], nothing)
        new_eqs[i] = D(dict[eq.lhs]) ~ total_rhs
        _collect_params!(ps_set, rhs, dict, iv_uw)
        _collect_params!(ps_set, noise_rhs, dict, iv_uw)
    end
    ps = [MTK.toparam(p) for p in ps_set]
    return MTK.System(new_eqs, iv, dvs, ps, [w]; name = name)
end

function to_system(eqs::MeanFieldEquations; name::Symbol)
    iv = eqs.iv
    iv_uw = SymbolicUtils.unwrap(iv)
    D = Symbolics.Differential(iv)
    dict, dvs = _avg_to_var_dict(eqs)
    conj_dict = _conj_substitution_dict(eqs, dict)
    new_eqs = Vector{Symbolics.Equation}(undef, length(eqs.equations))
    ps_set = Set{SymbolicUtils.BasicSymbolic}()
    @inbounds for (i, eq) in enumerate(eqs.equations)
        rhs_conj = _substitute_conj_avgs(eq.rhs, conj_dict)
        rhs = _safe_substitute(rhs_conj, dict)
        new_eqs[i] = D(dict[eq.lhs]) ~ rhs
        _collect_params!(ps_set, rhs, dict, iv_uw)
    end
    ps_old = collect(ps_set)
    ps = [MTK.toparam(p) for p in ps_old]
    return MTK.System(new_eqs, iv, dvs, ps; name = name)
end

# Build a substitution `⟨op†⟩ → conj(state_var(⟨op⟩))` for every state
# whose conjugate is *not* itself a state. This lets `to_system` codegen
# substitute the conjugate of a state without needing the conjugate to
# also be added to `eqs.states`. (See completion.jl::find_missing: by
# default, conjugates of states are excluded from missing-state scans.)
function _conj_substitution_dict(
        eqs::AbstractMeanFieldEquations,
        var_dict::AbstractDict
    )
    states = Set(eqs.states)
    conj_dict = Dict{SymbolicUtils.BasicSymbolic, Any}()
    for s in eqs.states
        cs = _avg_conj_for_codegen(s)
        cs === s && continue
        cs in states && continue
        haskey(var_dict, s) || continue
        conj_dict[cs] = conj(var_dict[s])
    end
    return conj_dict
end

function _avg_conj_for_codegen(x::SymbolicUtils.BasicSymbolic)
    SQA.is_average(x) || return x
    SymbolicUtils.iscall(x) || return x
    SymbolicUtils.operation(x) === SQA.sym_average || return x
    op = SQA.undo_average(x)
    return average(adjoint(op))
end

function _substitute_conj_avgs(x, conj_dict)
    isempty(conj_dict) && return x
    x isa SymbolicUtils.BasicSymbolic || return x
    if haskey(conj_dict, x)
        return conj_dict[x]
    end
    SymbolicUtils.iscall(x) || return x
    op = SymbolicUtils.operation(x)
    op === SQA.sym_average && return x
    args = SymbolicUtils.arguments(x)
    new_args = Any[_substitute_conj_avgs(a, conj_dict) for a in args]
    # `complex(re_sym, im_sym)` literals occasionally survive in the RHS
    # (e.g. from `im * H` constructing `Complex{Num}` coefficients). Calling
    # `op(new_args...)` with `op === complex` and symbolic args has no method;
    # rewrite to additive form `re + im_part * im` instead.
    # TODO fix https://github.com/JuliaSymbolics/SymbolicUtils.jl/pull/922
    op === complex && length(new_args) == 2 && return new_args[1] + new_args[2] * Symbolics.IM
    # Some interior operators (e.g. Σ-sums whose symtype isn't <: Number)
    # refuse `+` via SymbolicUtils's worker dispatch. Rebuild the node via
    # `maketerm` in those cases, preserving the original symtype/metadata.
    try
        return op(new_args...)
    catch err
        err isa MethodError && err.f === op || rethrow()
        return TermInterface.maketerm(
            typeof(x), op, new_args, TermInterface.metadata(x),
        )
    end
end

"""
    initial_values(eqs::AbstractMeanFieldEquations; defaults=Dict())

Return `Dict{Symbolics.Num, ComplexF64}` mapping each state's u(t) variable to
its initial value. Unspecified averages default to `zero(ComplexF64)`.
"""
function initial_values(
        eqs::AbstractMeanFieldEquations;
        defaults::AbstractDict = Dict()
    )
    dict, _ = _avg_to_var_dict(eqs)
    u0 = Dict{Symbolics.Num, ComplexF64}()
    for avg in eqs.states
        u0[dict[avg]] = ComplexF64(get(defaults, avg, 0))
    end
    return u0
end

"""
    initial_values(eqs::AbstractMeanFieldEquations, state; kwargs...)

For a set of symbolic equations `eqs` compute the initial state-average values
corresponding to the numeric quantum state `state` of the system. `state` can
be a `QuantumOpticsBase.StateVector` (Ket) or `Operator` (density matrix);
indexed/scaled equations are supported when `state` is a tensor product or
`LazyKet`.

Returns a `Vector{ComplexF64}` aligned with `eqs.states`.

See also: [`to_numeric`](@ref), [`numeric_average`](@ref).
"""
function initial_values(
        eqs::AbstractMeanFieldEquations, state;
        kwargs...
    )
    vals = ComplexF64[]
    for v in eqs.states
        push!(vals, ComplexF64(SQA.numeric_average(v, state; kwargs...)))
    end
    return vals
end

"""
    initial_values(eqs::AbstractMeanFieldEquations, u0::AbstractVector{<:Number})

Map a numeric initial-condition vector aligned with `eqs.states` to a
`Dict{Symbolics.Num, ComplexF64}` keyed by the MTK state variables produced
by `to_system`. Equivalent to building the dict via `unknowns(sys) .=> u0`.
"""
function initial_values(
        eqs::AbstractMeanFieldEquations, u0::AbstractVector{<:Number};
        kwargs...
    )
    length(u0) == length(eqs.states) || throw(DimensionMismatch(
        "initial value vector length $(length(u0)) does not match number of states $(length(eqs.states))",
    ))
    dict_var, _ = _avg_to_var_dict(eqs)
    out = Dict{Symbolics.Num, ComplexF64}()
    for (k, avg) in enumerate(eqs.states)
        out[dict_var[avg]] = ComplexF64(u0[k])
    end
    return out
end

"""
    parameter_map(eqs::AbstractMeanFieldEquations, pairs) -> Dict

Translate a user-facing `Pair`/`Dict` of parameter assignments into an
MTK-compatible parameter map. Accepts the original `IndexedVariable`
callable (e.g. `g` from `g(i) = IndexedVariable(:g, i)`) as a key and
matches it to the Symbolics-array parameter that `evaluate` synthesised
from the per-atom callable references. Scalar values are broadcast to
N-element vectors with N matching the array's concrete shape.

```julia
pmap = parameter_map(evaled, Dict(
    g(i)  => 0.1,           # broadcast to fill(0.1, N)
    Δa(i) => [0.0, 0.1],    # explicit per-atom values
    κ     => 1.0,           # scalar passes through
))
```
"""
function parameter_map(eqs::AbstractMeanFieldEquations, pairs)
    arr_by_name = Dict{Symbol, SymbolicUtils.BasicSymbolic}()
    scalar_by_name = Dict{Symbol, SymbolicUtils.BasicSymbolic}()
    for eq in eqs.equations
        _collect_named_params!(arr_by_name, scalar_by_name,
                               SymbolicUtils.unwrap(eq.rhs))
        _collect_named_params!(arr_by_name, scalar_by_name,
                               SymbolicUtils.unwrap(eq.lhs))
    end
    pmap = Dict{Any, Any}()
    for (k, v) in pairs
        ku = SymbolicUtils.unwrap(k)
        name = _param_name(ku)
        if name !== nothing && haskey(arr_by_name, name)
            arr = arr_by_name[name]
            shp = SymbolicUtils.shape(arr)
            n = (shp isa SymbolicUtils.SmallVec{UnitRange{Int}} && !isempty(shp)) ?
                length(shp[1]) : nothing
            if v isa AbstractArray
                pmap[arr] = v
            elseif n !== nothing
                pmap[arr] = fill(v, n)
            else
                pmap[arr] = v
            end
        elseif name !== nothing && haskey(scalar_by_name, name)
            pmap[scalar_by_name[name]] = v
        else
            pmap[k] = v
        end
    end
    return pmap
end

# Extract the user-visible name of a parameter expression:
# - bare Sym `:κ` → :κ
# - callable-Sym Term `g(i.sym)` (IndexedVariable shape) → :g
# - scalar after `unwrap(Num(...))` likewise
_param_name(x) = nothing
function _param_name(x::SymbolicUtils.BasicSymbolic)
    if !SymbolicUtils.iscall(x)
        return Base.nameof(x)
    end
    op = SymbolicUtils.operation(x)
    if op isa SymbolicUtils.BasicSymbolic && !SymbolicUtils.iscall(op)
        return Base.nameof(op)
    end
    return nothing
end

# Populate the two name → MTK-parameter dicts from a tree `x`.
# `arr_dict` collects array-typed params (from the `getindex` post-pass);
# `scalar_dict` collects bare-Sym scalar parameters (κ, Γ, …).
function _collect_named_params!(arr_dict, scalar_dict, x)
    x isa SymbolicUtils.BasicSymbolic || return
    if !SymbolicUtils.iscall(x)
        # Only Sym leaves have `nameof`; numeric Consts and the like don't.
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
end

"""
    get_solution(sol, avg_or_op, eqs)

Query an ODESolution `sol` for the trajectory of `avg_or_op`. Accepts either a
raw `Average` BasicSymbolic or a `QField` (which is averaged internally).
"""
function get_solution(
        sol, avg::SymbolicUtils.BasicSymbolic,
        eqs::AbstractMeanFieldEquations
    )
    dict, _ = _avg_to_var_dict(eqs)
    if haskey(dict, avg)
        var = dict[avg]
        return τ -> sol(τ; idxs = var)
    end
    # Conjugate fallback: ⟨op†⟩ is not stored explicitly because
    # `complete!` deduplicates conjugate pairs. Look up ⟨op⟩ instead and
    # return `conj` of its trajectory.
    conj_avg = _avg_conj_for_codegen(avg)
    if conj_avg !== avg && haskey(dict, conj_avg)
        var = dict[conj_avg]
        return τ -> conj.(sol(τ; idxs = var))
    end
    throw(KeyError(avg))
end
get_solution(sol, op::QField, eqs::AbstractMeanFieldEquations) =
    get_solution(sol, average(op), eqs)
