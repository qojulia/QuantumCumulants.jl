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
    ModelingToolkitBase.System(eqs::AbstractMeanFieldEquations; name::Symbol)
    ModelingToolkitBase.System(c::CorrelationFunction; name::Symbol)

Build a `ModelingToolkitBase.System` from the QC equation set. Substitutes
Averages with real-typed `u(t)` Num variables and passes `dvs`/`ps` explicitly.
For [`NoiseMeanFieldEquations`](@ref) the result is an SDE system whose
Brownian column is the aggregated per-jump noise drift. To compare against the
deterministic drift alone, pass `MeanFieldEquations(eqs)` instead.
"""
function MTK.System(
        eqs::NoiseMeanFieldEquations{O, H, Op, Jt, Jdt, R, E, S, Forward};
        name::Symbol,
    ) where {O, H, Op, Jt, Jdt, R, E, S}
    return _to_system_sde(eqs, name, +1)
end

function MTK.System(
        eqs::NoiseMeanFieldEquations{O, H, Op, Jt, Jdt, R, E, S, Backward};
        name::Symbol,
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
    w_uw = SymbolicUtils.unwrap(w)
    T = typeof(w_uw)
    new_eqs = Vector{Symbolics.Equation}(undef, length(eqs.equations))
    ps_set = Set{SymbolicUtils.BasicSymbolic}()
    merged = merge(conj_dict, dict)
    @inbounds for (i, eq) in enumerate(eqs.equations)
        rhs = _safe_substitute(eq.rhs, merged)
        noise_rhs = _safe_substitute(eqs.noise_equations[i].rhs, merged)
        # Substituted noise / drift trees can carry symtype `Any` (mixed
        # average products). SymbolicUtils refuses arithmetic between
        # mismatched symtypes, so build the product/sum nodes via `maketerm`,
        # which preserves structure without dispatching the worker buffer.
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

function MTK.System(eqs::MeanFieldEquations; name::Symbol)
    iv = eqs.iv
    iv_uw = SymbolicUtils.unwrap(iv)
    D = Symbolics.Differential(iv)
    dict, dvs = _avg_to_var_dict(eqs)
    conj_dict = _conj_substitution_dict(eqs, dict)
    new_eqs = Vector{Symbolics.Equation}(undef, length(eqs.equations))
    ps_set = Set{SymbolicUtils.BasicSymbolic}()
    merged = merge(conj_dict, dict)
    @inbounds for (i, eq) in enumerate(eqs.equations)
        rhs = _safe_substitute(eq.rhs, merged)
        new_eqs[i] = D(dict[eq.lhs]) ~ rhs
        _collect_params!(ps_set, rhs, dict, iv_uw)
    end
    ps_old = collect(ps_set)
    ps = [MTK.toparam(p) for p in ps_old]
    return MTK.System(new_eqs, iv, dvs, ps; name = name)
end

# Build a substitution `⟨op†⟩ → conj(state_var(⟨op⟩))` for every leaf on
# the RHS that is the conjugate of a tracked state. completion.jl emits
# one representative per conjugate-pair by default; this pass rewrites
# the missing partner before mtkcompile sees it.
#
# Direct adjoint (`_avg_conj_for_codegen(s)`) reverses operator order but
# does NOT re-sort cross-atom commuting ops on the same Hilbert subspace
# (per SQA's "undetermined free-index pairs stay put" invariant). So an
# RHS leaf like `⟨σ_k_1₂₁ * σ_k_2₂₂⟩` (constructed in atom-1-first order)
# does not literal-equal the adjoint of state `⟨σ_k_1₁₂ * σ_k_2₂₂⟩`,
# which lands as `⟨σ_k_2₂₂ * σ_k_1₂₁⟩` (atom-2-first). We walk RHS leaves
# and key the lookup on a permutation-canonical op signature so both
# orderings collapse to the same key.
function _conj_substitution_dict(
        eqs::AbstractMeanFieldEquations,
        var_dict::AbstractDict
    )
    state_set = Set(eqs.states)
    state_by_canon = Dict{String, Any}()
    for s in eqs.states
        haskey(var_dict, s) || continue
        v = SymbolicUtils.unwrap(var_dict[s])
        state_by_canon[_perm_canon_key(s)] = v
    end
    conj_dict = Dict{SymbolicUtils.BasicSymbolic, Any}()
    for eq in eqs.equations
        _collect_conj_subs!(conj_dict, eq.rhs, state_set, var_dict, state_by_canon)
    end
    if eqs isa NoiseMeanFieldEquations
        for eq in eqs.noise_equations
            _collect_conj_subs!(conj_dict, eq.rhs, state_set, var_dict, state_by_canon)
        end
    end
    return conj_dict
end

function _collect_conj_subs!(out, x, state_set, var_dict, state_by_canon)
    x isa SymbolicUtils.BasicSymbolic || return
    if SQA.is_average(x) && SymbolicUtils.iscall(x) &&
            SymbolicUtils.operation(x) === SQA.sym_average
        haskey(out, x) && return
        x in state_set && return       # state itself
        haskey(var_dict, x) && return  # literal match against the var dict
        # Direct match via permutation-canonical (NE-blind) key: leaves can
        # carry NE metadata or differ in operator order from the canonical
        # state form while still being the same physical average. Substitute
        # the state's var directly so the rhs lookup succeeds.
        k_direct = _perm_canon_key(x)
        if haskey(state_by_canon, k_direct)
            out[x] = state_by_canon[k_direct]
            return
        end
        # Conjugate match: same permutation-canonical key on the adjoint.
        cs = _avg_conj_for_codegen(x)
        ck = _perm_canon_key(cs)
        if haskey(state_by_canon, ck)
            v = state_by_canon[ck]
            # Build `conj(avg_var(t))` as a raw `SymbolicUtils.term` rather
            # than calling `Base.conj` directly. Julia's `conj` invokes
            # Symbolics' simplifier, which folds `conj(::SymReal)` to
            # identity and silently zeros every `⟨X⟩ - ⟨X†⟩` driving term
            # on the RHS. Using `term(...)` with `type=Number` skips the
            # simplifier so the symbolic `conj` node survives through
            # `_safe_substitute`, `mtkcompile`, and `build_function`.
            out[x] = SymbolicUtils.term(conj, v; type = Number)
        end
        return
    end
    SymbolicUtils.iscall(x) || return
    for a in SymbolicUtils.arguments(x)
        _collect_conj_subs!(out, a, state_set, var_dict, state_by_canon)
    end
    return
end

# Permutation-canonical key for an averaged operator. Stable-sort QTerm
# ops by their `acts_on` tuple then by `string(op)`; ops on the same
# Hilbert subspace (same acts_on) get a deterministic order based on
# their printed form, so the cross-atom adjoint reordering becomes
# invisible to dict-lookup.
function _perm_canon_key(avg::SymbolicUtils.BasicSymbolic)
    op = SQA.undo_average(avg)
    return _perm_canon_key(op)
end
function _perm_canon_key(op::SQA.QAdd)
    parts = String[]
    for (term, _) in op.arguments
        sorted = sort(
            collect(term.ops);
            by = o -> (Tuple(SQA.acts_on(o)), string(o)),
            alg = Base.Sort.MergeSort,
        )
        # NE-blind: matches `completion._canonical_dedup_key`. Two leaves that
        # differ only in NE metadata are the same physical state at codegen
        # time (the algebraic content has already been canonicalised
        # upstream via `_assume_distinct_atom_indices`).
        push!(parts, join(string.(sorted), '*'))
    end
    sort!(parts)
    idx_sig = sort(string.(op.indices))
    return string(join(parts, '+'), '|', idx_sig)
end

function _avg_conj_for_codegen(x::SymbolicUtils.BasicSymbolic)
    SQA.is_average(x) || return x
    SymbolicUtils.iscall(x) || return x
    SymbolicUtils.operation(x) === SQA.sym_average || return x
    op = SQA.undo_average(x)
    return average(adjoint(op))
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
    initial_values(eqs::AbstractMeanFieldEquations, state)

For a set of symbolic equations `eqs` compute the initial state-average values
corresponding to the numeric quantum state `state` of the system. `state` can
be a `QuantumOpticsBase.StateVector` (Ket) or `Operator` (density matrix);
indexed/scaled equations are supported when `state` is a tensor product or
`LazyKet`.

Returns a `Vector{ComplexF64}` aligned with `eqs.states`.

See also: [`to_numeric`](@ref), [`numeric_average`](@ref).
"""
function initial_values(eqs::AbstractMeanFieldEquations, state)
    vals = ComplexF64[]
    for v in eqs.states
        push!(vals, ComplexF64(SQA.numeric_average(v, state)))
    end
    return vals
end

"""
    initial_values(eqs::AbstractMeanFieldEquations, u0::AbstractVector{<:Number})

Map a numeric initial-condition vector aligned with `eqs.states` to a
`Dict{Symbolics.Num, ComplexF64}` keyed by the MTK state variables produced
by `System(eqs)`. Equivalent to building the dict via `unknowns(sys) .=> u0`.
"""
function initial_values(
        eqs::AbstractMeanFieldEquations, u0::AbstractVector{<:Number};
        kwargs...
    )
    length(u0) == length(eqs.states) || throw(
        DimensionMismatch(
            "initial value vector length $(length(u0)) does not match number of states $(length(eqs.states))",
        )
    )
    dict_var, _ = _avg_to_var_dict(eqs)
    out = Dict{Symbolics.Num, ComplexF64}()
    for (k, avg) in enumerate(eqs.states)
        out[dict_var[avg]] = ComplexF64(u0[k])
    end
    return out
end

"""
    parameter_map(sys::MTK.System, pairs)

Filter `pairs` (a `Pair` iterable or `Dict`) to entries whose key is a live
unknown or parameter of the compiled system `sys`. Use when a single
user-built parameter dict carries entries for several compiles of the same
physical model (e.g. a noisy / deterministic pair) and some compiles drop
parameters that the others keep. MTK v10 rejects superfluous entries with an
`Initial` parameter assertion; this strips them silently.

```julia
sys = mtkcompile(System(MeanFieldEquations(scaled_eqs); name=:sys))
dict = parameter_map(sys, merge(
    Dict(unknowns(sys) .=> u0),
    Dict(p .=> p0),
))
prob = ODEProblem(sys, dict, (0.0, T_end))
```
"""
function parameter_map(sys::MTK.System, pairs)
    live = Set{SymbolicUtils.BasicSymbolic}()
    for p in MTK.parameters(sys)
        push!(live, SymbolicUtils.unwrap(p))
    end
    for u in MTK.unknowns(sys)
        push!(live, SymbolicUtils.unwrap(u))
    end
    pmap = Dict{Any, Any}()
    for (k, v) in pairs
        if SymbolicUtils.unwrap(k) in live
            pmap[k] = v
        end
    end
    return pmap
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
        _collect_named_params!(
            arr_by_name, scalar_by_name,
            SymbolicUtils.unwrap(eq.rhs)
        )
        _collect_named_params!(
            arr_by_name, scalar_by_name,
            SymbolicUtils.unwrap(eq.lhs)
        )
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
    return
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
        return τ -> _eval_at(sol, var, τ)
    end
    # Conjugate fallback: ⟨op†⟩ is not stored explicitly because
    # `complete!` deduplicates conjugate pairs. Look up ⟨op⟩ instead and
    # return `conj` of its trajectory.
    conj_avg = _avg_conj_for_codegen(avg)
    if conj_avg !== avg && haskey(dict, conj_avg)
        var = dict[conj_avg]
        return τ -> conj.(_eval_at(sol, var, τ))
    end
    throw(KeyError(avg))
end

# `sol(τ; idxs=var)` for vector `τ` returns a plain Vector on ODE solutions
# but a `RecursiveArrayTools.DiffEqArray` on SDE solutions, which downstream
# `real.(...)` / broadcasting cannot consume. `Array(...)` materialises both
# uniformly in one batched interpolator pass.
_eval_at(sol, var, τ::AbstractVector) = Array(sol(τ; idxs = var))
_eval_at(sol, var, τ) = sol(τ; idxs = var)
get_solution(sol, op::QField, eqs::AbstractMeanFieldEquations) =
    get_solution(sol, average(op), eqs)
