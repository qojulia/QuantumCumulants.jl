"""
    evaluate(eqs::MeanFieldEquations; limits=nothing, h::Vector{Int}=Int[], kwargs...)
    evaluate(c::CorrelationFunction; limits=nothing, h::Vector{Int}=Int[], kwargs...)

Materialise symbolic indexed mean-field equations into concrete-size systems.

`limits` accepts a `Pair` (`N => 3`), tuple of pairs, or `Dict{BasicSymbolic,Int}`.
For every `Index` whose `range` matches a limits key:

  1. Sums `Σ_i^N f(i)` (carried as `QAdd.indices` metadata) are unrolled to
     `f(e_s_1) + f(e_s_2) + ... + f(e_s_n)`, where `e_s_k` is a fresh per-value
     `Index` shared across all enumerations of "position k in space s", so
     `i_k` and `j_k` from different source names collapse to the same physical
     atom. NE-violated terms drop out of `_canonicalize!`.

  2. Free LHS indices enumerate over their range, producing one equation per
     concrete value.

The result is a system whose `QAdd`s have empty `.indices`, ready for `System`
and `mtkcompile`.

`h` further restricts unrolling to Hilbert subspaces whose `space_index` is in
the vector (the 1-based tensor-product position). The empty default unrolls
every index covered by `limits`. Indices on unselected subspaces keep their
symbolic names and sum-scope, even if `limits` covers their range. Useful for
hybrid systems: unroll a filter cavity array via `evaluate(...; limits, h=[2])`
while leaving the atom ensemble symbolic for `scale`.

Implementation: this is a thin layer over `SQA.change_index`. All algebraic
normalisation (commuting-op sort, projector squashing, ne-violated drop,
diagonal splitting) is delegated to SQA primitives.
"""
function evaluate(
        eqs::MeanFieldEquations;
        limits = nothing,
        h::Vector{Int} = Int[],
        kwargs...,
    )
    sub_dict = _build_limits_dict(limits)
    isempty(sub_dict) && return _copy(eqs)
    h_set = Set{Int}(h)
    return _evaluate_unroll(eqs, sub_dict, h_set)
end

function evaluate(c::CorrelationFunction; limits = nothing, h::Vector{Int} = Int[], kwargs...)
    new_eqs = evaluate(c.eqs; limits = limits, h = h, kwargs...)
    return CorrelationFunction(
        c.op1, c.op2, c.op2_anc, c.aon_anc,
        new_eqs, c.eqs0, c.τ, c.steady_state,
    )
end

# ---------------------------------------------------------------------------
# Driver

function _evaluate_unroll(eqs::MeanFieldEquations, sub_dict, h_set::Set{Int})
    coeff_sub = Dict{Any, Any}(k => v for (k, v) in sub_dict)
    canon = _build_canonical_indices(eqs)
    # Normalise scope: lift SumIndices metadata from average leaves onto
    # their enclosing Mul parents so coefficient siblings (e.g. IndexedVariable
    # `g(i)`) are substituted together with the average leaf during unrolling.
    # See docs/superpowers/specs/2026-05-19-bound-index-scope-lift-design.md.
    lifted_rhs = [_lift_sum_scope(eq.rhs) for eq in eqs.equations]
    lifted_lhs = [_lift_sum_scope(eq.lhs) for eq in eqs.equations]
    lifted_op_eq_lhs = [_lift_sum_scope(opeq.lhs) for opeq in eqs.operator_equations]
    lifted_op_eq_rhs = [_lift_sum_scope(opeq.rhs) for opeq in eqs.operator_equations]
    new_eqs = Symbolics.Equation[]
    new_states = SymbolicUtils.BasicSymbolic[]
    new_ops = QAdd[]
    new_op_eqs = Symbolics.Equation[]
    seen = Set{Any}()

    for k in eachindex(eqs.equations)
        op_k = eqs.operators[k]
        free = _free_lhs_indices(op_k, sub_dict, h_set)
        ranges = [sub_dict[SymbolicUtils.unwrap(idx.range)] for idx in free]
        iter = Iterators.product([1:r for r in ranges]...)
        for tup in iter
            idx_sub = Dict{SQA.Index, SQA.Index}(
                zip(free, (_fresh_index(b, v, canon) for (b, v) in zip(free, tup)))
            )
            # Matching Symbolics-side substitution so `IndexedVariable`
            # references like `g(i.sym)` get renamed in lock-step with the
            # operator-side `change_index` calls. `_materialise` applies this
            # at every leaf Sym during its walk; pre-applying it via
            # `Symbolics.substitute` is too eager (it rebuilds Add nodes
            # whose shape-promote step fails on QC's indexed RHS shapes).
            sym_sub = isempty(idx_sub) ? _EMPTY_SYM_SUB :
                Dict{Any, Any}(
                    SymbolicUtils.unwrap(b.sym) => SymbolicUtils.unwrap(t.sym)
                    for (b, t) in idx_sub
                )
            new_lhs = _materialise(lifted_lhs[k], idx_sub, sub_dict, canon, sym_sub, h_set)
            # If the substitution made a same-name pair contradictory under
            # NE, `change_index` (via SQA's contradiction guard) drops the
            # term entirely. Skip the iteration in that case: an LHS of `0`
            # has no state semantics and would trip `_stable_avg_name`.
            _is_materialised_zero(new_lhs) && continue
            new_rhs = _materialise(lifted_rhs[k], idx_sub, sub_dict, canon, sym_sub, h_set)
            new_op = _materialise_qfield(op_k, idx_sub, sub_dict, canon, h_set)
            new_op_eq = _materialise(lifted_op_eq_lhs[k], idx_sub, sub_dict, canon, sym_sub, h_set) ~
                _materialise(lifted_op_eq_rhs[k], idx_sub, sub_dict, canon, sym_sub, h_set)
            key = _dedup_key(new_lhs)
            (key in seen) && continue
            push!(seen, key)
            push!(seen, _dedup_key_conj(new_lhs))
            final_rhs = _safe_substitute(new_rhs, coeff_sub)
            push!(new_states, new_lhs)
            push!(new_eqs, new_lhs ~ final_rhs)
            push!(new_ops, _as_qadd(new_op))
            push!(new_op_eqs, new_op_eq)
        end
    end

    # Canonicalise every average leaf in every RHS so the literal symbol
    # matches the state's literal symbol. `_canonical_key` folds
    # alpha-equivalent variants (e.g. `⟨σ_{j_2}^{11} σ_{j_1}^{21}⟩` and
    # `⟨σ_{j_1}^{11} σ_{j_2}^{21}⟩` for permutation-symmetric atoms) so
    # `find_missing` reports closed, but MTK's codegen uses literal equality
    # against `unknowns(sys)`. Without this pass the RHS-side leaf has a
    # different symbol than the state-side leaf and `_stable_avg_name`
    # treats the RHS leaf as a callable `AvgFunc` at runtime, which crashes.
    canon_post = _canon_from_states(canon, new_states)
    state_map = Dict{SQA.QAdd, SymbolicUtils.BasicSymbolic}()
    for s in new_states
        key = _canonical_key(s, canon_post)
        key isa SQA.QAdd && (state_map[key] = s)
    end
    for k in eachindex(new_eqs)
        new_eqs[k] = new_eqs[k].lhs ~
            _canonicalise_avg_leaves(new_eqs[k].rhs, canon_post, state_map)
    end

    # Convert surviving callable-Sym references (e.g. `g(i_k.sym)` from
    # `IndexedVariable`) into Symbolics-array `getindex` terms with concrete
    # shape. This makes the parameter MTK-native: a single array `g` of
    # size N rather than a callable Sym whose `FnType{Tuple{Int}}` signature
    # trips MTK's `_has_delays` precheck (which constructs `g(::Real)`).
    arr_sub = _build_callable_to_array_sub(new_eqs, new_states)
    if !isempty(arr_sub)
        for k in eachindex(new_eqs)
            new_eqs[k] = _safe_substitute(new_eqs[k].lhs, arr_sub) ~
                _safe_substitute(new_eqs[k].rhs, arr_sub)
        end
        for k in eachindex(new_states)
            new_states[k] = _safe_substitute(new_states[k], arr_sub)
        end
    end

    out = _copy(eqs)
    empty!(out.equations);          append!(out.equations, new_eqs)
    empty!(out.operator_equations); append!(out.operator_equations, new_op_eqs)
    empty!(out.states);              append!(out.states, new_states)
    empty!(out.operators);           append!(out.operators, new_ops)
    return out
end

# Collect every `f(arg_sym)` callable-Sym Term across the given equations
# whose `f` is a callable Sym (FnType) and whose args are leaf Int-typed
# Syms (the fresh per-atom indices minted by `_fresh_index`). Returns a
# substitution dict mapping each such Term to a `getindex(arr_f, slot)`
# Term, where `arr_f` is a freshly minted Symbolics-array parameter of
# the correct concrete size.
function _build_callable_to_array_sub(
        eqs::Vector{Symbolics.Equation},
        states::Vector{SymbolicUtils.BasicSymbolic},
    )
    uses = Dict{Symbol, Dict{Vector{Int}, Vector{SymbolicUtils.BasicSymbolic}}}()
    for eq in eqs
        _collect_callable_uses!(uses, SymbolicUtils.unwrap(eq.lhs))
        _collect_callable_uses!(uses, SymbolicUtils.unwrap(eq.rhs))
    end
    for s in states
        _collect_callable_uses!(uses, s)
    end
    isempty(uses) && return Dict{Any, Any}()
    sub = Dict{Any, Any}()
    for (name, slot_map) in uses
        n_args = length(first(keys(slot_map)))
        max_per_dim = [maximum(k[d] for k in keys(slot_map)) for d in 1:n_args]
        arr = SymbolicUtils.Sym{SymbolicUtils.SymReal}(
            name;
            type = Real,
            shape = SymbolicUtils.SmallVec{UnitRange{Int}}([1:m for m in max_per_dim]),
        )
        for (slots, terms) in slot_map
            indexed = TermInterface.maketerm(
                typeof(arr), getindex, Any[arr, slots...], nothing,
            )
            for t in terms
                sub[t] = indexed
            end
        end
    end
    return sub
end

# FnType callable with leaf-Sym arguments whose names encode slots
# (`g(i_3)`, `Γ(i_1, i_2)`, …): bucket by name + slot tuple so the emitter
# can mint a Symbolics array of the right rank. Non-matching shapes recurse.
function _collect_callable_uses!(uses, x)
    x isa SymbolicUtils.BasicSymbolic || return
    SymbolicUtils.iscall(x) || return
    op = SymbolicUtils.operation(x)
    args = SymbolicUtils.arguments(x)
    if op isa SymbolicUtils.BasicSymbolic &&
            !SymbolicUtils.iscall(op) &&
            _is_fntype(SymbolicUtils.symtype(op)) &&
            all(a -> a isa SymbolicUtils.BasicSymbolic && !SymbolicUtils.iscall(a), args)
        slots = [_parse_slot(Base.nameof(a)) for a in args]
        if all(!isnothing, slots)
            slot_key = Int[s for s in slots]
            name = Base.nameof(op)
            bucket = get!(uses, name) do
                Dict{Vector{Int}, Vector{SymbolicUtils.BasicSymbolic}}()
            end
            push!(get!(bucket, slot_key, SymbolicUtils.BasicSymbolic[]), x)
            return
        end
    end
    for a in args
        _collect_callable_uses!(uses, a)
    end
    return
end

_is_fntype(::Type) = false
_is_fntype(::Type{<:SymbolicUtils.FnType}) = true

# Extract the integer slot from a fresh-index symbol name like `:i_3`.
# Returns `nothing` if the name doesn't fit the `<base>_<int>` shape.
function _parse_slot(name::Symbol)
    s = String(name)
    idx = findlast('_', s)
    idx === nothing && return nothing
    tail = SubString(s, idx + 1)
    return tryparse(Int, tail)
end

# ---------------------------------------------------------------------------
# limits parsing

function _build_limits_dict(limits)
    sub = Dict{Any, Int}()
    limits === nothing && return sub
    if limits isa Pair
        _push_limit!(sub, limits)
    elseif limits isa Tuple
        for p in limits
            _push_limit!(sub, p)
        end
    elseif limits isa AbstractDict
        for p in limits
            _push_limit!(sub, p)
        end
    else
        throw(ArgumentError("`limits` must be a Pair, Tuple of Pairs, or Dict"))
    end
    return sub
end

function _push_limit!(sub::AbstractDict, p::Pair)
    sub[SymbolicUtils.unwrap(first(p))] = Int(last(p))
    return sub
end

# ---------------------------------------------------------------------------
# Free-index discovery (only those whose range is in `sub_dict`)

_free_lhs_indices(::Number, _, _) = SQA.Index[]
function _free_lhs_indices(op::SQA.QSym, sub_dict, h_set::Set{Int})
    idx = op.index
    SQA.has_index(idx) || return SQA.Index[]
    SymbolicUtils.unwrap(idx.range) in keys(sub_dict) || return SQA.Index[]
    _in_h(h_set, idx.space_index) || return SQA.Index[]
    return SQA.Index[idx]
end
function _free_lhs_indices(op::QAdd, sub_dict, h_set::Set{Int})
    out = SQA.Index[]
    bound = Set(op.indices)
    for (term, _) in op.arguments
        for o in term.ops
            idx = o.index
            SQA.has_index(idx) || continue
            idx in bound && continue
            SymbolicUtils.unwrap(idx.range) in keys(sub_dict) || continue
            _in_h(h_set, idx.space_index) || continue
            idx in out || push!(out, idx)
        end
    end
    return out
end

# Mint a per-(space, value) index. Name derives from the user's first-declared
# index for the source Hilbert subspace (SQA naming policy: algebra never
# invents indices outside the user's vocabulary; we trace back to one of the
# user's declared names and just suffix the position). Two free indices `i`
# and `j` on the same atom-space, after evaluate substitutes both to position
# 1, end up as the *same* `Index(:i_1, …)` (assuming `:i` is the user's
# alphabetically-first declared atom-space index), so the resulting averages
# dedup via op equality.

function _fresh_index(b::SQA.Index, k::Int, canon)
    return _canonical_base(b, canon)(k)
end

function _canonical_base(b::SQA.Index, canon)
    list = get(canon, b.space_index, SQA.Index[])
    isempty(list) && return b
    return list[1]
end

# ---------------------------------------------------------------------------
# Scope lift: SQA's `average(::QAdd)` stamps SumIndices metadata on each
# average leaf but flattens `c * avg` so the coefficient escapes the scope.
# This pass copies the metadata from each scoped average leaf onto its
# enclosing Mul parent so the whole product is recognised as scoped during
# unrolling. Walk is post-order; metadata stays on leaves as well (strictly
# additive, no consumer is affected).

function _lift_sum_scope(x)
    x isa SymbolicUtils.BasicSymbolic || return x
    SymbolicUtils.iscall(x) || return x
    args = SymbolicUtils.arguments(x)
    new_args = Any[_lift_sum_scope(a) for a in args]
    op = SymbolicUtils.operation(x)
    # Compare by object identity, not isequal: setmetadata returns a new
    # object that isequal-matches the original (metadata is invisible to
    # `isequal`). Identity-comparison forces the parent to rebuild so the
    # newly-tagged child is actually carried up the tree.
    needs_rebuild = any(((a, b),) -> a !== b, zip(args, new_args))
    new_x = needs_rebuild ?
        TermInterface.maketerm(typeof(x), op, new_args, TermInterface.metadata(x)) :
        x
    if op === (*)
        for a in new_args
            a isa SymbolicUtils.BasicSymbolic || continue
            SymbolicUtils.hasmetadata(a, SQA.SumIndices) || continue
            idxs = SymbolicUtils.getmetadata(a, SQA.SumIndices)
            new_x = SymbolicUtils.setmetadata(new_x, SQA.SumIndices, idxs)
            if SymbolicUtils.hasmetadata(a, SQA.SumNonEqual)
                ne = SymbolicUtils.getmetadata(a, SQA.SumNonEqual)
                new_x = SymbolicUtils.setmetadata(new_x, SQA.SumNonEqual, ne)
            end
            break
        end
    end
    return new_x
end

# ---------------------------------------------------------------------------
# Tree walker: substitute free indices, unroll sum-metadata averages.

const _EMPTY_SYM_SUB = Dict{Any, Any}()

function _materialise(x, idx_sub, sub_dict, canon, sym_sub, h_set::Set{Int})
    x isa SymbolicUtils.BasicSymbolic || return x
    if _is_leaf_average(x)
        return _materialise_avg(x, idx_sub, sub_dict, canon, h_set)
    end
    if !SymbolicUtils.iscall(x)
        # Leaf Sym: substitute Symbolics-side index references (e.g. the
        # `i.sym` arg of `IndexedVariable(:g, i)`) so they track the
        # operator-side `change_index` renames.
        return get(sym_sub, x, x)
    end
    if SymbolicUtils.operation(x) === (*) &&
            SymbolicUtils.hasmetadata(x, SQA.SumIndices)
        bound = SymbolicUtils.getmetadata(x, SQA.SumIndices)
        bound isa Vector{SQA.Index} && !isempty(bound) &&
            return _materialise_scoped(x, bound, idx_sub, sub_dict, canon, sym_sub, h_set)
    end
    op = SymbolicUtils.operation(x)
    args = SymbolicUtils.arguments(x)
    new_args = Any[_materialise(a, idx_sub, sub_dict, canon, sym_sub, h_set) for a in args]
    # Identity comparison: `isequal` ignores BasicSymbolic metadata and may
    # collapse averages whose inner QFields differ, causing the rebuilt
    # parent to silently reuse the stale child. See `_lift_sum_scope` for
    # the same trap.
    any(((a, b),) -> a !== b, zip(args, new_args)) || return x
    op === complex && length(new_args) == 2 &&
        return new_args[1] + new_args[2] * Symbolics.IM
    try
        return op(new_args...)
    catch err
        err isa MethodError || err isa ArgumentError || rethrow()
        return TermInterface.maketerm(
            typeof(x), op, new_args, TermInterface.metadata(x),
        )
    end
end

# Unroll a scoped `Mul` node: for each bound index whose range is in `sub_dict`
# and is not already covered by `idx_sub` (LHS-free substitution), enumerate
# 1:n and recurse into the stripped subtree with the per-enumeration
# `local_sub` / `local_sym` extensions. `_materialise` propagates `local_sub`
# to operator leaves via `change_index` and `local_sym` to Symbolics leaves
# (e.g. the `i.sym` arg inside `IndexedVariable(:g, i)`).
function _materialise_scoped(x, bound::Vector{SQA.Index}, idx_sub, sub_dict, canon, sym_sub, h_set::Set{Int})
    substituted = Set{SQA.Index}()
    for (from, to) in idx_sub
        push!(substituted, from); push!(substituted, to)
    end
    # SQA's `average(::QAdd)` stamps SumIndices on every term of a
    # `.indices != []` QAdd, including terms that don't actually reference
    # the bound index (e.g. `κ⟨a'a⟩` inside a sum-over-i Hamiltonian).
    # Filter `to_unroll` to bound indices whose name actually appears in
    # the subtree, otherwise enumeration produces N identical copies.
    # `h_set` skips indices on unselected subspaces so their sum stays
    # symbolic, matching `_materialise_qfield`'s gate.
    to_unroll = SQA.Index[
        b for b in bound
            if SymbolicUtils.unwrap(b.range) in keys(sub_dict) &&
            _in_h(h_set, b.space_index) &&
            !(b in substituted) &&
            _references_index(x, b)
    ]
    stripped = SymbolicUtils.setmetadata(x, SQA.SumIndices, SQA.Index[])
    if isempty(to_unroll)
        return _materialise(stripped, idx_sub, sub_dict, canon, sym_sub, h_set)
    end
    ranges = [sub_dict[SymbolicUtils.unwrap(b.range)] for b in to_unroll]
    terms = Any[]
    for tup in Iterators.product([1:r for r in ranges]...)
        local_sub = copy(idx_sub)
        local_sym = copy(sym_sub)
        for (b, v) in zip(to_unroll, tup)
            fresh = _fresh_index(b, v, canon)
            local_sub[b] = fresh
            local_sym[SymbolicUtils.unwrap(b.sym)] = SymbolicUtils.unwrap(fresh.sym)
        end
        push!(terms, _materialise(stripped, local_sub, sub_dict, canon, local_sym, h_set))
    end
    isempty(terms) && return 0
    length(terms) == 1 && return terms[1]
    # SymbolicUtils' `+(::T, ::T)` path enforces a numeric-symtype guard
    # that rejects Any-typed children (`g(i.sym::Int)` gives symtype Any).
    # Build the Add directly via `maketerm`.
    proto_idx = findfirst(t -> t isa SymbolicUtils.BasicSymbolic, terms)
    proto_idx === nothing && return sum(terms)
    return TermInterface.maketerm(typeof(terms[proto_idx]), +, terms, nothing)
end

# True if any leaf in `x` (a BasicSymbolic) actually references `b` -
# either via the Symbolics variable `b.sym` (e.g. an `IndexedVariable`
# `g(b.sym)`) or via a QField inside an average leaf whose operator carries
# `b` as its index. Deliberately does NOT consult SumIndices metadata,
# which SQA's `average(::QAdd)` stamps uniformly on every term of a
# scoped QAdd (including terms that don't actually depend on the bound
# index, e.g. `κ⟨a'a⟩` inside a sum-over-i Hamiltonian).
function _references_index(x, b::SQA.Index)
    x isa SymbolicUtils.BasicSymbolic || return false
    if _is_leaf_average(x)
        op = SQA.undo_average(x)
        return _qfield_references_index(op, b)
    end
    if SymbolicUtils.iscall(x)
        return any(a -> _references_index(a, b), SymbolicUtils.arguments(x))
    end
    return isequal(x, SymbolicUtils.unwrap(b.sym))
end

function _qfield_references_index(op::SQA.QAdd, b::SQA.Index)
    for (term, _) in op.arguments, o in term.ops
        SQA.has_index(o.index) || continue
        o.index === b && return true
    end
    return false
end
_qfield_references_index(op::SQA.QSym, b::SQA.Index) =
    SQA.has_index(op.index) && op.index === b
_qfield_references_index(_, ::SQA.Index) = false

function _materialise_avg(avg::SymbolicUtils.BasicSymbolic, idx_sub, sub_dict, canon, h_set::Set{Int})
    op = SQA.undo_average(avg)
    new_op = _materialise_qfield(op, idx_sub, sub_dict, canon, h_set)
    new_op isa SQA.QField && return average(new_op)
    return new_op  # already a numeric sum-of-averages from unrolling
end

# Materialise a `QField` (or `Number`). Two passes: first apply free-index
# substitution via SQA.change_index, then unroll any sum-metadata bound indices
# whose range matches `sub_dict` and whose Hilbert subspace is in `h_set`.

_materialise_qfield(x::Number, _, _, _, _) = x
function _materialise_qfield(op::SQA.QSym, idx_sub, _, _, _)
    haskey(idx_sub, op.index) || return op
    return SQA.change_index(op, op.index, idx_sub[op.index])
end
function _materialise_qfield(op::QAdd, idx_sub, sub_dict, canon, h_set::Set{Int})
    after_free = _apply_free_sub(op, idx_sub)
    # Two sources of "bound indices to unroll" in this leaf:
    #  (a) explicit sum-scope `.indices` metadata that actually appears in
    #      some operator inside `term.ops`. Bound indices left in
    #      `.indices` whose name no longer appears in any op are spurious
    #      metadata residue from commutator/factorisation (cumulant
    #      expansion of a sum over `i1,i2` leaves a `⟨S_i₂₂⟩` factor
    #      that inherits `.indices = [i1,i2]` even though `i₂₂` references
    #      the LHS-free index `i`, not `i1` or `i2`). Treat those as no-op
    #      and strip without unrolling, matching `scaling.jl`'s behaviour.
    #  (b) op-level indices referenced by an operator inside `term.ops` whose
    #      range is in `sub_dict` but is not bound by the surrounding LHS
    #      substitution. Cumulant expansion can break apart a `Σ_{i2}` and
    #      leave a bare `⟨S_i2⟩` leaf whose `.indices` is empty; semantically
    #      the implicit `Σ_{i2}` still needs to be unrolled by evaluate.
    # Indices we have ALREADY substituted (LHS-free path) must not be
    # re-unrolled. After `_apply_free_sub`, the operator carries the *target*
    # index from `idx_sub`, not the source, so we exclude both.
    substituted = Set{SQA.Index}()
    for (from, to) in idx_sub
        push!(substituted, from)
        push!(substituted, to)
    end
    op_indices = SQA.Index[]
    for (term, _) in after_free.arguments, o in term.ops
        idx = o.index
        SQA.has_index(idx) || continue
        idx in op_indices || push!(op_indices, idx)
    end
    # Strip stale `.indices` entries on two grounds:
    #   (a) The index no longer appears in any operator inside the QAdd
    #       (residue from SQA's commutator / cumulant-expansion pipeline;
    #       `scaling.jl` treats it the same way).
    #   (b) The index IS an LHS-free substitution target (`idx_sub` value).
    #       `_apply_free_sub` calls `change_index(qadd, from, to)`, which
    #       renames the index in BOTH operators AND `.indices`. The post-
    #       substitution QAdd then looks like `Σ_{j_k} σ_{j_k}₂₂`, but `j_k`
    #       is a CONCRETE atom (from the LHS enumeration), not a real sum.
    #       Drop the spurious sum scope.
    # Indices we KEEP: those whose range is outside `sub_dict` (different
    # bound subspace) AND which appear in some operator AND which are not
    # substituted.
    # "Targeted" = range covered by sub_dict AND on a subspace selected by
    # h_set. Indices on unselected subspaces stay in `.indices` as symbolic
    # sum scope, matching `_materialise_scoped`'s gate.
    _targeted(idx) = SymbolicUtils.unwrap(idx.range) in keys(sub_dict) &&
        _in_h(h_set, idx.space_index)
    stripped_indices = SQA.Index[
        b for b in after_free.indices
            if (!_targeted(b) || b in op_indices) && !(b in substituted)
    ]
    if length(stripped_indices) != length(after_free.indices)
        after_free = QAdd(after_free.arguments, stripped_indices)
    end
    bound = SQA.Index[]
    for b in stripped_indices
        _targeted(b) || continue
        b in substituted && continue
        b in bound || push!(bound, b)
    end
    for idx in op_indices
        idx in substituted && continue
        _targeted(idx) || continue
        idx in bound || push!(bound, idx)
    end
    isempty(bound) && return after_free
    return _unroll_bound_indices(after_free, bound, sub_dict, canon)
end

# Apply each free-index substitution via SQA.change_index. The QAdd's own
# .indices field may be non-empty (sum scope is preserved across this step);
# SQA.change_index handles diag-splitting correctly when needed.
function _apply_free_sub(op::QAdd, idx_sub::Dict{SQA.Index, SQA.Index})
    isempty(idx_sub) && return op
    result = op
    for (from, to) in idx_sub
        result = SQA.change_index(result, from, to)
    end
    return result
end

# Unroll one bound index at a time. Each unroll step:
#  1. Strips that index from .indices (drops sum scope).
#  2. For k = 1..n, calls SQA.change_index(bare, bound, fresh_k). Because bare
#     has empty .indices, change_index dispatches through _canonicalize!,
#     which handles projector squashing, commuting-op sort, and NE-violated
#     term drop natively.
#  3. Sums the n results.
# After unrolling all bound indices, the result is a sum of `average(qadd)` /
# Numeric scalar expressions, ready to be substituted by `coeff_sub`.
function _unroll_bound_indices(op::QAdd, bound::Vector{SQA.Index}, sub_dict, canon)
    cur::Any = op
    for b in bound
        n = sub_dict[SymbolicUtils.unwrap(b.range)]
        cur = _unroll_one(cur, b, n, canon)
    end
    return cur
end

_unroll_one(x::Number, _, _, _) = x
function _unroll_one(qadd::QAdd, b::SQA.Index, n::Int, canon)
    bare = QAdd(qadd.arguments, SQA.Index[idx for idx in qadd.indices if idx != b])
    terms = Any[]
    for k in 1:n
        push!(terms, average(SQA.change_index(bare, b, _fresh_index(b, k, canon))))
    end
    isempty(terms) && return 0
    return reduce(+, terms)
end
function _unroll_one(x::SymbolicUtils.BasicSymbolic, b::SQA.Index, n::Int, canon)
    if _is_leaf_average(x)
        op = SQA.undo_average(x)
        op isa QAdd || return x
        b in op.indices || return x
        return _unroll_one(op, b, n, canon)
    end
    SymbolicUtils.iscall(x) || return x
    op = SymbolicUtils.operation(x)
    args = SymbolicUtils.arguments(x)
    new_args = Any[_unroll_one(a, b, n, canon) for a in args]
    any(((a, b2),) -> !isequal(a, b2), zip(args, new_args)) || return x
    op === complex && length(new_args) == 2 &&
        return new_args[1] + new_args[2] * Symbolics.IM
    try
        return op(new_args...)
    catch err
        err isa MethodError || err isa ArgumentError || rethrow()
        return TermInterface.maketerm(
            typeof(x), op, new_args, TermInterface.metadata(x),
        )
    end
end

# ---------------------------------------------------------------------------
# Dedup key: the underlying operator. Because `_fresh_index` produces a
# single canonical Index per (space, value), two states for the same
# physical atoms have op-equal QAdds. Commuting ops with different concrete
# indices may still survive in non-canonical order (SQA's `Undetermined`
# regime); known limitation, see DESIGN.md "canonicalise_undetermined".

_dedup_key(avg::SymbolicUtils.BasicSymbolic) = SQA.undo_average(avg)

# Build a canonical-index registry that includes both the original user
# vocabulary (from the input `canon`) and every free index appearing in
# `new_states` (the materialised atom names like `j_1, j_2`). This ensures
# `_canonical_key` has enough canonical entries to rename all encountered
# positions in 2-atom and higher correlations.
function _canon_from_states(canon, states::Vector{SymbolicUtils.BasicSymbolic})
    out = Dict{Int, Vector{SQA.Index}}()
    for (sp, v) in canon
        out[sp] = copy(v)
    end
    for s in states
        op = SQA.undo_average(s)
        op isa SQA.QAdd || continue
        for (term, _) in op.arguments, o in term.ops
            SQA.has_index(o.index) || continue
            v = get!(out, o.index.space_index, SQA.Index[])
            o.index in v || push!(v, o.index)
        end
    end
    for v in values(out)
        sort!(v, by = idx -> idx.name)
    end
    return out
end

# Walk `x` and replace every leaf `average(...)` whose canonical key is in
# `state_map` with the corresponding state symbol. Leaves whose canonical
# form is not a state are left untouched.
function _canonicalise_avg_leaves(x, canon, state_map::Dict{SQA.QAdd, SymbolicUtils.BasicSymbolic})
    x isa SymbolicUtils.BasicSymbolic || return x
    if _is_leaf_average(x)
        key = _canonical_key(x, canon)
        key isa SQA.QAdd || return x
        rep = get(state_map, key, nothing)
        rep === nothing && return x
        return rep
    end
    SymbolicUtils.iscall(x) || return x
    op = SymbolicUtils.operation(x)
    args = SymbolicUtils.arguments(x)
    new_args = Any[_canonicalise_avg_leaves(a, canon, state_map) for a in args]
    all(((a, b),) -> a === b, zip(args, new_args)) && return x
    op === complex && length(new_args) == 2 &&
        return new_args[1] + new_args[2] * Symbolics.IM
    try
        return op(new_args...)
    catch err
        err isa MethodError || err isa ArgumentError || rethrow()
        return TermInterface.maketerm(typeof(x), op, new_args, TermInterface.metadata(x))
    end
end

_is_materialised_zero(x::Number) = iszero(x)
function _is_materialised_zero(x::SymbolicUtils.BasicSymbolic)
    SymbolicUtils.isconst(x) && return iszero(x.val)
    return false
end
_is_materialised_zero(_) = false
_dedup_key_conj(avg::SymbolicUtils.BasicSymbolic) = adjoint(SQA.undo_average(avg))

# Promote `x` (a QSym, QAdd, or Number) to a QAdd. SQA's `*(::QSym, ::Int)`
# does the wrapping via the canonicalisation pipeline; for Numbers we wrap
# into a zero-op QAdd (the constant 1 carrier).
_as_qadd(x::QAdd) = x
_as_qadd(x::SQA.QSym) = x * 1
_as_qadd(x::Number) = one(QAdd) * x
# After RHS unrolling, `new_op` may be a BasicSymbolic sum-of-averages
# expression rather than a single `QField`. We don't carry that through
# `eqs.operators` (which is `Vector{QAdd}`); placeholder with the empty
# QAdd so downstream code that doesn't use operators post-evaluate stays
# unaffected.
_as_qadd(x::SymbolicUtils.BasicSymbolic) = QAdd(SQA.QTermDict(), SQA.Index[])
