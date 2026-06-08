const TruncOrder = Union{Nothing, Int, Vector{Int}}

function average_and_truncate(R::QAdd, order::TruncOrder, mix_choice, ctx::CanonCtx)
    acc = 0
    for (term, c) in R.arguments
        # Index-dependent coefficients (Γ(i,j)/Ω(i,j)) must ride inside the sum so the
        # diagonal split substitutes them; scalar coefficients stay outside.
        if !isempty(term.ops) && !isempty(_coeff_scope_indices(c, R.indices))
            acc = acc + cumulant_expansion(
                _scoped_average_coeff(c, term.ops, term.ne, R.indices), order; mix_choice,
            )
        else
            acc = acc + _im_form(c) * _truncate_term(term.ops, term.ne, R.indices, order, mix_choice)
        end
    end
    return acc
end

"""Sum-scope indices the coefficient `c` depends on (empty for scalar coefficients)."""
function _coeff_scope_indices(c, scope)
    isempty(scope) && return SQA.Index[]
    vars = _coeff_vars(c)
    isempty(vars) && return SQA.Index[]
    out = SQA.Index[]
    for idx in scope
        isym = SymbolicUtils.unwrap(idx.sym)
        any(v -> isequal(v, isym), vars) && push!(out, idx)
    end
    return out
end

"""
Variables a coefficient depends on. `Num`/`Complex{Num}` are both `<: Number`, so test
the symbolic carrier explicitly rather than dispatching on `Number`.
"""
function _coeff_vars(c)
    cc = c isa Complex ? (real(c) + imag(c)) : c
    (cc isa Symbolics.Num || cc isa SymbolicUtils.BasicSymbolic) || return ()
    return Symbolics.get_variables(cc)
end

"""
Average a block with its coefficient inside the sum, scoping over indices used by the
operators or the coefficient so the diagonal split substitutes the coefficient.
"""
function _scoped_average_coeff(c, ops::AbstractVector{<:SQA.QSym}, ne, scope)
    cblock = c * reduce(*, ops)
    cblock isa QAdd || return average(cblock)
    # Carry NE on `cblock` (always a QAdd); a single-op block reduces to a QSym with
    # no NE slot, losing a `Σ_{j≠ext}` constraint.
    isempty(ne) || (cblock = _carry_ne(cblock, ne, scope))
    used = Set{SQA.Index}()
    for (t, _) in cblock.arguments, o in t.ops
        SQA.has_index(o.index) && push!(used, o.index)
    end
    for idx in _coeff_scope_indices(c, scope)
        push!(used, idx)
    end
    block_scope = SQA.Index[i for i in scope if i in used]
    isempty(block_scope) && return average(cblock)
    return average(SQA.Σ(cblock, block_scope[1], block_scope[2:end]...))
end

"""
Convert `Complex` coefficient literals to the factored `re + im*Symbolics.IM` form
(SymbolicUtils won't unify a raw `complex(0, …)` with what `average` emits).
"""
_im_form(c) = c
_im_form(c::Complex) = real(c) + imag(c) * Symbolics.IM
average_and_truncate(R::SQA.QField, order::TruncOrder, mix_choice, ::CanonCtx) =
    order === nothing ? average(R) : cumulant_expansion(average(R), order; mix_choice)

function _truncate_term(ops::AbstractVector{<:SQA.QSym}, ne, scope, order, mix_choice)
    # Average WITH the sum scope first so SQA's diagonal split collapses same-index
    # operator pairs, THEN cumulant-truncate; truncating the raw product first splits
    # them into different blocks and the collapse never fires.
    return cumulant_expansion(_scoped_average(ops, ne, scope), order; mix_choice)
end

"""
Average a single block, attaching only the scope indices its ops use, routed through
`SQA.Σ` so the off-diagonal/diagonal split runs.
"""
function _scoped_average(ops::AbstractVector{<:SQA.QSym}, ne, scope)
    isempty(ops) && return 1   # empty operator product is the identity, ⟨I⟩ = 1
    block = reduce(*, ops)
    block isa QAdd && !isempty(ne) && (block = _carry_ne(block, ne, scope))
    used = Set{SQA.Index}()
    for o in ops
        SQA.has_index(o.index) && push!(used, o.index)
    end
    block_scope = SQA.Index[i for i in scope if i in used]
    isempty(block_scope) && return average(block)
    return average(SQA.Σ(block, block_scope[1], block_scope[2:end]...))
end

"""
Reattach the term's NE pairs so the diagonal split sees them: keep a pair when both
indices are on the block, or one is on the block and its partner is external.
"""
function _carry_ne(block::QAdd, ne, scope)
    present = Set{SQA.Index}()
    for (term, _) in block.arguments, o in term.ops
        SQA.has_index(o.index) && push!(present, o.index)
    end
    scopeset = Set{SQA.Index}(scope)
    kept = SQA.NonEqualPair[
        p for p in ne if
            (p[1] in present && p[2] in present) ||
            (p[1] in present && !(p[2] in scopeset)) ||
            (p[2] in present && !(p[1] in scopeset))
    ]
    isempty(kept) && return block
    out = SQA.QTermDict()
    for (term, c) in block.arguments
        out[SQA.QTerm(copy(term.ops), vcat(term.ne, kept))] = c
    end
    return SQA.QAdd(out, block.indices)
end

"""
Moment-level ground-projector reduction: rewrite each averaged moment containing a
ground-state projector `σ^gg` via the completeness identity `σ^gg = 1 − Σ_{k≠g} σ^kk`,
canonicalising onto the minimal `N−1` population basis. A no-op without a ground projector.
"""
_reduce_ground_in_drift(x::Symbolics.Num) =
    Symbolics.Num(_reduce_ground_in_drift(SymbolicUtils.unwrap(x)))
function _reduce_ground_in_drift(x)
    return mapleaves(x) do leaf
        op = undo_average(leaf)
        op isa QAdd || return leaf
        expanded = SQA.expand_completeness(op)
        isequal(expanded, op) ? leaf : SymbolicUtils.unwrap(average(expanded))
    end
end

struct NodeData
    drift::Symbolics.Num                      # faithful averaged-and-truncated RHS
    op_drift::QAdd                            # operator RHS (latex / inspection / re-truncation)
    noise::Union{Nothing, Symbolics.Num}      # averaged + truncated noise drift (optional)
    op_noise::Union{Nothing, Symbolics.Num}   # operator-level noise form: deferred (always `nothing`)
    order::Int                                # cached
    aon::Vector{Int}                          # cached acts_on
end

"""
Derive the averaged-and-truncated drift (and noise) for one operator: build the
operator RHS via `_operator_rhs`, apply the `σ^gg` fold, then cross to moments
via `average_and_truncate`.
"""
function derive(op::QAdd, sys, ctx::CanonCtx)
    op_drift = _operator_rhs(
        sys.direction, op, im * sys.hamiltonian,
        sys.jumps, sys.jumps_dagger, sys.rates
    )
    op_drift = SQA.expand_completeness(op_drift)
    # Assert the LHS operator's free atom-space indices distinct (a multi-atom moment
    # ⟨σ_i σ_j⟩ is a distinct-slot cumulant, i≠j); else diagonal-split terms fail to
    # collapse and leak spurious higher-order cumulants.
    op_drift = _assume_distinct_atom_indices(op_drift, _distinct_atom_indices([op]))
    drift = Symbolics.Num(average_and_truncate(op_drift, sys.order, sys.mix_choice, ctx))
    drift = _reduce_ground_in_drift(drift)

    if sys.efficiencies === nothing
        op_noise = nothing
        noise = nothing
    else
        # `_noise_builder` returns (operator, averaged) noise eqs; the operator-level
        # drift is deferred (mixed operator/average), so truncate the averaged one.
        _, noise_eqs = _noise_builder(sys.direction)(
            [op], sys.jumps,
            sys.jumps_dagger, sys.rates, sys.efficiencies
        )
        op_noise = nothing
        noise_rhs = noise_eqs[1].rhs
        noise = Symbolics.Num(
            sys.order === nothing ? noise_rhs :
                cumulant_expansion(noise_rhs, sys.order; mix_choice = sys.mix_choice),
        )
        noise = _reduce_ground_in_drift(noise)
    end

    return NodeData(drift, op_drift, noise, op_noise, get_order(op), SQA.acts_on(op))
end
