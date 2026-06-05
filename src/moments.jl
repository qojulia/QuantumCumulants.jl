# The operator-to-moment crossing (spec Section 3.3). One fused pass: never a
# bare `average(_)` of the whole RHS. Walk term by term, holding the sum scope
# `R.indices` and each `term.ne`; for an above-cap block, Wick-expand and average
# each surviving block with its sum scope reattached via SQA.Σ (Correction 3).

const TruncOrder = Union{Nothing, Int, Vector{Int}}

function average_and_truncate(R::QAdd, order::TruncOrder, mix_choice, ctx::CanonCtx)
    acc = 0
    for (term, c) in R.arguments
        acc = acc + _im_form(c) * _truncate_term(term.ops, term.ne, R.indices, order, mix_choice)
    end
    return acc
end

# QAdd coefficients are `Complex` literals (e.g. `complex(0, g)` from `im*H`).
# SQA's `average` emits the factored `re + im*Symbolics.IM` form; SymbolicUtils
# does not unify a raw `complex(0, …)` literal with that form (see CLAUDE.md), so
# convert each coefficient the same way before it enters the c-number RHS.
_im_form(c) = c
_im_form(c::Complex) = real(c) + imag(c) * Symbolics.IM
# A bare operator RHS (QSym) with no QAdd structure: average it directly.
average_and_truncate(R::SQA.QField, order::TruncOrder, mix_choice, ::CanonCtx) =
    order === nothing ? average(R) : cumulant_expansion(average(R), order; mix_choice)

function _truncate_term(ops::AbstractVector{<:SQA.QSym}, ne, scope, order, mix_choice)
    # Average WITH the sum scope first, so SQA's diagonal split collapses same-
    # index operator pairs (e.g. the dissipator `σ_k₂₁ … σ_k₁₂` collapses for
    # k=i,j and vanishes off-diagonal), THEN cumulant-truncate the resulting
    # moments. Wicking the RAW operator product first (the earlier fused path)
    # split `σ_k₂₁`/`σ_k₁₂` into different cumulant blocks, so the collapse never
    # fired and multi-atom dissipative drifts grew spurious higher-order terms
    # (e.g. `⟨σ_i₂₂ σ_j₂₂⟩` decaying as cumulant junk instead of the exact -2γ).
    return cumulant_expansion(_scoped_average(ops, ne, scope), order; mix_choice)
end

# Average a single block, attaching ONLY the scope indices the block's ops use,
# routed through SQA.Σ so the off-diagonal/diagonal split runs (Correction 3).
function _scoped_average(ops::AbstractVector{<:SQA.QSym}, ne, scope)
    isempty(ops) && return 1   # empty operator product is the identity, ⟨I⟩ = 1
    block = reduce(*, ops)
    block isa QAdd && !isempty(ne) && (block = _carry_ne(block, ne))
    used = Set{SQA.Index}()
    for o in ops
        SQA.has_index(o.index) && push!(used, o.index)
    end
    block_scope = SQA.Index[i for i in scope if i in used]
    isempty(block_scope) && return average(block)
    return average(SQA.Σ(block, block_scope[1], block_scope[2:end]...))
end

# Reattach the term's NE pairs (restricted to indices present on the block) so
# the diagonal split sees them. Only the pairs both of whose indices the block
# carries are meaningful.
function _carry_ne(block::QAdd, ne)
    present = Set{SQA.Index}()
    for (term, _) in block.arguments, o in term.ops
        SQA.has_index(o.index) && push!(present, o.index)
    end
    kept = SQA.NonEqualPair[p for p in ne if p[1] in present && p[2] in present]
    isempty(kept) && return block
    out = SQA.QTermDict()
    for (term, c) in block.arguments
        out[SQA.QTerm(copy(term.ops), vcat(term.ne, kept))] = c
    end
    return SQA.QAdd(out, block.indices)
end

# Moment-level ground-projector reduction (spec Task 3). Replace each averaged
# moment ⟨O⟩ whose operator product contains a ground-state projector σ^gg by
# `average(expand_completeness(O))`, the EXACT completeness identity
# `σ^gg = 1 − Σ_{k≠g} σ^kk` (e.g. ⟨σ^11⟩ → 1 − ⟨σ^22⟩, and the multilevel
# analogue). Dynamics are unchanged; it only canonicalises the moment onto the
# minimal N−1 population basis, so a ground-projector moment never becomes its
# own state. Applied at the MOMENT (averaged) level so operator expressions keep
# σ^gg atomic (the SQA keep-atomic invariant, CLAUDE.md). Applied UNCONDITIONALLY
# (not gated on ctx.population): `expand_completeness` is identity on operators
# with no ground projector, so this is a no-op on every other leaf. This catches
# σ^gg terms surfaced by `average_and_truncate`'s diagonal split, which the
# one-shot operator-level `expand_completeness` in `derive` runs too early to see.
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
    op_noise::Union{Nothing, Symbolics.Num}   # operator-level noise form: deferred (always
    # `nothing` in Phase 0). The noise operator drift
    # is operator-valued with average coefficients, so
    # it fits neither QAdd nor Num; a later phase stores
    # it once the noise machinery yields a clean form.
    order::Int                                # cached
    aon::Vector{Int}                          # cached acts_on
end

# The single operator-algebra entry point. Reuses meanfield.jl's `_operator_rhs`
# (coherent + Lindblad recycling, carrying direction/collective/indexed decay)
# and noise.jl's noise builders verbatim; applies the σ^gg fold and the
# population closure policy, then crosses to moments via average_and_truncate.
function derive(op::QAdd, sys, ctx::CanonCtx)
    op_drift = _operator_rhs(
        sys.direction, op, im * sys.hamiltonian,
        sys.jumps, sys.jumps_dagger, sys.rates
    )
    op_drift = SQA.expand_completeness(op_drift)
    # Assert the LHS operator's free atom-space indices are mutually distinct
    # (a multi-atom moment ⟨σ_i σ_j⟩ is by construction a distinct-slot cumulant,
    # i≠j). Without this, SQA leaves the undetermined pair in physical order, so a
    # dissipator diagonal-split term like `σ_j σ_i σ_j` (the k=j contribution of
    # `Σ_k σ_k σ_i σ_j`) fails to collapse to `σ_i σ_j` and leaks spurious
    # higher-order cumulants into multi-atom drifts. Previously gated on
    # `ctx.population`; the distinctness holds for every multi-atom moment.
    op_drift = _assume_distinct_atom_indices(op_drift, _distinct_atom_indices([op]))
    drift = Symbolics.Num(average_and_truncate(op_drift, sys.order, sys.mix_choice, ctx))
    drift = _reduce_ground_in_drift(drift)

    if sys.efficiencies === nothing
        op_noise = nothing
        noise = nothing
    else
        # `_noise_builder` returns (operator noise eqs, averaged noise eqs). The
        # operator-level noise drift carries the measurement-backaction ⟨J†+J⟩
        # averages in its coefficients (an operator-valued mixed symbolic), so it
        # fits neither QAdd nor Num and is deferred; truncate the builder's
        # averaged drift, exactly as meanfield's `_finalize_noise_eqs`.
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
