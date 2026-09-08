using QuantumCumulants
using Symbolics
using SymbolicUtils

const QC = QuantumCumulants
const SQA = QuantumCumulants.SecondQuantizedAlgebra

function ising_graph(order)
    N = 6
    h = ⊗([PauliSpace(Symbol(:spin, i)) for i in 1:N]...)
    σx(i) = Pauli(h, :σ, 1, i)
    σy(i) = Pauli(h, :σ, 2, i)
    σz(i) = Pauli(h, :σ, 3, i)
    σm(i) = (σx(i) - 1im * σy(i)) / 2
    @variables J hₓ γ
    H = -J * sum(σz(i) * σz(i + 1) for i in 1:(N - 1)) - hₓ * sum(σx(i) for i in 1:N)
    c_ops = [σm(i) for i in 1:N]
    return meanfield(
        [σz(i) for i in 1:N], H, c_ops;
        rates = [γ for _ in 1:N], order = order,
    ).graph
end

function graphs_exact(a, b)
    collect(keys(a.nodes)) == collect(keys(b.nodes)) || return false
    for k in keys(a.nodes)
        x, y = a.nodes[k], b.nodes[k]
        isequal(x.drift, y.drift) || return false
        isequal(x.op_drift, y.op_drift) || return false
        isequal(x.noise, y.noise) || return false
        isequal(x.op_noise, y.op_noise) || return false
        x.order == y.order || return false
        x.aon == y.aon || return false
    end
    return true
end

function serial_frontier(derive_one)
    return function(keys, sys, ctx)
        out = Vector{QC.NodeData}(undef, length(keys))
        @inbounds for i in eachindex(keys)
            out[i] = derive_one(keys[i], sys, ctx)
        end
        return out
    end
end

# ---- empty-scope fast path -----------------------------------------------------

function scoped_average_fast(ops::AbstractVector{<:SQA.QSym}, non_equal, scope)
    isempty(ops) && return 1
    block = reduce(*, ops)
    block isa QC.QAdd && !isempty(non_equal) &&
        (block = QC._carry_non_equal(block, non_equal, scope))
    isempty(scope) && return average(block)
    used = Set{SQA.Index}()
    for op in ops
        SQA.has_index(op.index) && push!(used, op.index)
    end
    block_scope = SQA.Index[i for i in scope if i in used]
    isempty(block_scope) && return average(block)
    return average(SQA.Σ(block, block_scope[1], block_scope[2:end]...))
end

function average_truncate_scopefast(R::QC.QAdd, order, mix_choice, ctx::QC.CanonCtx)
    acc = 0
    for (term, coeff) in R.arguments
        c = QC._coeff_num(coeff)
        QC._iszero_coeff(c) && continue
        if !isempty(QC._coeff_scope_indices(c, R.indices))
            acc = acc + QC.cumulant_expansion(
                QC._scoped_average_coeff(c, term.ops, term.ne, R.indices), order; mix_choice,
            )
        else
            acc = acc + QC._im_form(QC._truncate_coeff(c, order, mix_choice)) *
                QC.cumulant_expansion(
                    scoped_average_fast(term.ops, term.ne, R.indices), order; mix_choice,
                )
        end
    end
    return acc
end

# ---- one memo across the whole derive -----------------------------------------

function expand_average_shared(ops, order::Vector{Int}, mix_choice, memo)
    ops isa QC.QAdd || return average(ops)
    terms = Any[]
    for (term, coeff) in ops.arguments
        c = QC._im_form(QC._coeff_num(coeff))
        expanded = QC._expand_product!(term.ops, order, mix_choice, memo)
        push!(terms, QC._reattach_scope(c, expanded, ops.indices, term.ne))
    end
    return QC._bulk_add(terms)
end

function cumulant_shared(x, order::Vector{Int}, mix_choice, memo)
    x isa Number && return x
    x isa Symbolics.Num && return cumulant_shared(SymbolicUtils.unwrap(x), order, mix_choice, memo)
    x isa SymbolicUtils.BasicSymbolic || return x
    QC.get_order(x) <= minimum(order) && return x
    if QC._is_moment_unit(x)
        ops = SQA.undo_average(x)
        aons = SQA.acts_on(ops)
        ord = isempty(aons) ? maximum(order) : mix_choice(order[k] for k in aons)
        QC.get_order(x) <= ord && return x
        return expand_average_shared(ops, order, mix_choice, memo)
    end
    if SymbolicUtils.iscall(x) && QC._has_average(x)
        op = SymbolicUtils.operation(x)
        args = SymbolicUtils.arguments(x)
        new_args = Any[cumulant_shared(a, order, mix_choice, memo) for a in args]
        all(i -> new_args[i] === args[i], eachindex(args)) && return x
        return op(new_args...)
    end
    return x
end

function average_truncate_memo(R::QC.QAdd, order::Vector{Int}, mix_choice, ctx::QC.CanonCtx; scopefast=false)
    acc = 0
    memo = Dict{Any, Any}()
    for (term, coeff) in R.arguments
        c = QC._coeff_num(coeff)
        QC._iszero_coeff(c) && continue
        if !isempty(QC._coeff_scope_indices(c, R.indices))
            avg = QC._scoped_average_coeff(c, term.ops, term.ne, R.indices)
            acc = acc + cumulant_shared(avg, order, mix_choice, memo)
        else
            truncated_coeff = QC._im_form(QC._truncate_coeff(c, order, mix_choice))
            avg = scopefast ? scoped_average_fast(term.ops, term.ne, R.indices) :
                QC._scoped_average(term.ops, term.ne, R.indices)
            acc = acc + truncated_coeff * cumulant_shared(avg, order, mix_choice, memo)
        end
    end
    return acc
end

function derive_with(op::QC.QAdd, sys, ctx::QC.CanonCtx, mode::Symbol)
    op_drift = QC._operator_rhs(
        sys.direction, op, im * sys.hamiltonian,
        sys.jumps, sys.jumps_dagger, sys.rates,
    )
    op_drift = SQA.expand_completeness(op_drift)
    op_drift = QC._assume_distinct_atom_indices(op_drift, QC._distinct_atom_indices([op]))
    raw = if mode === :scope
        average_truncate_scopefast(op_drift, sys.order, sys.mix_choice, ctx)
    elseif mode === :memo
        average_truncate_memo(op_drift, sys.order, sys.mix_choice, ctx)
    elseif mode === :combined
        average_truncate_memo(op_drift, sys.order, sys.mix_choice, ctx; scopefast=true)
    else
        error("unknown mode")
    end
    drift = QC._reduce_ground_in_drift(Symbolics.Num(raw))
    return QC.NodeData(drift, op_drift, nothing, nothing, QC.get_order(op), SQA.acts_on(op))
end

scope_derive(op, sys, ctx) = derive_with(op, sys, ctx, :scope)
memo_derive(op, sys, ctx) = derive_with(op, sys, ctx, :memo)
combined_derive(op, sys, ctx) = derive_with(op, sys, ctx, :combined)

frontiers = (
    baseline = serial_frontier(QC.derive),
    scope = serial_frontier(scope_derive),
    memo = serial_frontier(memo_derive),
    combined = serial_frontier(combined_derive),
)

# Warm all implementations before measuring.
for f in values(frontiers)
    QC._closure(ising_graph(2), f)
end

for order in (2, 3, 4)
    results = Dict{Symbol, Any}()
    for name in (:baseline, :scope, :memo, :combined)
        GC.gc()
        results[name] = @timed QC._closure(ising_graph(order), getproperty(frontiers, name))
    end
    base = results[:baseline]
    for name in (:scope, :memo, :combined)
        r = results[name]
        exact = graphs_exact(base.value, r.value)
        println(
            "MOMENT_EXPERIMENT order=$order mode=$name exact=$exact states=$(length(base.value.nodes)) " *
            "baseline_time=$(round(base.time; digits=6)) time=$(round(r.time; digits=6)) " *
            "speedup=$(round(base.time / r.time; digits=4)) baseline_bytes=$(base.bytes) " *
            "bytes=$(r.bytes) alloc_ratio=$(round(r.bytes / base.bytes; digits=5))",
        )
        exact || error("$name experiment changed graph semantics at order $order")
    end
end
