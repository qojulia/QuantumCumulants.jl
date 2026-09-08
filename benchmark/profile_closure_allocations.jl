using QuantumCumulants
using Symbolics
using SymbolicUtils

const QC = QuantumCumulants
const SQA = QuantumCumulants.SecondQuantizedAlgebra

mutable struct StageTotals
    calls::Int
    time_ns::Dict{Symbol, Float64}
    bytes::Dict{Symbol, Int}
end
StageTotals() = StageTotals(0, Dict{Symbol, Float64}(), Dict{Symbol, Int}())

function record!(s::StageTotals, name::Symbol, timed)
    s.time_ns[name] = get(s.time_ns, name, 0.0) + timed.time * 1e9
    s.bytes[name] = get(s.bytes, name, 0) + timed.bytes
    return timed.value
end

function profiled_derive(op::QC.QAdd, sys, ctx::QC.CanonCtx, totals::StageTotals)
    totals.calls += 1

    t = @timed QC._operator_rhs(
        sys.direction, op, im * sys.hamiltonian,
        sys.jumps, sys.jumps_dagger, sys.rates,
    )
    op_drift = record!(totals, :operator_rhs, t)

    t = @timed SQA.expand_completeness(op_drift)
    op_drift = record!(totals, :expand_completeness, t)

    t = @timed QC._distinct_atom_indices([op])
    distinct = record!(totals, :distinct_indices, t)

    t = @timed QC._assume_distinct_atom_indices(op_drift, distinct)
    op_drift = record!(totals, :assume_distinct, t)

    t = @timed Symbolics.Num(QC.average_and_truncate(op_drift, sys.order, sys.mix_choice, ctx))
    drift = record!(totals, :average_truncate, t)

    t = @timed QC._reduce_ground_in_drift(drift)
    drift = record!(totals, :ground_reduction, t)

    if sys.efficiencies === nothing
        op_noise = nothing
        noise = nothing
    else
        t = @timed QC._noise_builder(sys.direction)(
            [op], sys.jumps, sys.jumps_dagger, sys.rates, sys.efficiencies,
        )
        noise_pair = record!(totals, :noise_builder, t)
        _, noise_eqs = noise_pair
        op_noise = nothing
        noise_rhs = noise_eqs[1].rhs
        t = @timed Symbolics.Num(
            sys.order === nothing ? noise_rhs :
            QC.cumulant_expansion(noise_rhs, sys.order; mix_choice = sys.mix_choice),
        )
        noise = record!(totals, :noise_truncate, t)
        t = @timed QC._reduce_ground_in_drift(noise)
        noise = record!(totals, :noise_ground_reduction, t)
    end

    t = @timed (QC.get_order(op), SQA.acts_on(op))
    order, aon = record!(totals, :metadata, t)
    return QC.NodeData(drift, op_drift, noise, op_noise, order, aon)
end

# Experimental path: retain one product-expansion memo for the whole operator -> moment
# conversion of a single node. Production currently creates a fresh memo for each separately
# encountered average/product expansion. This experiment stays benchmark-local until measured.
function _expand_average_shared(ops, order::Int, memo)
    ops isa QC.QAdd || return average(ops)
    terms = Any[]
    for (term, coeff) in ops.arguments
        c = QC._im_form(QC._coeff_num(coeff))
        expanded = QC._expand_product!(term.ops, order, memo)
        push!(terms, QC._reattach_scope(c, expanded, ops.indices, term.ne))
    end
    return QC._bulk_add(terms)
end

function _cumulant_expansion_shared(x, order::Int, memo)
    x isa Number && return x
    x isa Symbolics.Num && return _cumulant_expansion_shared(SymbolicUtils.unwrap(x), order, memo)
    x isa SymbolicUtils.BasicSymbolic || return x
    QC.get_order(x) <= order && return x
    if QC._is_moment_unit(x)
        return _expand_average_shared(SQA.undo_average(x), order, memo)
    end
    if SymbolicUtils.iscall(x) && QC._has_average(x)
        op = SymbolicUtils.operation(x)
        args = SymbolicUtils.arguments(x)
        new_args = Any[_cumulant_expansion_shared(a, order, memo) for a in args]
        all(i -> new_args[i] === args[i], eachindex(args)) && return x
        return op(new_args...)
    end
    return x
end

function _average_and_truncate_shared(R::QC.QAdd, order::Int, mix_choice, ctx::QC.CanonCtx)
    acc = 0
    memo = Dict{Any, Any}()
    for (term, coeff) in R.arguments
        c = QC._coeff_num(coeff)
        QC._iszero_coeff(c) && continue
        if !isempty(QC._coeff_scope_indices(c, R.indices))
            avg = QC._scoped_average_coeff(c, term.ops, term.ne, R.indices)
            acc = acc + _cumulant_expansion_shared(avg, order, memo)
        else
            truncated_coeff = QC._im_form(QC._truncate_coeff(c, order, mix_choice))
            avg = QC._scoped_average(term.ops, term.ne, R.indices)
            acc = acc + truncated_coeff * _cumulant_expansion_shared(avg, order, memo)
        end
    end
    return acc
end

function memo_derive(op::QC.QAdd, sys, ctx::QC.CanonCtx)
    op_drift = QC._operator_rhs(
        sys.direction, op, im * sys.hamiltonian,
        sys.jumps, sys.jumps_dagger, sys.rates,
    )
    op_drift = SQA.expand_completeness(op_drift)
    op_drift = QC._assume_distinct_atom_indices(op_drift, QC._distinct_atom_indices([op]))
    drift = Symbolics.Num(_average_and_truncate_shared(op_drift, sys.order, sys.mix_choice, ctx))
    drift = QC._reduce_ground_in_drift(drift)
    return QC.NodeData(drift, op_drift, nothing, nothing, QC.get_order(op), SQA.acts_on(op))
end

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
    eqs = meanfield(
        [σz(i) for i in 1:N], H, c_ops;
        rates = [γ for _ in 1:N], order = order,
    )
    return eqs.graph
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

function run_profile(order)
    g = ising_graph(order)
    totals = StageTotals()
    frontier = function(keys, sys, ctx)
        out = Vector{QC.NodeData}(undef, length(keys))
        @inbounds for i in eachindex(keys)
            out[i] = profiled_derive(keys[i], sys, ctx, totals)
        end
        return out
    end
    GC.gc()
    total = @timed QC._closure(g, frontier)
    completed = total.value
    println("PROFILE order=$order states=$(length(completed.nodes)) closure_calls=$(totals.calls)")
    println("TOTAL time_s=$(round(total.time; digits=6)) bytes=$(total.bytes)")
    sum_stage_bytes = sum(values(totals.bytes); init = 0)
    sum_stage_time = sum(values(totals.time_ns); init = 0.0)
    for name in sort!(collect(keys(totals.bytes)); by = x -> -totals.bytes[x])
        b = totals.bytes[name]
        tns = totals.time_ns[name]
        println(
            "STAGE name=$name bytes=$b byte_fraction=$(round(b / max(sum_stage_bytes, 1); digits=4)) " *
            "time_s=$(round(tns / 1e9; digits=6)) time_fraction=$(round(tns / max(sum_stage_time, 1); digits=4))",
        )
    end
    println("STAGE_SUM bytes=$sum_stage_bytes time_s=$(round(sum_stage_time / 1e9; digits=6))")
    return completed
end

function run_memo_experiment(order)
    # Fresh contexts for both sides avoid warming the monotonic canonicalization cache.
    gb = ising_graph(order)
    gm = ising_graph(order)
    baseline_frontier = serial_frontier(QC.derive)
    memo_frontier = serial_frontier(memo_derive)

    # Compilation warmup on independent small graphs before the measured closure.
    order == 2 && begin
        QC._closure(ising_graph(2), baseline_frontier)
        QC._closure(ising_graph(2), memo_frontier)
    end

    GC.gc()
    bt = @timed QC._closure(gb, baseline_frontier)
    GC.gc()
    mt = @timed QC._closure(gm, memo_frontier)
    exact = graphs_exact(bt.value, mt.value)
    println(
        "MEMO order=$order exact=$exact baseline_time=$(round(bt.time; digits=6)) " *
        "memo_time=$(round(mt.time; digits=6)) speedup=$(round(bt.time / mt.time; digits=4)) " *
        "baseline_bytes=$(bt.bytes) memo_bytes=$(mt.bytes) alloc_ratio=$(round(mt.bytes / bt.bytes; digits=5))",
    )
    exact || error("memo experiment changed graph semantics at order $order")
end

for order in (2, 3, 4)
    run_profile(order)
end
for order in (2, 3, 4)
    run_memo_experiment(order)
end
