using QuantumCumulants
using Symbolics
using Polyester

const QC = QuantumCumulants
const SQA = QuantumCumulants.SecondQuantizedAlgebra

function has_transition(q::QC.QAdd)
    for (term, _) in q.arguments, op in term.ops
        SQA.is_transition(op) && return true
    end
    return false
end

function skip_ground_derive(op::QC.QAdd, sys, ctx::QC.CanonCtx)
    op_drift = QC._operator_rhs(
        sys.direction, op, im * sys.hamiltonian,
        sys.jumps, sys.jumps_dagger, sys.rates,
    )
    op_drift = SQA.expand_completeness(op_drift)
    op_drift = QC._assume_distinct_atom_indices(op_drift, QC._distinct_atom_indices([op]))
    drift = Symbolics.Num(QC.average_and_truncate(op_drift, sys.order, sys.mix_choice, ctx))
    has_transition(op_drift) && (drift = QC._reduce_ground_in_drift(drift))
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

function threaded_frontier(derive_one)
    return function(keys, sys, ctx)
        out = Vector{QC.NodeData}(undef, length(keys))
        Polyester.@batch per=thread for i in eachindex(keys)
            @inbounds out[i] = derive_one(keys[i], sys, ctx)
        end
        return out
    end
end

baseline = threaded_frontier(QC.derive)
skip = threaded_frontier(skip_ground_derive)
QC._closure(ising_graph(2), baseline)
QC._closure(ising_graph(2), skip)

for order in (2, 3, 4)
    GC.gc()
    bt = @timed QC._closure(ising_graph(order), baseline)
    GC.gc()
    st = @timed QC._closure(ising_graph(order), skip)
    exact = graphs_exact(bt.value, st.value)
    println(
        "GROUND_SKIP_THREADED order=$order exact=$exact threads=$(Threads.nthreads()) " *
        "states=$(length(bt.value.nodes)) baseline_time=$(round(bt.time; digits=6)) " *
        "skip_time=$(round(st.time; digits=6)) speedup=$(round(bt.time / st.time; digits=4)) " *
        "baseline_bytes=$(bt.bytes) skip_bytes=$(st.bytes) " *
        "alloc_ratio=$(round(st.bytes / bt.bytes; digits=5))",
    )
    exact || error("threaded ground-fold skip changed graph semantics at order $order")
end
