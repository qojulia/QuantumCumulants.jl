using QuantumCumulants
using Polyester
using Symbolics

const QC = QuantumCumulants
const SQA = QuantumCumulants.SecondQuantizedAlgebra
const GROUND_LOCK = ReentrantLock()

function may_ground(q::QC.QAdd)
    for (term, _) in q.arguments, op in term.ops
        SQA.is_transition(op) && return true
    end
    return false
end

locked_ground(x) = lock(GROUND_LOCK) do
    QC._reduce_ground_in_drift(x)
end

function safe_derive(op::QC.QAdd, sys, ctx::QC.CanonCtx)
    op_drift = QC._operator_rhs(
        sys.direction, op, im * sys.hamiltonian,
        sys.jumps, sys.jumps_dagger, sys.rates,
    )
    op_drift = SQA.expand_completeness(op_drift)
    op_drift = QC._assume_distinct_atom_indices(op_drift, QC._distinct_atom_indices([op]))
    drift = Symbolics.Num(QC.average_and_truncate(op_drift, sys.order, sys.mix_choice, ctx))
    may_ground(op_drift) && (drift = locked_ground(drift))
    return QC.NodeData(drift, op_drift, nothing, nothing, QC.get_order(op), SQA.acts_on(op))
end

function serial_frontier(keys, sys, ctx)
    out = Vector{QC.NodeData}(undef, length(keys))
    @inbounds for i in eachindex(keys)
        out[i] = safe_derive(keys[i], sys, ctx)
    end
    return out
end

function threaded_frontier(keys, sys, ctx)
    out = Vector{QC.NodeData}(undef, length(keys))
    Polyester.@batch per = thread for i in eachindex(keys)
        @inbounds out[i] = safe_derive(keys[i], sys, ctx)
    end
    return out
end

function ising_graph(order)
    N = 6
    h = ⊗([PauliSpace(Symbol(:spin, i)) for i in 1:N]...)
    sx(i) = Pauli(h, :σ, 1, i)
    sy(i) = Pauli(h, :σ, 2, i)
    sz(i) = Pauli(h, :σ, 3, i)
    sm(i) = (sx(i) - 1im * sy(i)) / 2
    @variables J hx γ
    H = -J * sum(sz(i) * sz(i + 1) for i in 1:(N - 1)) - hx * sum(sx(i) for i in 1:N)
    return meanfield([sz(i) for i in 1:N], H, [sm(i) for i in 1:N]; rates = fill(γ, N), order = order).graph
end

function indexed_transition_graph()
    h = NLevelSpace(:atom, 2)
    @variables N
    i = Index(h, :i, N, h)
    j = Index(h, :j, N, h)
    σ(a, b, k) = IndexedOperator(Transition(h, :σ, a, b), k)
    return meanfield(
        [σ(1, 2, i), σ(2, 2, i)], 0 * Σ(σ(2, 2, i), i), [σ(1, 2, i)];
        rates = [DoubleIndexedVariable(:Γ, i, j)], order = 2,
    ).graph
end

function exact_graph(a, b)
    ka, kb = collect(keys(a.nodes)), collect(keys(b.nodes))
    length(ka) == length(kb) || return false
    all(i -> isequal(ka[i], kb[i]), eachindex(ka)) || return false
    for i in eachindex(ka)
        x, y = a.nodes[ka[i]], b.nodes[kb[i]]
        isequal(x.drift, y.drift) || return false
        isequal(x.op_drift, y.op_drift) || return false
        isequal(x.noise, y.noise) || return false
        isequal(x.op_noise, y.op_noise) || return false
        x.order == y.order || return false
        x.aon == y.aon || return false
    end
    return true
end

# Warm worker code on a small closure.
QC._closure(ising_graph(2), threaded_frontier)

for order in (3, 4)
    GC.gc()
    s = @timed QC._closure(ising_graph(order), serial_frontier)
    GC.gc()
    t = @timed QC._closure(ising_graph(order), threaded_frontier)
    exact = exact_graph(s.value, t.value)
    println("GROUND_SAFE_THREADED order=$order exact=$exact states=$(length(s.value.nodes)) serial_time=$(round(s.time; digits=6)) threaded_time=$(round(t.time; digits=6)) speedup=$(round(s.time/t.time; digits=4)) serial_bytes=$(s.bytes) threaded_bytes=$(t.bytes) alloc_ratio=$(round(t.bytes/s.bytes; digits=5))")
    exact || error("guarded serial/threaded mismatch")
end

reference = QC._closure(indexed_transition_graph(), serial_frontier; get_adjoints = false)
for rep in 1:25
    candidate = QC._closure(indexed_transition_graph(), threaded_frontier; get_adjoints = false)
    exact_graph(reference, candidate) || error("transition threaded mismatch on repetition $rep")
end
println("GROUND_SAFE_TRANSITION exact=true repetitions=25 states=$(length(reference.nodes))")
