using QuantumCumulants
using Symbolics

const QC = QuantumCumulants
const SQA = QuantumCumulants.SecondQuantizedAlgebra

function bulk_average_and_truncate(R::QC.QAdd, order, mix_choice, ctx::QC.CanonCtx)
    terms = Any[]
    sizehint!(terms, length(R.arguments))
    for (term, coeff) in R.arguments
        c = QC._coeff_num(coeff)
        QC._iszero_coeff(c) && continue
        value = if !isempty(QC._coeff_scope_indices(c, R.indices))
            QC.cumulant_expansion(
                QC._scoped_average_coeff(c, term.ops, term.ne, R.indices), order; mix_choice,
            )
        else
            QC._im_form(QC._truncate_coeff(c, order, mix_choice)) *
                QC._truncate_term(term.ops, term.ne, R.indices, order, mix_choice)
        end
        push!(terms, value)
    end
    return QC._bulk_add(terms)
end

function bulk_derive(op::QC.QAdd, sys, ctx::QC.CanonCtx)
    op_drift = QC._operator_rhs(
        sys.direction, op, im * sys.hamiltonian,
        sys.jumps, sys.jumps_dagger, sys.rates,
    )
    op_drift = SQA.expand_completeness(op_drift)
    op_drift = QC._assume_distinct_atom_indices(op_drift, QC._distinct_atom_indices([op]))
    drift = Symbolics.Num(bulk_average_and_truncate(op_drift, sys.order, sys.mix_choice, ctx))
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

baseline = serial_frontier(QC.derive)
bulk = serial_frontier(bulk_derive)
QC._closure(ising_graph(2), baseline)
QC._closure(ising_graph(2), bulk)

for order in (2, 3, 4)
    GC.gc()
    bt = @timed QC._closure(ising_graph(order), baseline)
    GC.gc()
    xt = @timed QC._closure(ising_graph(order), bulk)
    exact = graphs_exact(bt.value, xt.value)
    println(
        "BULK_AVERAGE order=$order exact=$exact states=$(length(bt.value.nodes)) " *
        "baseline_time=$(round(bt.time; digits=6)) bulk_time=$(round(xt.time; digits=6)) " *
        "speedup=$(round(bt.time / xt.time; digits=4)) baseline_bytes=$(bt.bytes) " *
        "bulk_bytes=$(xt.bytes) alloc_ratio=$(round(xt.bytes / bt.bytes; digits=5))",
    )
    exact || error("bulk average changed graph semantics at order $order")
end
