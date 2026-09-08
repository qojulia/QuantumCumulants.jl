using QuantumCumulants
using Statistics: median
using Symbolics: @variables

const QC = QuantumCumulants
const SQA = QuantumCumulants.SecondQuantizedAlgebra
const N = 6
const NSAMPLES = 3

function ising_graph(order)
    h = ⊗([PauliSpace(Symbol(:spin, i)) for i in 1:N]...)
    σx(i) = Pauli(h, :σ, 1, i)
    σy(i) = Pauli(h, :σ, 2, i)
    σz(i) = Pauli(h, :σ, 3, i)
    σm(i) = (σx(i) - 1im * σy(i)) / 2
    @variables J hx γ
    H = -J * sum(σz(i) * σz(i + 1) for i in 1:(N - 1)) - hx * sum(σx(i) for i in 1:N)
    eqs = meanfield(
        [σz(i) for i in 1:N],
        H,
        [σm(i) for i in 1:N];
        rates = fill(γ, N),
        order,
    )
    return eqs.graph
end

function derive_with_iH(op::SQA.QAdd, sys, ctx::QC.CanonCtx, iH)
    op_drift = QC._operator_rhs(
        sys.direction, op, iH,
        sys.jumps, sys.jumps_dagger, sys.rates,
    )
    op_drift = SQA.expand_completeness(op_drift)
    op_drift = QC._assume_distinct_atom_indices(
        op_drift,
        QC._distinct_atom_indices([op]),
    )
    drift = Symbolics.Num(QC.average_and_truncate(op_drift, sys.order, sys.mix_choice, ctx))
    drift = QC._reduce_ground_in_drift(drift)

    if sys.efficiencies === nothing
        op_noise = nothing
        noise = nothing
    else
        _, noise_eqs = QC._noise_builder(sys.direction)(
            [op], sys.jumps,
            sys.jumps_dagger, sys.rates, sys.efficiencies,
        )
        op_noise = nothing
        noise_rhs = noise_eqs[1].rhs
        noise = Symbolics.Num(
            sys.order === nothing ? noise_rhs :
                QC.cumulant_expansion(noise_rhs, sys.order; mix_choice = sys.mix_choice),
        )
        noise = QC._reduce_ground_in_drift(noise)
    end

    return QC.NodeData(
        drift,
        op_drift,
        noise,
        op_noise,
        QC.get_order(op),
        SQA.acts_on(op),
    )
end

canon_default(op, g) = QC.canon_key(op, g.ctx)
canon_reuse_treatments(op::SQA.QAdd, g) = QC._treatment_key(op, g.ctx, g.treatments)
canon_reuse_treatments(op, g) = op

function close_custom(
        g::QC.MomentGraph;
        hoist_iH::Bool = false,
        reuse_treatments::Bool = false,
        filter = QC._alltrue,
        get_adjoints::Bool = false,
        foldable = QC._alltrue,
        max_iter::Int = 100_000,
    )
    ctx = g.ctx
    nodes = copy(g.nodes)
    seen = Set(keys(nodes))
    pending = collect(keys(nodes))
    cursor = 1
    iH = hoist_iH ? im * g.sys.hamiltonian : nothing

    keyfn = reuse_treatments ? canon_reuse_treatments : canon_default

    while cursor <= length(pending)
        cursor > max_iter && error(
            "closure did not close the hierarchy within $max_iter iterations",
        )
        nd = nodes[pending[cursor]]
        cursor += 1
        for leaf in QC._drift_leaves(nd)
            op = SQA.undo_average(leaf)
            k = keyfn(op, g)
            k in seen && continue
            kc = keyfn(adjoint(op), g)
            if kc in seen && foldable(op)
                push!(seen, k)
                continue
            end
            filter(SQA.average(k)) || continue

            nodes[k] = hoist_iH ? derive_with_iH(k, g.sys, ctx, iH) : QC.derive(k, g.sys, ctx)
            push!(seen, k)
            push!(pending, k)

            if kc != k && !(kc in seen)
                if get_adjoints
                    nodes[kc] = hoist_iH ?
                        derive_with_iH(kc, g.sys, ctx, iH) : QC.derive(kc, g.sys, ctx)
                    push!(seen, kc)
                    push!(pending, kc)
                elseif foldable(op)
                    push!(seen, kc)
                end
            end
        end
    end

    return QC.MomentGraph(nodes, g.sys, ctx, g.treatments)
end

function assert_same_graph(reference, candidate)
    rkeys = collect(keys(reference.nodes))
    ckeys = collect(keys(candidate.nodes))
    @assert length(rkeys) == length(ckeys)
    @assert all(i -> isequal(rkeys[i], ckeys[i]), eachindex(rkeys))
    for i in eachindex(rkeys)
        r = reference.nodes[rkeys[i]]
        c = candidate.nodes[ckeys[i]]
        @assert isequal(r.drift, c.drift)
        @assert isequal(r.op_drift, c.op_drift)
        @assert isequal(r.noise, c.noise)
        @assert isequal(r.op_noise, c.op_noise)
        @assert r.order == c.order
        @assert r.aon == c.aon
    end
    return nothing
end

function timed(f)
    GC.gc()
    result = @timed f()
    return result.value, result.time, result.bytes
end

run_baseline(g) = QC.closure(g; get_adjoints = false)
run_iH(g) = close_custom(g; hoist_iH = true, reuse_treatments = false)
run_treatments(g) = close_custom(g; hoist_iH = false, reuse_treatments = true)
run_both(g) = close_custom(g; hoist_iH = true, reuse_treatments = true)

const METHODS = (
    ("baseline cursor", run_baseline),
    ("hoist iH", run_iH),
    ("reuse free treatments", run_treatments),
    ("hoist both", run_both),
)

println("samples per method: ", NSAMPLES)
println("warming on fresh order-2 graphs")
warm_ref = run_baseline(ising_graph(2))
for (name, method) in METHODS
    candidate = method(ising_graph(2))
    assert_same_graph(warm_ref, candidate)
    println("  warm $name: exact")
end

for order in (3, 4)
    println("\norder $order")
    reference = run_baseline(ising_graph(order))
    println("  states=$(length(reference.nodes))")
    results = Dict{String, Tuple{Float64, Float64}}()

    for (name, method) in METHODS
        times = Float64[]
        bytes = Float64[]
        for _ in 1:NSAMPLES
            graph = ising_graph(order)
            candidate, elapsed, allocated = timed(() -> method(graph))
            assert_same_graph(reference, candidate)
            push!(times, elapsed)
            push!(bytes, allocated)
        end
        mt = median(times)
        mb = median(bytes)
        results[name] = (mt, mb)
        println("  $name: median_time=$(mt)s median_bytes=$(mb) exact=true")
    end

    baseline_time, baseline_bytes = results["baseline cursor"]
    println("  relative_to_baseline")
    for (name, _) in METHODS
        t, b = results[name]
        println("    $name: speedup=$(baseline_time / t)x bytes_ratio=$(b / baseline_bytes)")
    end
end
