using QuantumCumulants
using OhMyThreads
using Polyester
using Symbolics: @variables

const QC = QuantumCumulants
const SQA = QuantumCumulants.SecondQuantizedAlgebra
const N = 6

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

function derive_serial(keys, sys, ctx)
    return map(k -> QC.derive(k, sys, ctx), keys)
end

function derive_threads_static(keys, sys, ctx)
    out = Vector{QC.NodeData}(undef, length(keys))
    Threads.@threads :static for i in eachindex(keys)
        @inbounds out[i] = QC.derive(keys[i], sys, ctx)
    end
    return out
end

function derive_threads_dynamic(keys, sys, ctx)
    out = Vector{QC.NodeData}(undef, length(keys))
    Threads.@threads :dynamic for i in eachindex(keys)
        @inbounds out[i] = QC.derive(keys[i], sys, ctx)
    end
    return out
end

function derive_polyester(keys, sys, ctx)
    out = Vector{QC.NodeData}(undef, length(keys))
    Polyester.@batch per=thread for i in eachindex(keys)
        @inbounds out[i] = QC.derive(keys[i], sys, ctx)
    end
    return out
end

function derive_omt_static(keys, sys, ctx)
    return OhMyThreads.tmap(keys; scheduler = :static) do k
        QC.derive(k, sys, ctx)
    end
end

function derive_omt_dynamic(keys, sys, ctx)
    return OhMyThreads.tmap(keys; scheduler = :dynamic) do k
        QC.derive(k, sys, ctx)
    end
end

function closure_frontier(
        g::QC.MomentGraph,
        derive_frontier;
        filter = QC._alltrue,
        get_adjoints::Bool = false,
        foldable = QC._alltrue,
        max_iter::Int = 100_000,
    )
    ctx = g.ctx
    nodes = copy(g.nodes)
    seen = Set(keys(nodes))
    frontier = collect(keys(nodes))
    KeyT = eltype(frontier)
    processed = 0

    while !isempty(frontier)
        newkeys = KeyT[]
        for key in frontier
            processed >= max_iter && error(
                "closure did not close the hierarchy within $max_iter iterations",
            )
            processed += 1
            nd = nodes[key]
            for leaf in QC._drift_leaves(nd)
                op = SQA.undo_average(leaf)
                k = QC.canon_key(op, ctx)
                k in seen && continue
                kc = QC.canon_key(adjoint(op), ctx)
                if kc in seen && foldable(op)
                    push!(seen, k)
                    continue
                end
                filter(SQA.average(k)) || continue

                # Decide membership and ordering serially, exactly in frontier/node/leaf order.
                # Only the already-decided `derive` work below is parallelized.
                push!(seen, k)
                push!(newkeys, k)
                if kc != k && !(kc in seen)
                    if get_adjoints
                        push!(seen, kc)
                        push!(newkeys, kc)
                    elseif foldable(op)
                        push!(seen, kc)
                    end
                end
            end
        end

        derived = derive_frontier(newkeys, g.sys, ctx)
        @inbounds for i in eachindex(newkeys)
            nodes[newkeys[i]] = derived[i]
        end
        frontier = newkeys
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

const METHODS = (
    ("frontier serial", derive_serial),
    ("Threads static", derive_threads_static),
    ("Threads dynamic", derive_threads_dynamic),
    ("Polyester", derive_polyester),
    ("OhMyThreads static", derive_omt_static),
    ("OhMyThreads dynamic", derive_omt_dynamic),
)

println("Julia threads: ", Threads.nthreads())
println("warming schedulers on order 2")
warm = ising_graph(2)
warm_ref = QC.closure(warm; get_adjoints = false)
for (name, derive_frontier) in METHODS
    candidate = closure_frontier(warm, derive_frontier; get_adjoints = false)
    assert_same_graph(warm_ref, candidate)
    println("  warm $name: exact")
end

for order in (3, 4)
    println("\norder $order")
    graph = ising_graph(order)
    reference, queue_time, queue_bytes =
        timed(() -> QC.closure(graph; get_adjoints = false))
    println(
        "  queue cursor: time=$(queue_time)s bytes=$(queue_bytes) states=$(length(reference.nodes))",
    )

    for (name, derive_frontier) in METHODS
        candidate, elapsed, bytes = timed(
            () -> closure_frontier(graph, derive_frontier; get_adjoints = false),
        )
        assert_same_graph(reference, candidate)
        speedup = queue_time / elapsed
        println("  $name: time=$(elapsed)s bytes=$(bytes) speedup=$(speedup)x exact=true")
    end
end
