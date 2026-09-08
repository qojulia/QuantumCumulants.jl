using QuantumCumulants
using OhMyThreads
using Polyester
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

function indexed_collective_graph()
    h = NLevelSpace(:atom, 2)
    @variables Nat
    i = Index(h, :i, Nat, h)
    j = Index(h, :j, Nat, h)
    σ(x, y, k) = IndexedOperator(Transition(h, :σ, x, y), k)
    eqs = meanfield(
        [σ(1, 2, i), σ(2, 2, i)],
        0 * Σ(σ(2, 2, i), i),
        [σ(1, 2, i)];
        rates = [DoubleIndexedVariable(:Γ, i, j)],
        order = 2,
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

function closure_popfirst(
        g::QC.MomentGraph;
        filter = QC._alltrue,
        get_adjoints::Bool = false,
        foldable = QC._alltrue,
        max_iter::Int = 100_000,
    )
    ctx = g.ctx
    nodes = copy(g.nodes)
    seen = Set(keys(nodes))
    pending = collect(keys(nodes))
    iters = 0
    while !isempty(pending)
        iters >= max_iter && error(
            "closure did not close the hierarchy within $max_iter iterations " *
                "($(length(pending)) moments still pending); the system may not close.",
        )
        iters += 1
        nd = nodes[popfirst!(pending)]
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
            nodes[k] = QC.derive(k, g.sys, ctx)
            push!(seen, k)
            push!(pending, k)
            if kc != k && !(kc in seen)
                if get_adjoints
                    nodes[kc] = QC.derive(kc, g.sys, ctx)
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

function closure_frontier(
        g::QC.MomentGraph,
        derive_frontier;
        filter = QC._alltrue,
        get_adjoints::Bool = false,
        foldable = QC._alltrue,
        max_iter::Int = 100_000,
        assert_intern_stable::Bool = true,
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

                # Discovery, deduplication, and ordering remain serial and deterministic.
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

        # The worker phase must not mint SQA names/ranges. SQA's intern tables allow
        # concurrent canonicalisation only after construction has populated them.
        names_before = length(SQA.NAME_BY_ID)
        ranges_before = length(SQA.RANGE_BY_ID)
        derived = derive_frontier(newkeys, g.sys, ctx)
        if assert_intern_stable
            @assert length(SQA.NAME_BY_ID) == names_before
            @assert length(SQA.RANGE_BY_ID) == ranges_before
        end

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

run_popfirst(g) = closure_popfirst(g; get_adjoints = false)
run_cursor(g) = QC.closure(g; get_adjoints = false)
run_frontier_serial(g) = closure_frontier(g, derive_serial; get_adjoints = false)
run_threads_static(g) = closure_frontier(g, derive_threads_static; get_adjoints = false)
run_threads_dynamic(g) = closure_frontier(g, derive_threads_dynamic; get_adjoints = false)
run_polyester(g) = closure_frontier(g, derive_polyester; get_adjoints = false)
run_omt_static(g) = closure_frontier(g, derive_omt_static; get_adjoints = false)
run_omt_dynamic(g) = closure_frontier(g, derive_omt_dynamic; get_adjoints = false)

const METHODS = (
    ("queue popfirst", run_popfirst),
    ("queue cursor", run_cursor),
    ("frontier serial", run_frontier_serial),
    ("Threads static", run_threads_static),
    ("Threads dynamic", run_threads_dynamic),
    ("Polyester", run_polyester),
    ("OhMyThreads static", run_omt_static),
    ("OhMyThreads dynamic", run_omt_dynamic),
)

println("Julia threads: ", Threads.nthreads())
println("samples per method: ", NSAMPLES)
println("warming methods on fresh order-2 graphs")
warm_ref = run_cursor(ising_graph(2))
for (name, method) in METHODS
    candidate = method(ising_graph(2))
    assert_same_graph(warm_ref, candidate)
    println("  warm $name: exact")
end

for order in (3, 4)
    println("\norder $order")
    reference = run_cursor(ising_graph(order))
    println("  states=$(length(reference.nodes))")
    results = Dict{String, Tuple{Float64, Float64}}()

    for (name, method) in METHODS
        times = Float64[]
        bytes = Float64[]
        for sample in 1:NSAMPLES
            graph = ising_graph(order) # fresh CanonCtx/cache for every timed closure
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

    cursor_time = results["queue cursor"][1]
    println("  speedups_vs_cursor")
    for (name, _) in METHODS
        println("    $name=$(cursor_time / results[name][1])x")
    end
end

println("\nindexed collective-rate stress")
indexed_reference = run_cursor(indexed_collective_graph())
for rep in 1:20
    candidate = run_polyester(indexed_collective_graph())
    assert_same_graph(indexed_reference, candidate)
end
println("  Polyester: 20/20 fresh-context exact; worker derive did not grow SQA name/range intern tables")
