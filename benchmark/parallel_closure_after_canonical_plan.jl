using QuantumCumulants
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

function derive_polyester(keys, sys, ctx)
    out = Vector{QC.NodeData}(undef, length(keys))
    Polyester.@batch per=thread for i in eachindex(keys)
        @inbounds out[i] = QC.derive(keys[i], sys, ctx)
    end
    return out
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
    free_treatments = QC.all_free_treatments(ctx)
    free_fp = QC.treatment_fp(free_treatments)
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
                k = QC._treatment_key(op, ctx, free_treatments, free_fp)
                k in seen && continue
                kc = QC._treatment_key(adjoint(op), ctx, free_treatments, free_fp)
                if kc in seen && foldable(op)
                    push!(seen, k)
                    continue
                end
                filter(SQA.average(k)) || continue

                # Discovery, de-duplication and final insertion order remain serial.
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

        names_before = length(SQA._NAME_BY_ID)
        ranges_before = length(SQA._RANGE_BY_ID)
        derived = derive_frontier(newkeys, g.sys, ctx)
        if assert_intern_stable
            @assert length(SQA._NAME_BY_ID) == names_before
            @assert length(SQA._RANGE_BY_ID) == ranges_before
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

run_production(g) = QC.closure(g; get_adjoints = false)
run_frontier_serial(g) = closure_frontier(g, derive_serial; get_adjoints = false)
run_polyester(g) = closure_frontier(g, derive_polyester; get_adjoints = false)

const METHODS = (
    ("production serial", run_production),
    ("frontier serial", run_frontier_serial),
    ("Polyester frontier", run_polyester),
)

println("Julia threads: ", Threads.nthreads())
println("samples per method: ", NSAMPLES)
println("warming on fresh order-2 graphs")
warm_ref = run_production(ising_graph(2))
for (name, method) in METHODS
    candidate = method(ising_graph(2))
    assert_same_graph(warm_ref, candidate)
    println("  warm $name: exact")
end

for order in (2, 3, 4)
    println("\norder $order")
    reference = run_production(ising_graph(order))
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

    baseline_time, baseline_bytes = results["production serial"]
    println("  relative_to_production")
    for (name, _) in METHODS
        t, b = results[name]
        println("    $name: speedup=$(baseline_time / t)x bytes_ratio=$(b / baseline_bytes)")
    end
end

println("\nindexed collective-rate stress")
indexed_reference = run_production(indexed_collective_graph())
for _ in 1:20
    candidate = run_polyester(indexed_collective_graph())
    assert_same_graph(indexed_reference, candidate)
end
println("  Polyester: 20/20 fresh-context exact; worker derive did not grow SQA name/range intern tables")
