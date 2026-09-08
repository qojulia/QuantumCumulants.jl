using QuantumCumulants
using Statistics: median
using Symbolics: @variables

const QC = QuantumCumulants
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

function assert_same_graph(reference, candidate)
    rkeys = collect(keys(reference.nodes))
    ckeys = collect(keys(candidate.nodes))
    @assert length(rkeys) == length(ckeys)
    for i in eachindex(rkeys)
        @assert isequal(rkeys[i], ckeys[i])
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

run_serial(g) = QC._closure(g, QC._derive_frontier_serial; get_adjoints = false)
run_production(g) = QC.closure(g; get_adjoints = false)

println("Julia threads: ", Threads.nthreads())
println("samples per method: ", NSAMPLES)

for order in (2, 3, 4)
    println("\norder $order")
    reference = run_serial(ising_graph(order))
    println("  states=$(length(reference.nodes))")
    for (name, method) in (("serial", run_serial), ("production", run_production))
        times = Float64[]
        bytes = Float64[]
        for _ in 1:NSAMPLES
            candidate, elapsed, allocated = timed(() -> method(ising_graph(order)))
            assert_same_graph(reference, candidate)
            push!(times, elapsed)
            push!(bytes, allocated)
        end
        println("  $name: median_time=$(median(times))s median_bytes=$(median(bytes)) exact=true")
    end
end
