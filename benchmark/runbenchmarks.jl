using BenchmarkTools
using QuantumCumulants

const SUITE = BenchmarkGroup()

BenchmarkTools.tune!(SUITE)
results = BenchmarkTools.run(SUITE; verbose = true)
display(median(results))

BenchmarkTools.save("benchmarks_output.json", median(results))
