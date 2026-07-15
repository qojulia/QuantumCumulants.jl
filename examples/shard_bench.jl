# Static analysis for issue #294: why a coupling-aware shard partition does not beat a blind
# equal split. This file does NOT time compilation: first-call/first-solve compile time is
# shared across builds within one Julia session (and perturbed by Revise), so any two builds
# timed back-to-back are not comparable. Use shard_measure.jl for timing, one path per fresh
# session. See docs/superpowers/specs/2026-07-15-coupling-aware-shard-findings.md.
#
# Cold, isolated first-solve times measured via shard_measure.jl (order 3, N=6, 693 eqs,
# Symbolics 7.32 / MTKBase 1.52):
#
#   ODEProblem default                 60.5 s
#   ODEProblem eval_expression=true    60.5 s
#   legacy build_function serial       72.7 s
#   legacy ShardedForm(10,4) blind     28.1 s   <- sharding is the win (2.6x vs serial)
#   legacy + support-locality reorder  28.3 s   <- coupling-aware: no improvement
#   legacy + RCM reorder               28.5 s   <- coupling-aware: no improvement
#
# This file shows WHY the reorderings can't help: the heavily shared subexpressions are global.

using QuantumCumulants
using ModelingToolkitBase
using ModelingToolkitBase: unknowns, equations, parameters, get_iv
using Symbolics
include(joinpath(@__DIR__, "shard_partition.jl"))

const N = 6
h = ⊗([PauliSpace(Symbol(:spin, i)) for i in 1:N]...)
σx(i) = Pauli(h, :σ, 1, i); σy(i) = Pauli(h, :σ, 2, i); σz(i) = Pauli(h, :σ, 3, i)
σm(i) = (σx(i) - 1im * σy(i)) / 2
@variables J hx γ
H = -J * sum(σz(i) * σz(i + 1) for i in 1:(N - 1)) - hx * sum(σx(i) for i in 1:N)
eqs = meanfield([σz(i) for i in 1:N], H, [σm(i) for i in 1:N]; rates = [γ for i in 1:N], order = 3)
complete!(eqs)
sys = mtkcompile(System(eqs; name = :sys))
dvs = unknowns(sys); ps = parameters(sys); tiv = get_iv(sys)
rhss = [eq.rhs for eq in equations(sys)]

expr_nodes(ex) = (n = 0; walk(x) = (n += 1; x isa Expr && foreach(walk, x.args)); walk(ex); n)

# Code duplication under CSE: a partition duplicates any subexpression shared across groups.
# Sum of per-group CSE'd code size; serial (one group) is the no-duplication floor. If a
# coupling-aware ordering kept shared terms together, its excess over serial would be smaller
# than the blind split's. It is not.
println("=== code duplication, cse=true (serial = no-duplication floor) ===")
serial_code = expr_nodes(build_function(rhss, dvs, ps, tiv; expression = Val{true}, cse = true)[2])
println("serial (one fused function): ", serial_code, " nodes\n")
function dup(perm, ncalls)
    per = ceil(Int, length(perm) / ncalls); tot = 0
    for chunk in Iterators.partition(perm, per)
        tot += expr_nodes(build_function(rhss[collect(chunk)], dvs, ps, tiv; expression = Val{true}, cse = true)[2])
    end
    return tot
end
perm_supp = support_ordering(eqs.graph)
perm_rcm = rcm_ordering(eqs.graph)
for nc in (4, 8), (name, perm) in (("blind", collect(1:length(rhss))), ("support", perm_supp), ("rcm", perm_rcm))
    d = dup(perm, nc)
    println(rpad("$name  nc=$nc", 16), ": +", round(100 * (d - serial_code) / serial_code, digits = 1), "% over serial")
end

# The mechanism: the most-shared subexpressions are the low-order moments, each on the RHS of
# a large fraction of all equations, spread across every support. No contiguous partition can
# localise a term referenced by hundreds of equations, so every split duplicates it regardless
# of ordering.
adj, nodekeys = coupling_adjacency(eqs.graph)
refcount = zeros(Int, length(nodekeys))
for i in eachindex(adj), j in adj[i]
    refcount[j] += 1
end
println("\n=== most-referenced moments (why coupling-awareness can't localise them) ===")
for j in sortperm(refcount; rev = true)[1:6]
    println(
        "  ", rpad(string(nodekeys[j]), 12), " on RHS of ", refcount[j], "/", length(rhss),
        " eqs, support=", eqs.graph.nodes[nodekeys[j]].aon
    )
end
