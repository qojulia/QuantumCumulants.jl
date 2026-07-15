# Coupling-aware shard partitioner (prototype for issue #294)
#
# Reads the coupling structure of a completed cumulant hierarchy (QC's `MomentGraph`)
# and produces an equation ordering intended to keep tightly-coupled moments in the same
# contiguous block. Because `ShardedForm` chunks the RHS vector by position
# (`make_array(::ShardedForm, ...)` in Symbolics `build_function.jl`), permuting the
# equation order is all that is needed to turn a blind equal-chunk split into a
# coupling-aware one: no new codegen.
#
# Conclusion of the benchmark (see shard_bench.jl): for the N=6 order=3 transverse-field
# Ising chain the coupling-aware orderings do NOT beat a blind equal split. The heavily
# shared subexpressions are the low-order moments (e.g. ⟨σz⟩ on a single site is referenced
# by ~255 of 693 equations), which are global across every support, so no contiguous
# partition can localise them. Kept as a reproducible negative result.

using QuantumCumulants
using QuantumCumulants: MomentGraph, _drift_leaves, canon_key, undo_average

"""
    coupling_adjacency(g::MomentGraph)

Directed adjacency of the coupling graph: `adj[i]` is the set of node indices whose moment
appears on the RHS of equation `i`. Nodes are indexed by their order in `g.nodes`, i.e. the
same order as the assembled equation / state vector. Falls back to the conjugate key when a
drift leaf's own key is folded out of the graph.
"""
function coupling_adjacency(g::MomentGraph)
    ctx = g.ctx
    nodekeys = collect(keys(g.nodes))
    key2idx = Dict(k => i for (i, k) in enumerate(nodekeys))
    n = length(nodekeys)
    adj = [Set{Int}() for _ in 1:n]
    for (i, k) in enumerate(nodekeys)
        for leaf in _drift_leaves(g.nodes[k])
            op = undo_average(leaf)
            j = get(key2idx, canon_key(op, ctx), 0)
            j == 0 && (j = get(key2idx, canon_key(adjoint(op), ctx), 0))
            j != 0 && push!(adj[i], j)
        end
    end
    return adj, nodekeys
end

"""
    support_ordering(g)

QC-native ordering: cluster moments by support (the sites the operator acts on, `NodeData.aon`)
then by cumulant order. Moments touching the same site block, which is where shared RHS terms
come from in a local model, land in the same contiguous chunk.
"""
function support_ordering(g::MomentGraph)
    nodekeys = collect(keys(g.nodes))
    supp(i) = sort(collect(g.nodes[nodekeys[i]].aon))
    ordr(i) = g.nodes[nodekeys[i]].order
    return sort(1:length(nodekeys); by = i -> (supp(i), ordr(i), string(nodekeys[i])))
end

"""
    rcm_ordering(g)

Reverse Cuthill-McKee ordering of the (symmetrised) coupling graph: makes directly-coupled
moments contiguous, minimising bandwidth. This is the best case for contiguous positional
shards, so it bounds how much any coupling-aware ordering could help.
"""
function rcm_ordering(g::MomentGraph)
    adj, nodekeys = coupling_adjacency(g)
    n = length(adj)
    sym = [Set{Int}() for _ in 1:n]
    for i in 1:n, j in adj[i]
        push!(sym[i], j); push!(sym[j], i)
    end
    deg = [length(sym[i]) for i in 1:n]
    visited = falses(n); order = Int[]
    for s in sortperm(deg)
        visited[s] && continue
        q = [s]; visited[s] = true
        while !isempty(q)
            v = popfirst!(q); push!(order, v)
            for w in sort([w for w in sym[v] if !visited[w]]; by = w -> deg[w])
                visited[w] = true; push!(q, w)
            end
        end
    end
    return reverse(order)
end
