struct SystemSpec{H, Jt, Jdt, R, E, MC, D <: EvolutionDirection}
    hamiltonian::H
    jumps::Jt
    jumps_dagger::Jdt
    rates::R
    efficiencies::E # Nothing means deterministic
    iv::Symbolics.Num
    order::TruncOrder
    mix_choice::MC
    direction::D
end

const NodeKey = QAdd

struct MomentGraph{S <: SystemSpec}
    nodes::OrderedCollections.OrderedDict{NodeKey, NodeData}
    sys::S
    ctx::CanonCtx
    treatments::Dict{Int, SubspaceTreatment}
end

"""A graph is "quotiented" once any subspace has been reduced to a `Scaled` treatment."""
quotiented(g::MomentGraph) = any(==(Scaled), values(g.treatments))

"""True when `g` keeps at most one member of each conjugate pair (closed `get_adjoints=false`)."""
function conj_folded(g::MomentGraph)
    for k in keys(g.nodes)
        kc = canon_key(adjoint(k), g.ctx)
        kc != k && haskey(g.nodes, kc) && return false
    end
    return true
end

"""
The key function matching the graph's current level: `scaled_key` once any subspace is
`Scaled`, otherwise `canon_key`. Kept for callers that want the treatment-matched key;
`closure` keys with `canon_key` directly.
"""
nodekey(g::MomentGraph) = quotiented(g) ? scaled_key : canon_key

"""
Derive each seed op into a graph node, keyed by `canon_key` (the graph starts
un-quotiented). Each op is promoted to a `QAdd` first, exactly as `meanfield` does.
"""
function seed(ops::AbstractVector, sys::SystemSpec, ctx::CanonCtx)
    nodes = OrderedCollections.OrderedDict{NodeKey, NodeData}()
    for op in ops
        opqa = op isa QAdd ? op : op * 1
        k = canon_key(opqa, ctx)
        haskey(nodes, k) && continue
        nodes[k] = derive(k, sys, ctx)
    end
    return MomentGraph(nodes, sys, ctx, all_free_treatments(ctx))
end

"""Leaves of a node's drift (and noise drift, when present): the moments it couples to."""
_drift_leaves(nd::NodeData) = nd.noise === nothing ? eachleaf(nd.drift) :
    vcat(eachleaf(nd.drift), eachleaf(nd.noise))

"""
Close the cumulant hierarchy: starting from the seeded moments, repeatedly derive the
equation of motion for every moment appearing on a right-hand side until the set is
self-contained. Moments are keyed by `canon_key`, NOT `scaled_key`: symmetric reduction
here would collapse the unscaled count. `filter` zeroes unwanted moments (e.g. ancilla /
phase-invariant); `get_adjoints=false` keeps one moment per conjugate pair, the partner
recovered via `conj` when the numerical system is built.
"""
function closure(
        g::MomentGraph; filter = _alltrue, get_adjoints::Bool = true,
        max_iter::Int = 100_000,
    )
    ctx = g.ctx
    nodes = copy(g.nodes)   # shallow copy: NodeData values are shared (immutable), new moments appended here
    seen = Set(keys(nodes))
    pending = collect(keys(nodes))
    iters = 0
    while !isempty(pending)
        # `max_iter` is a runaway backstop, NOT a closure limiter. Hitting it means the
        # hierarchy did not close; ERROR rather than silently return a truncated
        # (non-closed) system, which the numerical-system build would mask (it would
        # look closed but leak/drop moments).
        iters >= max_iter && error(
            "closure did not close the hierarchy within $max_iter iterations " *
                "($(length(pending)) moments still pending); the system may not close.",
        )
        iters += 1
        nd = nodes[popfirst!(pending)]
        for leaf in _drift_leaves(nd)
            op = undo_average(leaf)
            k = canon_key(op, ctx)
            k in seen && continue
            kc = canon_key(adjoint(op), ctx)
            # A moment and its conjugate are ONE physical unknown: if the conjugate
            # partner is already tracked, this moment is covered (the numerical system
            # resolves ⟨X⟩ = conj⟨X†⟩). Unconditional in get_adjoints.
            if kc in seen
                push!(seen, k)
                continue
            end
            filter(average(k)) || continue
            # Genuinely new moment (neither rep seen yet).
            nodes[k] = derive(k, g.sys, ctx)
            push!(seen, k)
            push!(pending, k)
            if get_adjoints && kc != k
                # Track the conjugate partner as its own moment.
                nodes[kc] = derive(kc, g.sys, ctx)
                push!(seen, kc)
                push!(pending, kc)
            else
                push!(seen, kc)
            end
        end
    end
    return MomentGraph(nodes, g.sys, ctx, g.treatments)
end
_alltrue(_) = true

"""
Copy a `MomentGraph`'s mutable parts (the `nodes` and `treatments` dicts) while sharing the
immutable `sys` and `ctx`. Backs the non-mutating `_copy` of a wrapper.
"""
copy_graph(g::MomentGraph) =
    MomentGraph(copy(g.nodes), g.sys, g.ctx, copy(g.treatments))

"""
Rewrite each node's `drift` (and, when `noise = true` and the node carries one, its `noise`)
via `f(key, expr)`, returning a new graph with the same keys, `sys`, `ctx`, and `treatments`.
The single primitive behind the drift-rewrite transforms (simplify, modify, translate,
filter, arrayize). `f` receives the node key (the moment's operator) so key-aware rewrites
(`modify_equations`, `translate_W_to_Y`) can use the LHS operator.
"""
function map_drifts(g::MomentGraph, f; noise::Bool = true)
    nodes = OrderedCollections.OrderedDict{NodeKey, NodeData}()
    for (k, nd) in g.nodes
        new_drift = Symbolics.Num(SymbolicUtils.unwrap(f(k, nd.drift)))
        new_noise = (noise && nd.noise !== nothing) ?
            Symbolics.Num(SymbolicUtils.unwrap(f(k, nd.noise))) : nd.noise
        nodes[k] = NodeData(new_drift, nd.op_drift, new_noise, nd.op_noise, nd.order, nd.aon)
    end
    return MomentGraph(nodes, g.sys, g.ctx, g.treatments)
end

"""
Assemble a graph into the array-backed wrapper. Dispatches on `g.sys.efficiencies`; the
wrapper stores `g` as its source of truth and regenerates its view from it.
"""
assemble_equations(g::MomentGraph) =
    g.sys.efficiencies === nothing ? MeanfieldEquations(g) : NoiseMeanfieldEquations(g)
