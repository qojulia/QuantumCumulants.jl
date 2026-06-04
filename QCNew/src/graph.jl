# The IR: the moment dependency graph (spec Section 3.4). Nodes are canonical
# QAdds, each annotated with its NodeData. Edges stay implicit (a node's
# out-edges are the keyed leaves of its drift). Phase 0 builds the types + seed;
# closure!/quotient!/specialize/lower/jacobian are Phases 3 to 6.

struct SystemSpec{H, Jt, Jdt, R, E, MC, D <: EvolutionDirection}
    hamiltonian::H
    jumps::Jt
    jumps_dagger::Jdt
    rates::R
    efficiencies::E                  # Nothing means deterministic
    iv::Symbolics.Num
    order::TruncOrder
    mix_choice::MC                   # parametric so dispatch stays static (not abstract `Function`)
    direction::D                     # parametric `D<:EvolutionDirection` (matches MeanFieldEquations)
end

const NodeKey = QAdd

struct MomentGraph{S <: SystemSpec}
    nodes::OrderedCollections.OrderedDict{NodeKey, NodeData}
    sys::S
    ctx::CanonCtx
    coords::Dict{Int, Coordinate}     # per-subspace coordinate; all-Free after seed/closure!
end

# Backward-compatible helper: a graph is "quotiented" once any subspace is Scaled.
quotiented(g::MomentGraph) = any(==(Scaled), values(g.coords))

# The key function matching the graph's current level. Every leaf-to-node
# resolution that needs a single key function reads the coordinate: Scaled
# present means orbit_key, otherwise canon_key. (closure! keys with canon_key
# directly; this is kept for callers that want the coordinate-matched key.)
nodekey(g::MomentGraph) = quotiented(g) ? orbit_key : canon_key

# Derive each seed op into a node, keyed by canon_key (the graph starts
# un-quotiented). Promote each op to a QAdd first, exactly as meanfield does.
function seed(ops::AbstractVector, sys::SystemSpec, ctx::CanonCtx)
    nodes = OrderedCollections.OrderedDict{NodeKey, NodeData}()
    for op in ops
        opqa = op isa QAdd ? op : op * 1
        k = canon_key(opqa, ctx)
        haskey(nodes, k) && continue
        nodes[k] = derive(k, sys, ctx)
    end
    return MomentGraph(nodes, sys, ctx, all_free_coords(ctx))
end

# Leaves of a node's drift (and noise drift, when present) to scan for edges.
_drift_leaves(nd::NodeData) = nd.noise === nothing ? eachleaf(nd.drift) :
    vcat(eachleaf(nd.drift), eachleaf(nd.noise))

# frontier: edge-target node keys not yet in the graph (the public find_missing).
# Scans each node's drift/noise leaves, keys via canon_key, dedups; folds the
# conjugate partner when get_adjoints=false.
function frontier(g::MomentGraph; get_adjoints::Bool = true)
    ctx = g.ctx
    coords = g.coords
    # Match in the system's coordinate (spec "Closure and the system coordinate"):
    # fold every node key AND every leaf through the SAME `canonical_rep(·; coords)`,
    # so a permutation-image / conjugate leaf folds to the rep its stored node was
    # keyed to. This is the SAME code path the codegen resolver uses, so
    # `find_missing == 0` means exactly "every leaf resolves at codegen".
    seen = Set{NodeKey}(canonical_rep(k, ctx; coords)[1] for k in keys(g.nodes))
    missing = NodeKey[]
    seen_missing = Set{NodeKey}()
    for nd in values(g.nodes)
        for leaf in _drift_leaves(nd)
            op = undo_average(leaf)
            op isa QAdd || continue
            rep, _ = canonical_rep(op, ctx; coords)
            (rep in seen || rep in seen_missing) && continue
            push!(missing, rep)
            push!(seen_missing, rep)
            if get_adjoints
                repc, _ = canonical_rep(adjoint(op), ctx; coords)
                if !isequal(repc, rep) && !(repc in seen) && !(repc in seen_missing)
                    push!(missing, repc)
                    push!(seen_missing, repc)
                end
            end
        end
    end
    return missing
end

# closure!: BFS to fixpoint over the implicit edge set. Each discovered key is
# derived on enqueue (the node key is canon_key, NOT orbit_key: symmetric
# reduction here would collapse the unscaled count). `filter` zeroes unwanted
# leaves (ancilla/phase-invariant); `get_adjoints=false` tracks one rep per
# conjugate pair (the partner resolved via conj at codegen). Spec section 3.4.
function closure!(
        g::MomentGraph; filter = _alltrue, get_adjoints::Bool = true,
        max_iter::Int = 100_000,
    )
    ctx = g.ctx
    seen = Set(keys(g.nodes))
    worklist = collect(keys(g.nodes))
    iters = 0
    while !isempty(worklist)
        # `max_iter` is a runaway backstop, NOT a closure limiter. Hitting it
        # means the BFS did not reach the fixpoint; ERROR rather than silently
        # return a truncated (non-closed) system, which downstream codegen would
        # mask (the system would look closed but leak/drop leaves).
        iters >= max_iter && error(
            "closure! did not reach the fixpoint within $max_iter iterations " *
            "($(length(worklist)) nodes still pending); the system may not close.",
        )
        iters += 1
        nd = g.nodes[popfirst!(worklist)]
        for leaf in _drift_leaves(nd)
            op = undo_average(leaf)
            k = canon_key(op, ctx)
            k in seen && continue
            kc = canon_key(adjoint(op), ctx)
            # A state and its conjugate are ONE physical unknown: if the
            # conjugate partner is already a node, this leaf is covered
            # (codegen resolves ⟨X⟩ = conj⟨X†⟩). Unconditional in
            # get_adjoints, matching master's conjugate dedup.
            if kc in seen
                push!(seen, k)
                continue
            end
            filter(average(k)) || continue
            # Genuinely new state (neither rep seen yet).
            g.nodes[k] = derive(k, g.sys, ctx)
            push!(seen, k)
            push!(worklist, k)
            if get_adjoints && kc != k
                # Track the conjugate partner as its own node.
                g.nodes[kc] = derive(kc, g.sys, ctx)
                push!(seen, kc)
                push!(worklist, kc)
            else
                push!(seen, kc)
            end
        end
    end
    return g
end
_alltrue(_) = true

# Rebuild a MomentGraph from an existing equations struct: every current state
# becomes a node (re-derived on its canon key), ready for closure! to extend.
function _graph_from_eqs(eqs::AbstractMeanFieldEquations; mix_choice = maximum)
    ctx = build_ctx(eqs.operators, eqs.hamiltonian, eqs.jumps, eqs.jumps_dagger)
    eff = eqs isa NoiseMeanFieldEquations ? eqs.efficiencies : nothing
    sys = SystemSpec(
        eqs.hamiltonian, eqs.jumps, eqs.jumps_dagger, eqs.rates, eff,
        eqs.iv, eqs.order, mix_choice, eqs.direction,
    )
    nodes = OrderedCollections.OrderedDict{NodeKey, NodeData}()
    for op in eqs.operators
        opqa = op isa QAdd ? op : op * 1
        k = canon_key(opqa, ctx)
        haskey(nodes, k) && continue
        nodes[k] = derive(k, sys, ctx)
    end
    return MomentGraph(nodes, sys, ctx, all_free_coords(ctx))
end

# Lower a graph's nodes back into the array-backed MeanFieldEquations /
# NoiseMeanFieldEquations struct (the frozen Phase-1 interface). Nodes are
# canon-rep QAdds; LHS and drift share that rep so the system is self-consistent.
function lower_to_eqs(g::MomentGraph)
    sys = g.sys
    ks = collect(keys(g.nodes))
    states = SymbolicUtils.BasicSymbolic[average(k) for k in ks]
    operators = QAdd[k for k in ks]
    operator_eqs = Symbolics.Equation[k ~ g.nodes[k].op_drift for k in ks]
    avg_eqs = Symbolics.Equation[average(k) ~ g.nodes[k].drift for k in ks]
    coords_int = Dict{Int, Int}(sp => Int(c) for (sp, c) in g.coords)
    if sys.efficiencies === nothing
        return MeanFieldEquations(
            avg_eqs, operator_eqs, states, operators,
            sys.hamiltonian, collect(sys.jumps), collect(sys.jumps_dagger),
            collect(sys.rates), sys.iv, sys.order, sys.direction;
            coords = coords_int,
        )
    else
        noise_eqs = Symbolics.Equation[average(k) ~ g.nodes[k].noise for k in ks]
        # op_noise is deferred (Phase 0); rebuild the operator noise eqs from the
        # builder so the struct's operator_noise_equations column is populated.
        op_noise = Symbolics.Equation[]
        for k in ks
            on, _ = _noise_builder(sys.direction)(
                [k], sys.jumps, sys.jumps_dagger,
                sys.rates, sys.efficiencies
            )
            push!(op_noise, on[1])
        end
        return NoiseMeanFieldEquations(
            avg_eqs, noise_eqs, operator_eqs, op_noise,
            states, operators, sys.hamiltonian,
            collect(sys.jumps), collect(sys.jumps_dagger),
            collect(sys.rates), collect(sys.efficiencies),
            sys.iv, sys.order, sys.direction;
            coords = coords_int,
        )
    end
end
