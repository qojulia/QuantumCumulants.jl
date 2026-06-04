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
    quotiented::Bool                  # false after seed/closure! (canon_key); true after quotient! (orbit_key)
end

# The key function matching the graph's current level. Every leaf-to-node
# resolution (closure!, lower, jacobian, specialize) goes through this so levels
# never mix.
nodekey(g::MomentGraph) = g.quotiented ? orbit_key : canon_key

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
    return MomentGraph(nodes, sys, ctx, false)
end

# Leaves of a node's drift (and noise drift, when present) to scan for edges.
_drift_leaves(nd::NodeData) = nd.noise === nothing ? eachleaf(nd.drift) :
    vcat(eachleaf(nd.drift), eachleaf(nd.noise))

# frontier: edge-target node keys not yet in the graph (the public find_missing).
# Scans each node's drift/noise leaves, keys via canon_key, dedups; folds the
# conjugate partner when get_adjoints=false.
function frontier(g::MomentGraph; get_adjoints::Bool = true)
    ctx = g.ctx
    seen = Set(keys(g.nodes))
    missing = NodeKey[]
    for nd in values(g.nodes)
        for leaf in _drift_leaves(nd)
            op = undo_average(leaf)
            k = canon_key(op, ctx)
            (k in seen || k in missing) && continue
            kc = canon_key(adjoint(op), ctx)
            # A state and its conjugate are ONE physical unknown: if the
            # conjugate partner is already present (node or pending), this
            # leaf is covered. Unconditional in get_adjoints (master's
            # find_missing pre-marks every state's conjugate as seen).
            (kc in seen || kc in missing) && continue
            push!(missing, k)
            # get_adjoints=true tracks both reps of a genuinely-new pair.
            get_adjoints && kc != k && push!(missing, kc)
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
        max_iter::Int = 1000,
    )
    ctx = g.ctx
    seen = Set(keys(g.nodes))
    worklist = collect(keys(g.nodes))
    iters = 0
    while !isempty(worklist) && iters < max_iter
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
    return MomentGraph(nodes, sys, ctx, false)
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
    if sys.efficiencies === nothing
        return MeanFieldEquations(
            avg_eqs, operator_eqs, states, operators,
            sys.hamiltonian, collect(sys.jumps), collect(sys.jumps_dagger),
            collect(sys.rates), sys.iv, sys.order, sys.direction,
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
            sys.iv, sys.order, sys.direction,
        )
    end
end
