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
    treatments::Dict{Int, SubspaceTreatment}     # per-subspace treatment; all-Free after seed/closure!
end

"""A graph is "quotiented" once any subspace has been reduced to a `Scaled` treatment."""
quotiented(g::MomentGraph) = any(==(Scaled), values(g.treatments))

"""
The key function matching the graph's current level: `scaled_key` once any subspace is
`Scaled`, otherwise `canon_key`. Kept for callers that want the treatment-matched key;
`closure!` keys with `canon_key` directly.
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
function closure!(
        g::MomentGraph; filter = _alltrue, get_adjoints::Bool = true,
        max_iter::Int = 100_000,
    )
    ctx = g.ctx
    seen = Set(keys(g.nodes))
    pending = collect(keys(g.nodes))
    iters = 0
    while !isempty(pending)
        # `max_iter` is a runaway backstop, NOT a closure limiter. Hitting it means the
        # hierarchy did not close; ERROR rather than silently return a truncated
        # (non-closed) system, which the numerical-system build would mask (it would
        # look closed but leak/drop moments).
        iters >= max_iter && error(
            "closure! did not close the hierarchy within $max_iter iterations " *
                "($(length(pending)) moments still pending); the system may not close.",
        )
        iters += 1
        nd = g.nodes[popfirst!(pending)]
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
            g.nodes[k] = derive(k, g.sys, ctx)
            push!(seen, k)
            push!(pending, k)
            if get_adjoints && kc != k
                # Track the conjugate partner as its own moment.
                g.nodes[kc] = derive(kc, g.sys, ctx)
                push!(seen, kc)
                push!(pending, kc)
            else
                push!(seen, kc)
            end
        end
    end
    return g
end
_alltrue(_) = true

"""
Rebuild a `MomentGraph` from an existing equations struct: every current state becomes a
node, re-derived on its canon key, ready for `closure!` to extend.
"""
function _graph_from_eqs(eqs::AbstractMeanFieldEquations; mix_choice = maximum)
    ctx = build_ctx(eqs)
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
    return MomentGraph(nodes, sys, ctx, all_free_treatments(ctx))
end

"""
Assemble the graph's moments into the array-backed `MeanFieldEquations` /
`NoiseMeanFieldEquations` struct. Each moment is a canon-rep `QAdd`; LHS and drift share
that rep so the system is self-consistent.
"""
function assemble_equations(g::MomentGraph)
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
            collect(sys.rates), sys.iv, sys.order, sys.direction;
            treatments = copy(g.treatments),
        )
    else
        noise_eqs = Symbolics.Equation[average(k) ~ g.nodes[k].noise for k in ks]
        # Nodes don't store op_noise; rebuild it from the noise builder so the struct's
        # operator_noise_equations column is populated.
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
            treatments = copy(g.treatments),
        )
    end
end
