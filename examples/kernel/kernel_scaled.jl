# Treatments-aware kernel lowering for scaled and evaluated systems (spec finding 3).
#
# The prototype `statevars` resolves conjugate partners with an all-Free
# `canon_key(adjoint(k))`, which is wrong-by-construction once `scale` records a Scaled
# treatment (permutation-symmetry reduction) or `evaluate` records Concrete sites. The
# fix resolves every drift average leaf through the SAME `MomentMap`/`_treatments`
# machinery the MTK path uses (`_state_registry` in src/mtk.jl): one resolution code
# path, two consumers. Validations here:
#   1. identity on an unscaled system (bit-identical IR tables vs the all-Free path)
#   2. the old all-Free path fails concretely on the scaled superradiant laser
#   3. du vs a substitution reference: scaled + evaluated superradiant laser and the
#      evaluated cavity antiresonance (1D and 2D array parameters)
#   4. trajectory agreement vs the MTK path (independent anchor for the resolution)

using QuantumCumulants
using QuantumCumulants: build_ctx, _treatments, MomentMap, match_moment, undo_average,
    QAdd, eachleaf, _param_name, _param_slots
using ModelingToolkitBase
using OrdinaryDiffEqLowOrderRK
using Symbolics
using SymbolicUtils
using Random
import QuantumCumulants.SecondQuantizedAlgebra as SQA
include(joinpath(@__DIR__, "kernel_lower.jl"))
include(joinpath(@__DIR__, "kernel_eval.jl"))

# ---- treatments-aware state resolution -----------------------------------------------

"""
Resolution of every drift average leaf through the system's recorded treatments, the
`_state_registry` pattern with the state INDEX as the `MomentMap` payload. Returns
`(vars, idx)` in the `_lower_ir` contract: `vars` are the distinct leaf forms as they
appear in the drifts (handed to `polynomial_coeffs`), `idx[leaf]` the signed state index
(negative = conjugate side of the stored representative).
"""
function statevars_resolved(eqs)
    g = eqs.graph
    ctx = build_ctx(eqs)
    treatments = _treatments(eqs, ctx)
    ops = QAdd[(o = undo_average(s); o isa QAdd ? o : o * 1) for s in eqs.states]
    moments = MomentMap(ctx, treatments, ops, collect(Int32, 1:length(ops)))
    idx = Dict{Any, Int32}()
    vars = Any[]
    for nd in values(g.nodes), leaf in eachleaf(Symbolics.unwrap(nd.drift))
        haskey(idx, leaf) && continue
        op = undo_average(leaf)
        r = match_moment(moments, op isa QAdd ? op : op * 1)
        r === nothing && throw(UnresolvedMomentError(leaf))
        i, same = r
        idx[leaf] = same ? i : Int32(-i)
        push!(vars, leaf)
    end
    return vars, idx
end

lower_resolved(eqs) =
    _lower_ir(eqs.graph, statevars_resolved(eqs)..., Symbolics.unwrap(eqs.iv))

# ---- array-aware parameter values ----------------------------------------------------

_pname(p) = SymbolicUtils.iscall(p) && SymbolicUtils.operation(p) === getindex ?
    Base.nameof(SymbolicUtils.arguments(p)[1]) : _param_name(p)
function _pslots(p)
    if SymbolicUtils.iscall(p) && SymbolicUtils.operation(p) === getindex
        return Int[
            a isa Number ? Int(a) : Int(SymbolicUtils.unwrap_const(a))
                for a in SymbolicUtils.arguments(p)[2:end]
        ]
    end
    return _param_slots(p)
end

"""
Substitution dict for `coefficient_values` from a `parameter_map(eqs, ...)` result:
scalar entries pass through keyed by their unwrapped sym; each discovered kernel
parameter that is an array access (`g[1]`, `Γ[2,1]`, or a callable indexed variable)
is matched by (name, concrete slots) against the array values.
"""
function kernel_pdict(ir, pmap)
    pd = Dict{Any, Any}()
    arrs = Dict{Symbol, Any}()
    named = Dict{Any, Any}()
    for (k, v) in pmap
        ku = Symbolics.unwrap(k)
        if v isa AbstractArray
            arrs[_pname(ku)] = v
        else
            pd[ku] = v
            n = _pname(ku)
            n === nothing || (named[(n, _pslots(ku))] = v)
        end
    end
    unmatched = Any[]
    for p in ir.params
        haskey(pd, p) && continue
        name = _pname(p)
        slots = _pslots(p)
        if name !== nothing && haskey(arrs, name) && slots !== nothing
            pd[p] = arrs[name][slots...]
        elseif name !== nothing && haskey(named, (name, slots))
            pd[p] = named[(name, slots)]
        elseif name !== nothing && haskey(arrs, name)
            pd[p] = arrs[name]
        else
            push!(unmatched, p)
        end
    end
    isempty(unmatched) || error("unmatched kernel parameters: $(unmatched)")
    return pd
end

# ---- validations ---------------------------------------------------------------------

"""du of the resolved kernel vs numeric substitution into the symbolic drifts."""
function check_du_resolved(eqs, pmap; npoints = 3, seed = 1)
    ir = lower_resolved(eqs)
    pd = kernel_pdict(ir, pmap)
    k = MomentKernel(ir, coefficient_values(ir, pd))
    vars, idx = statevars_resolved(eqs)
    drifts = [Symbolics.unwrap(nd.drift) for nd in values(eqs.graph.nodes)]
    n = length(drifts)
    rng = MersenneTwister(seed)
    maxrel = 0.0
    for _ in 1:npoints
        u = randn(rng, ComplexF64, n)
        du = similar(u)
        k(du, u, nothing, 0.0)
        subs = Dict{Any, Any}(pd)
        for v in vars
            j = idx[v]
            subs[v] = j > 0 ? u[j] : conj(u[-j])
        end
        for d in drifts, v in Symbolics.get_variables(d)
            uu = SymbolicUtils.unwrap(v)
            SymbolicUtils.issym(uu) && SymbolicUtils.nameof(uu) === :im && (subs[uu] = im)
        end
        for i in 1:n
            ref = ComplexF64(SymbolicUtils.unwrap_const(Symbolics.substitute(drifts[i], subs)))
            maxrel = max(maxrel, abs(du[i] - ref) / max(abs(ref), 1e-12))
        end
    end
    return maxrel
end

"""Max deviation over all states and save points between the kernel and MTK trajectories."""
function traj_agreement(eqs, pmap, tspan; saveat)
    ir = lower_resolved(eqs)
    k = MomentKernel(ir, coefficient_values(ir, kernel_pdict(ir, pmap)))
    u0 = zeros(ComplexF64, ir.nstates)
    sol_k = solve(
        ODEProblem(ODEFunction(k), u0, tspan), RK4();
        saveat, abstol = 1.0e-10, reltol = 1.0e-10,
    )
    sys = mtkcompile(System(eqs; name = :sys))
    prob_m = ODEProblem(sys, merge(initial_values(eqs, u0), Dict(pmap)), tspan)
    sol_m = solve(prob_m, RK4(); saveat, abstol = 1.0e-10, reltol = 1.0e-10)
    ts = sol_k.t
    maxdev = 0.0
    for (i, st) in enumerate(eqs.states)
        ref = get_solution(sol_m, SymbolicUtils.unwrap(st), eqs).(ts)
        mine = [sol_k.u[j][i] for j in eachindex(ts)]
        maxdev = max(maxdev, maximum(abs.(ref .- mine)))
    end
    return maxdev
end

# ---- 1. identity on an unscaled system ------------------------------------------------

let
    Nc = 3
    h = ⊗([PauliSpace(Symbol(:spin, i)) for i in 1:Nc]...)
    sz(i) = Pauli(h, :σ, 3, i); sx(i) = Pauli(h, :σ, 1, i); sy(i) = Pauli(h, :σ, 2, i)
    sm(i) = (sx(i) - 1im * sy(i)) / 2
    @variables Jc hc γc
    Hc = -Jc * sum(sz(i) * sz(i + 1) for i in 1:(Nc - 1)) - hc * sum(sx(i) for i in 1:Nc)
    eqs = meanfield([sz(i) for i in 1:Nc], Hc, [sm(i) for i in 1:Nc]; rates = [γc for i in 1:Nc], order = 2)
    complete!(eqs)
    a = lower(eqs); b = lower_resolved(eqs)
    same = a.parent == b.parent && a.leaf == b.leaf && a.coo_i == b.coo_i &&
        a.coo_j == b.coo_j && a.coo_c == b.coo_c &&
        all(isequal(x, y) for (x, y) in zip(a.coeffs, b.coeffs)) &&
        all(isequal(x, y) for (x, y) in zip(a.params, b.params))
    println("SCALED identity on unscaled system (IR tables equal): ", same)
end

# ---- superradiant laser model (order 2, phase-invariant completion) --------------------

hc = FockSpace(:cavity)
ha = NLevelSpace(:atom, 2)
h = hc ⊗ ha
@qnumbers a::Destroy(h)
σ(α, β, k) = IndexedOperator(Transition(h, :σ, α, β), k)
@variables N Δ κ Γ R ν
g(k) = IndexedVariable(:g, k)
i = Index(h, :i, N, ha)
j = Index(h, :j, N, ha)
H = -Δ * a'a + Σ(g(i) * (a' * σ(1, 2, i) + a * σ(2, 1, i)), i)
J = [a, σ(1, 2, i), σ(2, 1, i), σ(2, 2, i)]
rates = [κ, Γ, R, ν]

φ(x) = 0
function φ(op::SQA.Op)
    SQA.is_destroy(op) && return -1
    SQA.is_create(op) && return 1
    SQA.is_transition(op) && return op.l1 - op.l2
    return 0
end
function φ(q::SQA.QAdd)
    for (term, _) in q.arguments
        return sum(φ(op) for op in term.ops)
    end
    return 0
end
function φ(avg)
    SQA.is_average(avg) || return 0
    return φ(SQA.undo_average(avg))
end
phase_invariant(x) = iszero(φ(x))

eqs_sr = meanfield([a' * a, σ(2, 2, j)], H, J; rates = rates, order = 2)
eqs_src = complete(eqs_sr; filter_func = phase_invariant)

pvals = Dict(Δ => 2.0, κ => 8.0, Γ => 1.0, R => 2.0, ν => 1.0)

# ---- 2 + 3 + 4a. scaled superradiant laser --------------------------------------------

"""Concrete old-path probe: does the all-Free lowering error, and if not, does it build
the same tables as the treatments-aware one?"""
function probe_old_path(label, eqs)
    r = try
        ir_o = lower(eqs)
        ir_n = lower_resolved(eqs)
        same = ir_o.parent == ir_n.parent && ir_o.leaf == ir_n.leaf &&
            ir_o.coo_i == ir_n.coo_i && ir_o.coo_j == ir_n.coo_j && ir_o.coo_c == ir_n.coo_c
        same ? "lowered, tables identical (coincides here)" :
            "lowered but tables DIFFER (silent wrong RHS)"
    catch e
        "threw $(typeof(e))"
    end
    println("SCALED old all-Free path on $label: ", r)
end

eqs_sc = scale(eqs_src)
probe_old_path("scaled superradiant", eqs_sc)

pmap_sc = parameter_map(eqs_sc, merge(pvals, Dict(N => 20.0, g(i) => 1.5)))
ir_sc = lower_resolved(eqs_sc)
println("SCALED scaled superradiant: nstates=$(ir_sc.nstates) params=$(ir_sc.params)")
e_sc = check_du_resolved(eqs_sc, pmap_sc)
println("SCALED scaled superradiant du vs substitution: max rel err = ", e_sc)
d_sc = traj_agreement(eqs_sc, pmap_sc, (0.0, 2.0); saveat = 0.2)
println("SCALED scaled superradiant trajectory vs MTK: max dev = ", d_sc)

# ---- 3 + 4b. evaluated superradiant laser ---------------------------------------------

Nev = 4
eqs_ev = evaluate(eqs_src; limits = (N => Nev))
probe_old_path("evaluated superradiant", eqs_ev)
pmap_ev = parameter_map(eqs_ev, merge(pvals, Dict(g(i) => [1.5 - 0.1k for k in 1:Nev])))
ir_ev = lower_resolved(eqs_ev)
println("SCALED evaluated superradiant: nstates=$(ir_ev.nstates) nparams=$(length(ir_ev.params))")
e_ev = check_du_resolved(eqs_ev, pmap_ev)
println("SCALED evaluated superradiant du vs substitution: max rel err = ", e_ev)
d_ev = traj_agreement(eqs_ev, pmap_ev, (0.0, 2.0); saveat = 0.2)
println("SCALED evaluated superradiant trajectory vs MTK: max dev = ", d_ev)

# ---- 3 + 4c. evaluated cavity antiresonance (2D array parameters) ----------------------

let
    hc2 = FockSpace(:cavity)
    ha2 = NLevelSpace(:atom2, 2)
    h2 = hc2 ⊗ ha2
    @variables N2 Δc η Δa κ2
    g2(k) = IndexedVariable(:g2, k)
    Γ2(k, l) = DoubleIndexedVariable(:Γ2, k, l)
    Ω2(k, l) = DoubleIndexedVariable(:Ω2, k, l; identical = false)
    i2 = Index(h2, :i, N2, ha2)
    j2 = Index(h2, :j, N2, ha2)
    @qnumbers b::Destroy(h2)
    s(x, y, k) = IndexedOperator(Transition(h2, :σ, x, y), k)
    Hc2 = Δc * b'b + η * (b' + b)
    Ha2 = Δa * Σ(s(2, 2, i2), i2) + Σ(Σ(Ω2(i2, j2) * s(2, 1, i2) * s(1, 2, j2), j2, [i2]), i2)
    Hi2 = Σ(g2(i2) * (b' * s(1, 2, i2) + b * s(2, 1, i2)), i2)
    eqs_a = meanfield(b, Hc2 + Ha2 + Hi2, [b, s(1, 2, i2)]; rates = [κ2, Γ2(i2, j2)], order = 1)
    complete!(eqs_a)
    Na = 2
    eqs_ae = evaluate(eqs_a; limits = (N2 => Na))
    probe_old_path("evaluated antiresonance", eqs_ae)
    Γm = [k == l ? 1.0 : 0.4 for k in 1:Na, l in 1:Na]
    Ωm = [k == l ? 0.0 : 1.1 for k in 1:Na, l in 1:Na]
    pmap_a = parameter_map(
        eqs_ae, Dict(
            Δc => 1.0, η => 0.2, Δa => 0.5, κ2 => 20.0,
            g2(i2) => [2.0, -2.0], Γ2(i2, j2) => Γm, Ω2(i2, j2) => Ωm,
        )
    )
    ir_a = lower_resolved(eqs_ae)
    println("SCALED evaluated antiresonance: nstates=$(ir_a.nstates) params=$(ir_a.params)")
    e_a = check_du_resolved(eqs_ae, pmap_a)
    println("SCALED evaluated antiresonance du vs substitution: max rel err = ", e_a)
    d_a = traj_agreement(eqs_ae, pmap_a, (0.0, 5.0); saveat = 0.5)
    println("SCALED evaluated antiresonance trajectory vs MTK: max dev = ", d_a)
end
