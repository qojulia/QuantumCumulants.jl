using QuantumCumulants
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

# Same canonical key calculation as `_treatment_key`, but with the treatment
# fingerprint supplied by the caller so it is not reconstructed for every probe.
function treatment_key_pre_fp(
        op::SQA.QAdd,
        ctx::QC.CanonCtx,
        treatments::Dict{Int, QC.SubspaceTreatment},
        fp::QC.TreatmentFP,
    )
    return get!(ctx.cache.key, (op, fp)) do
        scaled = Set{Int}()
        concrete = Set{Int}()
        for (sp, t) in treatments
            t == QC.Scaled && push!(scaled, sp)
            t == QC.Concrete && push!(concrete, sp)
        end
        rename_spaces = Set{Int}()
        for sp in ctx.symmetric
            sp in concrete || push!(rename_spaces, sp)
        end

        base = QC._reorder_commuting(
            QC._drop_scope_non_equal(SQA.QAdd(op.arguments, SQA.Index[])),
        )
        base = QC._relabel_spaces(base, ctx, rename_spaces)
        key = QC._drop_all_non_equal(SQA.QAdd(base.arguments, SQA.Index[]))
        isempty(scaled) && return key
        return QC.symmetric_min(key, ctx, scaled)
    end
end
treatment_key_pre_fp(op, ::QC.CanonCtx, treatments, fp) = op

# Closure always uses the all-Free treatment. In that specialization `scaled` and
# `concrete` are empty and every symmetric subspace is relabelled. Supply both the
# fingerprint and rename set once for the whole construction pass.
function treatment_key_free_plan(
        op::SQA.QAdd,
        ctx::QC.CanonCtx,
        fp::QC.TreatmentFP,
        rename_spaces::Set{Int},
    )
    return get!(ctx.cache.key, (op, fp)) do
        base = QC._reorder_commuting(
            QC._drop_scope_non_equal(SQA.QAdd(op.arguments, SQA.Index[])),
        )
        base = QC._relabel_spaces(base, ctx, rename_spaces)
        return QC._drop_all_non_equal(SQA.QAdd(base.arguments, SQA.Index[]))
    end
end
treatment_key_free_plan(op, ::QC.CanonCtx, fp, rename_spaces) = op

function close_custom(
        g::QC.MomentGraph,
        mode::Symbol;
        filter = QC._alltrue,
        get_adjoints::Bool = false,
        foldable = QC._alltrue,
        max_iter::Int = 100_000,
    )
    ctx = g.ctx
    nodes = copy(g.nodes)
    seen = Set(keys(nodes))
    pending = collect(keys(nodes))
    cursor = 1

    # Preserve closure semantics: this pass is always all-Free, independent of any
    # treatment map stored on a transformed graph.
    free_treatments = mode === :baseline ? nothing : QC.all_free_treatments(ctx)
    fp = mode in (:fp_once, :free_plan) ? QC.treatment_fp(free_treatments) : nothing
    rename_spaces = mode === :free_plan ? copy(ctx.symmetric) : nothing

    keyfn = if mode === :baseline
        op -> QC.canon_key(op, ctx)
    elseif mode === :map_once
        op -> op isa SQA.QAdd ? QC._treatment_key(op, ctx, free_treatments) : op
    elseif mode === :fp_once
        op -> treatment_key_pre_fp(op, ctx, free_treatments, fp)
    elseif mode === :free_plan
        op -> treatment_key_free_plan(op, ctx, fp, rename_spaces)
    else
        error("unknown mode $mode")
    end

    while cursor <= length(pending)
        cursor > max_iter && error(
            "closure did not close the hierarchy within $max_iter iterations",
        )
        nd = nodes[pending[cursor]]
        cursor += 1
        for leaf in QC._drift_leaves(nd)
            op = SQA.undo_average(leaf)
            k = keyfn(op)
            k in seen && continue
            kc = keyfn(adjoint(op))
            if kc in seen && foldable(op)
                push!(seen, k)
                continue
            end
            filter(SQA.average(k)) || continue

            nodes[k] = QC.derive(k, g.sys, ctx)
            push!(seen, k)
            push!(pending, k)

            if kc != k && !(kc in seen)
                if get_adjoints
                    nodes[kc] = QC.derive(kc, g.sys, ctx)
                    push!(seen, kc)
                    push!(pending, kc)
                elseif foldable(op)
                    push!(seen, kc)
                end
            end
        end
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

run_baseline(g) = QC.closure(g; get_adjoints = false)
run_map_once(g) = close_custom(g, :map_once; get_adjoints = false)
run_fp_once(g) = close_custom(g, :fp_once; get_adjoints = false)
run_free_plan(g) = close_custom(g, :free_plan; get_adjoints = false)

const METHODS = (
    ("baseline cursor", run_baseline),
    ("free map once", run_map_once),
    ("free map + fp once", run_fp_once),
    ("all-Free canonical plan", run_free_plan),
)

println("samples per method: ", NSAMPLES)
println("warming on fresh order-2 graphs")
warm_ref = run_baseline(ising_graph(2))
for (name, method) in METHODS
    candidate = method(ising_graph(2))
    assert_same_graph(warm_ref, candidate)
    println("  warm $name: exact")
end

for order in (3, 4)
    println("\norder $order")
    reference = run_baseline(ising_graph(order))
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

    baseline_time, baseline_bytes = results["baseline cursor"]
    println("  relative_to_baseline")
    for (name, _) in METHODS
        t, b = results[name]
        println("    $name: speedup=$(baseline_time / t)x bytes_ratio=$(b / baseline_bytes)")
    end
end
