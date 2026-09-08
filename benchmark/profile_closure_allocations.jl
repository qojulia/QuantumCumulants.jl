using QuantumCumulants
using Symbolics

const QC = QuantumCumulants
const SQA = QuantumCumulants.SecondQuantizedAlgebra

mutable struct StageTotals
    calls::Int
    time_ns::Dict{Symbol, Float64}
    bytes::Dict{Symbol, Int}
end
StageTotals() = StageTotals(0, Dict{Symbol, Float64}(), Dict{Symbol, Int}())

function record!(s::StageTotals, name::Symbol, timed)
    s.time_ns[name] = get(s.time_ns, name, 0.0) + timed.time * 1e9
    s.bytes[name] = get(s.bytes, name, 0) + timed.bytes
    return timed.value
end

function profiled_derive(op::QC.QAdd, sys, ctx::QC.CanonCtx, totals::StageTotals)
    totals.calls += 1

    t = @timed QC._operator_rhs(
        sys.direction, op, im * sys.hamiltonian,
        sys.jumps, sys.jumps_dagger, sys.rates,
    )
    op_drift = record!(totals, :operator_rhs, t)

    t = @timed SQA.expand_completeness(op_drift)
    op_drift = record!(totals, :expand_completeness, t)

    t = @timed QC._distinct_atom_indices([op])
    distinct = record!(totals, :distinct_indices, t)

    t = @timed QC._assume_distinct_atom_indices(op_drift, distinct)
    op_drift = record!(totals, :assume_distinct, t)

    t = @timed Symbolics.Num(QC.average_and_truncate(op_drift, sys.order, sys.mix_choice, ctx))
    drift = record!(totals, :average_truncate, t)

    t = @timed QC._reduce_ground_in_drift(drift)
    drift = record!(totals, :ground_reduction, t)

    if sys.efficiencies === nothing
        op_noise = nothing
        noise = nothing
    else
        t = @timed QC._noise_builder(sys.direction)(
            [op], sys.jumps, sys.jumps_dagger, sys.rates, sys.efficiencies,
        )
        noise_pair = record!(totals, :noise_builder, t)
        _, noise_eqs = noise_pair
        op_noise = nothing
        noise_rhs = noise_eqs[1].rhs
        t = @timed Symbolics.Num(
            sys.order === nothing ? noise_rhs :
            QC.cumulant_expansion(noise_rhs, sys.order; mix_choice = sys.mix_choice),
        )
        noise = record!(totals, :noise_truncate, t)
        t = @timed QC._reduce_ground_in_drift(noise)
        noise = record!(totals, :noise_ground_reduction, t)
    end

    t = @timed (QC.get_order(op), SQA.acts_on(op))
    order, aon = record!(totals, :metadata, t)
    return QC.NodeData(drift, op_drift, noise, op_noise, order, aon)
end

function ising_graph(order)
    N = 6
    h = ⊗([PauliSpace(Symbol(:spin, i)) for i in 1:N]...)
    σx(i) = Pauli(h, :σ, 1, i)
    σy(i) = Pauli(h, :σ, 2, i)
    σz(i) = Pauli(h, :σ, 3, i)
    σm(i) = (σx(i) - 1im * σy(i)) / 2
    @variables J hₓ γ
    H = -J * sum(σz(i) * σz(i + 1) for i in 1:(N - 1)) - hₓ * sum(σx(i) for i in 1:N)
    c_ops = [σm(i) for i in 1:N]
    eqs = meanfield(
        [σz(i) for i in 1:N], H, c_ops;
        rates = [γ for _ in 1:N], order = order,
    )
    return eqs.graph
end

function run_profile(order)
    g = ising_graph(order)
    totals = StageTotals()
    frontier = function(keys, sys, ctx)
        out = Vector{QC.NodeData}(undef, length(keys))
        @inbounds for i in eachindex(keys)
            out[i] = profiled_derive(keys[i], sys, ctx, totals)
        end
        return out
    end
    GC.gc()
    total = @timed QC._closure(g, frontier)
    completed = total.value
    println("PROFILE order=$order states=$(length(completed.nodes)) closure_calls=$(totals.calls)")
    println("TOTAL time_s=$(round(total.time; digits=6)) bytes=$(total.bytes)")
    sum_stage_bytes = sum(values(totals.bytes); init = 0)
    sum_stage_time = sum(values(totals.time_ns); init = 0.0)
    for name in sort!(collect(keys(totals.bytes)); by = x -> -totals.bytes[x])
        b = totals.bytes[name]
        tns = totals.time_ns[name]
        println(
            "STAGE name=$name bytes=$b byte_fraction=$(round(b / max(sum_stage_bytes, 1); digits=4)) " *
            "time_s=$(round(tns / 1e9; digits=6)) time_fraction=$(round(tns / max(sum_stage_time, 1); digits=4))",
        )
    end
    println("STAGE_SUM bytes=$sum_stage_bytes time_s=$(round(sum_stage_time / 1e9; digits=6))")
end

for order in (2, 3, 4)
    run_profile(order)
end
