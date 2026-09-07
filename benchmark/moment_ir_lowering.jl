# Stage-separated benchmark for the evaluator-independent MomentIR compiler.
#
# Run from the repository root in a fresh Julia process:
#   QC_BENCH_ORDER=3 julia --project=benchmark benchmark/moment_ir_lowering.jl
#   QC_BENCH_ORDER=4 julia --project=benchmark benchmark/moment_ir_lowering.jl

using QuantumCumulants
using Symbolics: @variables

const QC = QuantumCumulants
const N = 6

function ising_model(order)
    h = ⊗([PauliSpace(Symbol(:spin, i)) for i in 1:N]...)
    σx(i) = Pauli(h, :σ, 1, i)
    σy(i) = Pauli(h, :σ, 2, i)
    σz(i) = Pauli(h, :σ, 3, i)
    σm(i) = (σx(i) - 1im * σy(i)) / 2
    @variables J hx γ
    H = -J * sum(σz(i) * σz(i + 1) for i in 1:(N - 1)) - hx * sum(σx(i) for i in 1:N)
    return meanfield(
        [σz(i) for i in 1:N],
        H,
        [σm(i) for i in 1:N];
        rates = fill(γ, N),
        order,
    )
end

function timed(f)
    GC.gc()
    result = @timed f()
    return result.value, result.time, result.bytes
end

function report(order)
    open, meanfield_time, meanfield_bytes = timed(() -> ising_model(order))
    eqs, complete_time, complete_bytes = timed(() -> complete(open))
    ir, lowering_time, lowering_bytes = timed(() -> QC._lower_moment_ir(eqs))

    println("MomentIR lowering order=$order equations=$(length(eqs.states))")
    println("  meanfield time=$(meanfield_time)s allocations=$(meanfield_bytes)")
    println("  complete time=$(complete_time)s allocations=$(complete_bytes)")
    println("  lowering time=$(lowering_time)s allocations=$(lowering_bytes)")
    println("  states=$(length(ir.states)) monomials=$(length(ir.parent)) terms=$(length(ir.monomial))")
    println("  IR bytes=$(Base.summarysize(ir))")
    return nothing
end

report(parse(Int, get(ENV, "QC_BENCH_ORDER", "3")))
