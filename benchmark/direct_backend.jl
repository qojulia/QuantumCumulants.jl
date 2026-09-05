# Reproducible stage-separated benchmark for the first direct structured backend.
#
# Run from the repository root with:
#   julia --project=benchmark benchmark/direct_backend.jl
#
# The N=6 transverse-field Ising hierarchy is the issue #294 fixture: order 3 is about
# 693 equations and order 4 is about 1900. This script deliberately keeps the stages separate
# so lowering, table construction, SciML construction, and runtime are not conflated.

using BenchmarkTools
using ModelingToolkitBase: mtkcompile
using OrdinaryDiffEqTsit5: Tsit5, solve
using QuantumCumulants
using SciMLBase: ODEProblem

const QC = QuantumCumulants
const N = 6
const TSPAN = (0.0, 0.01)

function ising_model(order)
    h = ⊗([PauliSpace(Symbol(:spin, i)) for i = 1:N]...)
    σx(i) = Pauli(h, :σ, 1, i)
    σy(i) = Pauli(h, :σ, 2, i)
    σz(i) = Pauli(h, :σ, 3, i)
    σm(i) = (σx(i) - 1im * σy(i)) / 2
    @variables J hx γ
    H = -J * sum(σz(i) * σz(i + 1) for i = 1:(N-1)) - hx * sum(σx(i) for i = 1:N)
    eqs = meanfield(
        [σz(i) for i = 1:N],
        H,
        [σm(i) for i = 1:N];
        rates = [γ for _ = 1:N],
        order,
    )
    return eqs, Dict(J => 1.0, hx => 1.0, γ => 0.2)
end

function timed(f)
    GC.gc()
    result = @timed f()
    return result.value, result.time, result.bytes
end

ns_to_seconds(ns) = ns / 1e9

function direct_evaluator(eqs, ir, ps)
    values = QC.kernel_pdict(ir.params, parameter_map(eqs, ps))
    cvals = QC.coefficient_values(ir, values)
    kernel = QC.MomentKernel(ir, cvals; parallel = false)
    return QC.KernelRHS(kernel, QC.KernelParameters(ir, kernel.Mt, values))
end

function report(order)
    open, ps = ising_model(order)
    _, meanfield_time, meanfield_bytes = timed(() -> ising_model(order)[1])
    eqs, complete_time, complete_bytes = timed(() -> complete(open))
    nstates = length(eqs.states)
    u0 = zeros(ComplexF64, nstates)

    ir, lowering_time, lowering_bytes = timed(() -> QC._lower_moment_ir(eqs))
    rhs, evaluator_time, evaluator_bytes = timed(() -> direct_evaluator(eqs, ir, ps))
    direct_problem, ode_time, ode_bytes = timed(
        () -> ODEProblem(eqs, u0, TSPAN, ps; backend = KernelBackend(parallel = false)),
    )

    du = similar(u0)
    _, first_rhs_time, first_rhs_bytes =
        timed(() -> direct_problem.f(du, direct_problem.u0, direct_problem.p, 0.0))
    direct_f, direct_p = direct_problem.f, direct_problem.p
    warm_rhs = @benchmark $direct_f($du, $u0, $direct_p, 0.0) samples = 100 evals = 1
    _, first_solve_time, first_solve_bytes =
        timed(() -> solve(direct_problem, Tsit5(); saveat = TSPAN[2]))
    update_values = Dict(first(keys(ps)) => 1.1)
    _, update_time, update_bytes =
        timed(() -> update_parameters!(direct_problem, update_values))

    println("direct order=$order equations=$nstates")
    println("  meanfield time=$(meanfield_time)s allocations=$(meanfield_bytes)")
    println("  complete time=$(complete_time)s allocations=$(complete_bytes)")
    println("  direct lowering time=$(lowering_time)s allocations=$(lowering_bytes)")
    println("  direct evaluator time=$(evaluator_time)s allocations=$(evaluator_bytes)")
    println("  direct ODEProblem time=$(ode_time)s allocations=$(ode_bytes)")
    println("  direct first RHS time=$(first_rhs_time)s allocations=$(first_rhs_bytes)")
    println(
        "  direct warm RHS median=$(ns_to_seconds(median(warm_rhs).time))s allocations=$(median(warm_rhs).allocs)",
    )
    println(
        "  direct first solve time=$(first_solve_time)s allocations=$(first_solve_bytes)",
    )
    println("  direct parameter update time=$(update_time)s allocations=$(update_bytes)")
    println("  direct evaluator bytes=$(Base.summarysize(rhs))")
    flush(stdout)
    get(ENV, "QC_BENCH_SKIP_MTK", "false") == "true" && return nothing

    sys, system_time, system_bytes =
        timed(() -> System(eqs; name = Symbol(:benchmark_, order)))
    compiled, compile_time, compile_bytes = timed(() -> mtkcompile(sys))
    mtk_values = parameter_map(compiled, merge(initial_values(eqs, u0), ps))
    mtk_problem, mtk_ode_time, mtk_ode_bytes =
        timed(() -> ODEProblem(compiled, mtk_values, TSPAN))
    mtk_du = similar(u0)
    _, mtk_rhs_time, mtk_rhs_bytes =
        timed(() -> mtk_problem.f(mtk_du, mtk_problem.u0, mtk_problem.p, 0.0))
    mtk_f, mtk_p = mtk_problem.f, mtk_problem.p
    mtk_warm_rhs =
        @benchmark $mtk_f($mtk_du, $mtk_problem.u0, $mtk_p, 0.0) samples = 100 evals = 1
    _, mtk_solve_time, mtk_solve_bytes =
        timed(() -> solve(mtk_problem, Tsit5(); saveat = TSPAN[2]))

    println("order=$order equations=$nstates")
    println("  meanfield time=$(meanfield_time)s allocations=$(meanfield_bytes)")
    println("  complete time=$(complete_time)s allocations=$(complete_bytes)")
    println("  direct lowering time=$(lowering_time)s allocations=$(lowering_bytes)")
    println("  direct evaluator time=$(evaluator_time)s allocations=$(evaluator_bytes)")
    println("  direct ODEProblem time=$(ode_time)s allocations=$(ode_bytes)")
    println("  direct first RHS time=$(first_rhs_time)s allocations=$(first_rhs_bytes)")
    println(
        "  direct warm RHS median=$(ns_to_seconds(median(warm_rhs).time))s allocations=$(median(warm_rhs).allocs)",
    )
    println(
        "  direct first solve time=$(first_solve_time)s allocations=$(first_solve_bytes)",
    )
    println("  direct parameter update time=$(update_time)s allocations=$(update_bytes)")
    println("  direct evaluator bytes=$(Base.summarysize(rhs))")
    println("  MTK System time=$(system_time)s allocations=$(system_bytes)")
    println("  MTK compile time=$(compile_time)s allocations=$(compile_bytes)")
    println("  MTK ODEProblem time=$(mtk_ode_time)s allocations=$(mtk_ode_bytes)")
    println("  MTK first RHS time=$(mtk_rhs_time)s allocations=$(mtk_rhs_bytes)")
    println(
        "  MTK warm RHS median=$(ns_to_seconds(median(mtk_warm_rhs).time))s allocations=$(median(mtk_warm_rhs).allocs)",
    )
    println("  MTK first solve time=$(mtk_solve_time)s allocations=$(mtk_solve_bytes)")
    println("  MTK evaluator bytes=$(Base.summarysize(mtk_problem.f))")
    return nothing
end

orders = haskey(ENV, "QC_BENCH_ORDER") ? (parse(Int, ENV["QC_BENCH_ORDER"]),) : (3, 4)
for order in orders
    report(order)
end
