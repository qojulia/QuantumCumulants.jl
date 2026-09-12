# Reproducible stage-separated benchmark for the direct structured backend.
#
# Run from the repository root in a fresh Julia process:
#   QC_BENCH_ORDER=3 julia --project=benchmark benchmark/direct_backend.jl
#   QC_BENCH_ORDER=4 julia --project=benchmark benchmark/direct_backend.jl
#   QC_BENCH_MODE=direct QC_BENCH_ORDER=4 julia --project=benchmark benchmark/direct_backend.jl
#   QC_BENCH_MODE=mtk QC_BENCH_ORDER=4 julia --project=benchmark benchmark/direct_backend.jl
#
# The N=6 transverse-field Ising hierarchy is the issue #294 fixture. Keep stages separate so
# lowering, evaluator construction, SciML construction, first solve, and warm RHS cost are not
# conflated. Run each mode/order in a fresh Julia process for cold construct-to-solve numbers.

using BenchmarkTools
using ModelingToolkitBase: mtkcompile
using OrdinaryDiffEqTsit5: Tsit5, solve
using QuantumCumulants

const QC = QuantumCumulants
const N = 6
const TSPAN = (0.0, 0.01)

function ising_model(order)
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
    return eqs, Dict(J => 1.0, hx => 1.0, γ => 0.2)
end

function timed(f)
    GC.gc()
    result = @timed f()
    return result.value, result.time, result.bytes
end

ns_to_seconds(ns) = ns / 1.0e9

function direct_evaluator(eqs, ir, ps)
    kernel = QC.MomentKernel(ir, ComplexF64)
    plan = QC.KernelParameterPlan(ir)
    p = QC.KernelParameters(plan, parameter_map(eqs, ps))
    return QC.KernelRHS(kernel, plan), p
end

function report(order)
    open, ps = ising_model(order)
    _, meanfield_time, meanfield_bytes = timed(() -> ising_model(order)[1])
    eqs, complete_time, complete_bytes = timed(() -> complete(open))
    nstates = length(eqs.states)
    u0 = zeros(ComplexF64, nstates)

    mode = get(ENV, "QC_BENCH_MODE", "both")
    mode in ("direct", "mtk", "both") ||
        throw(ArgumentError("QC_BENCH_MODE must be `direct`, `mtk`, or `both`; got $mode"))

    println("  meanfield time=$(meanfield_time)s allocations=$(meanfield_bytes)")
    println("  complete time=$(complete_time)s allocations=$(complete_bytes)")

    if mode == "direct" || mode == "both"
        ir, lowering_time, lowering_bytes = timed(() -> QC._lower_moment_ir(eqs))
        (rhs, direct_p), evaluator_time, evaluator_bytes = timed(
            () -> direct_evaluator(eqs, ir, ps),
        )
        direct_f, function_time, function_bytes =
            timed(() -> QC.SciMLBase.ODEFunction{true}(rhs))
        direct_problem, ode_time, ode_bytes = timed(
            () -> QC.SciMLBase.ODEProblem(direct_f, u0, TSPAN, direct_p),
        )
        _, first_solve_time, first_solve_bytes =
            timed(() -> solve(direct_problem, Tsit5(); saveat = TSPAN[2]))

        du = similar(u0)
        direct_f, direct_p = direct_problem.f, direct_problem.p
        post_solve_rhs = @benchmark $direct_f($du, $u0, $direct_p, 0.0) samples = 100 evals = 1
        update_values = Dict(first(keys(ps)) => 1.1)
        _, update_time, update_bytes =
            timed(() -> update_parameters!(direct_problem, update_values))

        println("direct order=$order equations=$nstates")
        println("  direct lowering time=$(lowering_time)s allocations=$(lowering_bytes)")
        println("  direct evaluator time=$(evaluator_time)s allocations=$(evaluator_bytes)")
        println("  direct ODEFunction time=$(function_time)s allocations=$(function_bytes)")
        println("  direct ODEProblem time=$(ode_time)s allocations=$(ode_bytes)")
        println("  direct first solve time=$(first_solve_time)s allocations=$(first_solve_bytes)")
        println(
            "  direct post-solve RHS median=$(ns_to_seconds(median(post_solve_rhs).time))s allocations=$(median(post_solve_rhs).allocs)",
        )
        println("  direct parameter update time=$(update_time)s allocations=$(update_bytes)")
        println("  direct evaluator bytes=$(Base.summarysize((rhs, direct_p)))")
    end

    if mode == "mtk" || mode == "both"
        sys, system_time, system_bytes =
            timed(() -> System(eqs; name = Symbol(:benchmark_, order)))
        compiled, compile_time, compile_bytes = timed(() -> mtkcompile(sys))
        mtk_values = parameter_map(compiled, merge(initial_values(eqs, u0), ps))
        mtk_problem, mtk_ode_time, mtk_ode_bytes =
            timed(() -> QC.SciMLBase.ODEProblem(compiled, mtk_values, TSPAN))
        _, mtk_solve_time, mtk_solve_bytes =
            timed(() -> solve(mtk_problem, Tsit5(); saveat = TSPAN[2]))
        mtk_du = similar(u0)
        mtk_f, mtk_p = mtk_problem.f, mtk_problem.p
        mtk_post_solve_rhs =
            @benchmark $mtk_f($mtk_du, $mtk_problem.u0, $mtk_p, 0.0) samples = 100 evals = 1

        println("MTK order=$order equations=$nstates")
        println("  MTK System time=$(system_time)s allocations=$(system_bytes)")
        println("  MTK compile time=$(compile_time)s allocations=$(compile_bytes)")
        println("  MTK ODEProblem time=$(mtk_ode_time)s allocations=$(mtk_ode_bytes)")
        println(
            "  MTK post-solve RHS median=$(ns_to_seconds(median(mtk_post_solve_rhs).time))s allocations=$(median(mtk_post_solve_rhs).allocs)",
        )
        println("  MTK first solve time=$(mtk_solve_time)s allocations=$(mtk_solve_bytes)")
        println("  MTK evaluator bytes=$(Base.summarysize(mtk_problem.f))")
    end
    return nothing
end

report(parse(Int, get(ENV, "QC_BENCH_ORDER", "3")))
