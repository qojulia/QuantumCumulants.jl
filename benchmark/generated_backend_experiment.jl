# PROTOTYPE: generated/native evaluator experiment from the existing MomentIR.
#
# This file is deliberately outside the production backend. It answers one question:
# can bounded native units improve the warm RHS or total time-to-solution enough to
# justify their construction and compilation cost? Delete it if the answer is no.
#
# Run each mode and order in a fresh Julia process:
#   QC_GEN_MODE=direct QC_GEN_ORDER=3 julia --project=benchmark benchmark/generated_backend_experiment.jl
#   QC_GEN_MODE=generated QC_GEN_ORDER=3 QC_GEN_UNIT=128 julia --project=benchmark benchmark/generated_backend_experiment.jl
#   QC_GEN_MODE=mtk QC_GEN_ORDER=3 julia --project=benchmark benchmark/generated_backend_experiment.jl
#
# The generated path preserves the MomentIR monomial IDs and parameter table. Only the
# monomial products and sparse row sums are emitted as bounded native functions.

using BenchmarkTools
using ModelingToolkitBase: mtkcompile
using OrdinaryDiffEqTsit5: Tsit5, solve
using QuantumCumulants
using RuntimeGeneratedFunctions
using SciMLBase: ODEProblem

const QC = QuantumCumulants
const N = 6
const TSPAN = (0.0, 0.01)

RuntimeGeneratedFunctions.init(@__MODULE__)

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
        rates = [γ for _ in 1:N],
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

struct GeneratedKernel{M, R, MS, RS, S}
    monomials::M
    rows::R
    monomial_units::MS
    row_units::RS
    scratch::S
end

function (k::GeneratedKernel)(du, u, p, t)
    v = QC._vbuf(k.scratch)
    k.monomials(k.monomial_units, v, u)
    k.rows(k.row_units, du, p, v)
    return nothing
end

function _define_monomial_unit!(lo, hi, parent, leaf)
    statements = Expr[]
    for m in lo:hi
        j = leaf[m]
        factor = j > 0 ? :(u[$j]) : :(Base.conj(u[$(-j)]))
        push!(statements, :(v[$m] = v[$(parent[m])] * $factor))
    end
    expr = :((v, u) -> begin
        @inbounds begin
            $(statements...)
        end
        return nothing
    end)
    return RuntimeGeneratedFunction(@__MODULE__, @__MODULE__, expr)
end

function _define_row_unit!(lo, hi, pattern)
    statements = Expr[]
    for i in lo:hi
        terms = Expr[]
        for k in pattern.colptr[i]:(pattern.colptr[i + 1] - 1)
            push!(terms, :(p.nzval[$k] * v[$(pattern.rowval[k])]))
        end
        value = isempty(terms) ? :(zero(ComplexF64)) : foldl((a, b) -> :($a + $b), terms)
        push!(statements, :(du[$i] = $value))
    end
    expr = :((du, p, v) -> begin
        @inbounds begin
            $(statements...)
        end
        return nothing
    end)
    return RuntimeGeneratedFunction(@__MODULE__, @__MODULE__, expr)
end

function _define_dispatch!(units, args)
    calls = [Expr(:call, Expr(:ref, :units, i), args...) for i in eachindex(units)]
    expr = :((units, $(args...)) -> begin
        $(calls...)
        return nothing
    end)
    return RuntimeGeneratedFunction(@__MODULE__, @__MODULE__, expr)
end

function generated_kernel(ir, ps, unit_size)
    unit_size > 0 || throw(ArgumentError("QC_GEN_UNIT must be positive"))
    pattern = QC.assemble_pattern(ir)
    values = QC.kernel_pdict(ir.params, ps)
    parameters = QC.KernelParameters(ir, pattern, values)
    monomial_units = Function[]
    for lo in 2:unit_size:length(ir.parent)
        hi = min(lo + unit_size - 1, length(ir.parent))
        push!(monomial_units, _define_monomial_unit!(lo, hi, ir.parent, ir.leaf))
    end
    row_units = Function[]
    for lo in 1:unit_size:ir.nstates
        hi = min(lo + unit_size - 1, ir.nstates)
        push!(row_units, _define_row_unit!(lo, hi, pattern))
    end
    monomial_dispatch = _define_dispatch!(monomial_units, (:v, :u))
    row_dispatch = _define_dispatch!(row_units, (:du, :p, :v))
    scratch = QC._make_vbufs(length(ir.parent))
    return GeneratedKernel(monomial_dispatch, row_dispatch, Tuple(monomial_units), Tuple(row_units), scratch), parameters, pattern,
        length(monomial_units), length(row_units)
end

function direct_kernel(ir, ps)
    values = QC.kernel_pdict(ir.params, ps)
    kernel = QC.MomentKernel(ir; parallel = false)
    return QC.KernelRHS(kernel), QC.KernelParameters(ir, kernel.pattern, values)
end

function mtk_problem(eqs, ps, u0)
    sys, system_time, system_bytes = timed(() -> System(eqs; name = :generated_experiment))
    compiled, compile_time, compile_bytes = timed(() -> mtkcompile(sys))
    values = parameter_map(compiled, merge(initial_values(eqs, u0), ps))
    problem, problem_time, problem_bytes = timed(() -> ODEProblem(compiled, values, TSPAN))
    return problem, system_time, system_bytes, compile_time, compile_bytes, problem_time, problem_bytes
end

function report(order, mode, unit_size)
    open, ps = ising_model(order)
    _, meanfield_time, meanfield_bytes = timed(() -> ising_model(order)[1])
    eqs, complete_time, complete_bytes = timed(() -> complete(open))
    ir, lowering_time, lowering_bytes = timed(() -> QC._lower_moment_ir(eqs))
    u0 = zeros(ComplexF64, ir.nstates)

    println("order=$order equations=$(ir.nstates) monomials=$(length(ir.parent)) nnz=$(length(ir.coo_i)) mode=$mode unit=$unit_size")
    println("  common meanfield time=$(meanfield_time)s allocations=$(meanfield_bytes)")
    println("  common complete time=$(complete_time)s allocations=$(complete_bytes)")
    println("  common MomentIR time=$(lowering_time)s allocations=$(lowering_bytes)")

    if mode == "direct" || mode == "all"
        rhs, direct_p = direct_kernel(ir, ps)
        direct_f, function_time, function_bytes = timed(() -> QC._ode_function(rhs))
        direct_problem, problem_time, problem_bytes = timed(() -> ODEProblem(direct_f, u0, TSPAN, direct_p))
        _, solve_time, solve_bytes = timed(() -> solve(direct_problem, Tsit5(); saveat = TSPAN[2]))
        du = similar(u0)
        warm = @benchmark $direct_f($du, $u0, $direct_p, 0.0) samples = 100 evals = 1
        println("direct function=$(function_time)s problem=$(problem_time)s first_solve=$(solve_time)s total=$(function_time + problem_time + solve_time)s")
        println("  allocations function=$(function_bytes) problem=$(problem_bytes) solve=$(solve_bytes)")
        println("  warm_rhs median=$(ns_to_seconds(median(warm).time))s allocs=$(median(warm).allocs) bytes=$(Base.summarysize(rhs))")
    end

    if mode == "generated" || mode == "all"
        (generated, generated_p, pattern, nmonomial_units, nrow_units), build_time, build_bytes =
            timed(() -> generated_kernel(ir, ps, unit_size))
        generated_f, function_time, function_bytes = timed(() -> QC._ode_function(generated))
        generated_problem, problem_time, problem_bytes =
            timed(() -> ODEProblem(generated_f, u0, TSPAN, generated_p))
        du = similar(u0)
        compile_time, compile_bytes = let
            _, t, b = timed(() -> generated_f(du, u0, generated_problem.p, 0.0))
            t, b
        end
        _, solve_time, solve_bytes = timed(() -> solve(generated_problem, Tsit5(); saveat = TSPAN[2]))
        warm = @benchmark $generated_f($du, $u0, $generated_problem.p, 0.0) samples = 100 evals = 1
        println("generated units monomial=$nmonomial_units rows=$nrow_units")
        println("generated build=$(build_time)s compile=$(compile_time)s function=$(function_time)s problem=$(problem_time)s first_solve=$(solve_time)s total=$(build_time + compile_time + function_time + problem_time + solve_time)s")
        println("  allocations build=$(build_bytes) compile=$(compile_bytes) function=$(function_bytes) problem=$(problem_bytes) solve=$(solve_bytes)")
        println("  warm_rhs median=$(ns_to_seconds(median(warm).time))s allocs=$(median(warm).allocs) bytes=$(Base.summarysize(generated))")
        # The generated path and compact path share the same numeric parameter payload shape.
        reference_rhs, reference_p = direct_kernel(ir, ps)
        reference_du = similar(u0)
        generated_f(du, u0, generated_problem.p, 0.0)
        reference_rhs(reference_du, u0, reference_p, 0.0)
        println("  single_rhs_max_error=$(maximum(abs.(du .- reference_du)))")
        println("  generated_pattern_nnz=$(length(pattern.nzval))")
    end

    if mode == "mtk" || mode == "all"
        problem, system_time, system_bytes, compile_time, compile_bytes, problem_time, problem_bytes =
            mtk_problem(eqs, ps, u0)
        _, solve_time, solve_bytes = timed(() -> solve(problem, Tsit5(); saveat = TSPAN[2]))
        du = similar(u0)
        warm = @benchmark $problem.f($du, $problem.u0, $problem.p, 0.0) samples = 100 evals = 1
        println("MTK system=$(system_time)s compile=$(compile_time)s problem=$(problem_time)s first_solve=$(solve_time)s total=$(system_time + compile_time + problem_time + solve_time)s")
        println("  allocations system=$(system_bytes) compile=$(compile_bytes) problem=$(problem_bytes) solve=$(solve_bytes)")
        println("  warm_rhs median=$(ns_to_seconds(median(warm).time))s allocs=$(median(warm).allocs) bytes=$(Base.summarysize(problem.f))")
    end
end

mode = get(ENV, "QC_GEN_MODE", "all")
mode in ("direct", "generated", "mtk", "all") ||
    throw(ArgumentError("QC_GEN_MODE must be `direct`, `generated`, `mtk`, or `all`; got $mode"))
order = parse(Int, get(ENV, "QC_GEN_ORDER", "3"))
unit_size = parse(Int, get(ENV, "QC_GEN_UNIT", "128"))
report(order, mode, unit_size)
