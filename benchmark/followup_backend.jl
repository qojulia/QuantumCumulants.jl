# Reproducible follow-up measurements for the structured direct backend.
#
# Every invocation is intended for a fresh Julia process. Package loading and
# precompilation happen before the measured model pipeline and are reported as
# excluded. Results are appended to QC_FOLLOWUP_RESULTS as each stage finishes,
# so a long MTK run can be polled after the Julia tool call returns.
#
# Examples:
#   QC_FOLLOWUP_MODE=compact QC_FOLLOWUP_ORDER=3 \
#     julia --project=benchmark benchmark/followup_backend.jl
#   QC_FOLLOWUP_MODE=mtk QC_FOLLOWUP_ORDER=4 \
#     julia --project=benchmark benchmark/followup_backend.jl
#   QC_FOLLOWUP_MODE=generated QC_FOLLOWUP_ORDER=3 QC_GEN_TARGET=256 \
#     julia --project=benchmark benchmark/followup_backend.jl

using BenchmarkTools
using ModelingToolkitBase: mtkcompile
using OrdinaryDiffEqTsit5: Tsit5, solve
using Printf
using QuantumCumulants
using SciMLBase: ODEProblem
using SparseArrays
using RuntimeGeneratedFunctions

const QC = QuantumCumulants
const N = 6
const SHORT_TSPAN = (0.0, 0.01)
const LONG_TEND = parse(Float64, get(ENV, "QC_FOLLOWUP_LONG_TEND", "1.0"))
const RESULTS = get(
    ENV,
    "QC_FOLLOWUP_RESULTS",
    joinpath(@__DIR__, "results", "followup-$(get(ENV, "QC_FOLLOWUP_MODE", "compact"))-order$(get(ENV, "QC_FOLLOWUP_ORDER", "3")).log"),
)

RuntimeGeneratedFunctions.init(@__MODULE__)

mkpath(dirname(RESULTS))
open(RESULTS, "w") do io
    println(io, "# package-load and precompilation time excluded")
end

function emit(line)
    println(line)
    open(RESULTS, "a") do io
        println(io, line)
    end
    flush(stdout)
    return line
end

function timed(f)
    GC.gc()
    result = @timed f()
    return result.value, result.time, result.bytes
end

function benchmark_stats(f; samples = 100)
    trial = @benchmark $f() samples = samples evals = 1
    med = median(trial)
    return (;
        median_seconds = med.time / 1.0e9, allocations = minimum(trial).allocs,
        allocated_bytes = minimum(trial).memory,
    )
end

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

function nonzero_state(seed, n)
    return ComplexF64[
        0.11 * sin(seed + 0.17i) + 0.07im * cos(0.31seed + 0.11i) for i in 1:n
    ]
end

function parameter_update(ps)
    key = first(keys(ps))
    return Dict(key => ps[key] + 0.37)
end

function stage_header(order, mode)
    emit("metadata julia=$(VERSION) threads=$(Threads.nthreads()) order=$order mode=$mode long_tend=$LONG_TEND")
    return emit("metadata package_load=excluded warm_runs_are_not_used_for_cold_stages=true")
end

function common_pipeline(order)
    open_eqs, ps, meanfield_time, meanfield_bytes = let
        result, t, b = timed(() -> ising_model(order))
        result[1], result[2], t, b
    end
    eqs, complete_time, complete_bytes = timed(() -> complete(open_eqs))
    emit("stage meanfield seconds=$(meanfield_time) bytes=$(meanfield_bytes)")
    emit("stage complete seconds=$(complete_time) bytes=$(complete_bytes) equations=$(length(eqs.states))")
    return eqs, ps
end

function direct_setup(eqs, ps)
    ir, lowering_time, lowering_bytes = timed(() -> QC._lower_moment_ir(eqs))
    kernel, kernel_time, kernel_bytes = timed(() -> QC.MomentKernel(ir; parallel = false))
    parameters, parameter_time, parameter_bytes = timed(
        () -> QC.KernelParameters(ir, kernel.pattern, QC.kernel_pdict(ir.params, parameter_map(eqs, ps))),
    )
    rhs = QC.KernelRHS(kernel)
    direct_f, function_time, function_bytes = timed(() -> QC._ode_function(rhs))
    emit("stage moment_ir seconds=$(lowering_time) bytes=$(lowering_bytes) monomials=$(length(ir.parent)) coo=$(length(ir.coo_i))")
    emit("stage compact_kernel seconds=$(kernel_time) bytes=$(kernel_bytes)")
    emit("stage parameter_payload seconds=$(parameter_time) bytes=$(parameter_bytes)")
    emit("stage ode_function seconds=$(function_time) bytes=$(function_bytes)")
    return (;
        ir, kernel, parameters, rhs, direct_f,
        lowering_time, kernel_time, parameter_time, function_time,
    )
end

function compact_problem(eqs, setup)
    u0 = nonzero_state(0, setup.ir.nstates)
    problem, problem_time, problem_bytes = timed(
        () -> ODEProblem(setup.direct_f, u0, SHORT_TSPAN, setup.parameters),
    )
    solution, solve_time, solve_bytes = timed(
        () -> solve(problem, Tsit5(); saveat = SHORT_TSPAN[2]),
    )
    emit("stage compact_problem seconds=$(problem_time) bytes=$(problem_bytes)")
    emit("stage compact_first_solution seconds=$(solve_time) bytes=$(solve_bytes) saved=$(length(solution.t))")
    emit("compact_total_after_complete seconds=$(setup.lowering_time + setup.kernel_time + setup.parameter_time + setup.function_time + problem_time + solve_time)")
    return problem, solution
end

function compact_profile(setup)
    u = nonzero_state(1, setup.ir.nstates)
    du = similar(u)
    v = QC._vbuf(setup.kernel.v)
    scratch = benchmark_stats(() -> QC._vbuf(setup.kernel.v))
    monomial = benchmark_stats(() -> QC.update_v!(v, setup.kernel.parent, setup.kernel.leaf, u))
    accumulation = benchmark_stats(
        () -> QC.spmv!(du, setup.kernel.pattern, setup.parameters.nzval, v, false),
    )
    rhs = benchmark_stats(() -> setup.direct_f(du, u, setup.parameters, 0.0))
    for (name, result) in (
            ("scratch_lookup", scratch), ("monomial_pass", monomial),
            ("sparse_accumulation", accumulation), ("compact_rhs", rhs),
        )
        emit("profile $name median_seconds=$(result.median_seconds) allocations=$(result.allocations) bytes=$(result.allocated_bytes)")
    end
    return nothing
end

function compact_long_integration(problem)
    u0 = nonzero_state(2, length(problem.u0))
    long_problem = ODEProblem(problem.f, u0, (0.0, LONG_TEND), problem.p)
    solution, solve_time, solve_bytes = timed(
        () -> solve(long_problem, Tsit5(); saveat = range(0.0, LONG_TEND; length = 101)),
    )
    emit("long_integration seconds=$(solve_time) bytes=$(solve_bytes) saved=$(length(solution.t))")
    return solution
end

const GeneratedStage = Symbol

struct GeneratedKernel{M, R, S, K}
    monomial_units::Vector{M}
    row_units::Vector{R}
    monomial_runner::Function
    row_runner::Function
    scratch::S
    compact::K
    stage::GeneratedStage
end

Base.@noinline function call_generated_unit(f, args...)
    f(args...)
    return nothing
end

function (k::GeneratedKernel)(du, u, p, t)
    v = QC._vbuf(k.scratch)
    if k.stage === :monomial || k.stage === :both
        call_generated_unit(k.monomial_runner, k.monomial_units, v, u)
    else
        QC.update_v!(v, k.compact.parent, k.compact.leaf, u)
    end
    if k.stage === :rows || k.stage === :both
        call_generated_unit(k.row_runner, k.row_units, du, p, v)
    else
        QC.spmv!(du, k.compact.pattern, p.nzval, v, false)
    end
    return nothing
end

function operation_ranges(costs::Vector{Int}, target::Int)
    target > 0 || throw(ArgumentError("QC_GEN_TARGET must be positive"))
    ranges = UnitRange{Int}[]
    lo = firstindex(costs)
    while lo <= lastindex(costs)
        hi = lo
        total = costs[lo]
        while hi < lastindex(costs) && total + costs[hi + 1] <= target
            hi += 1
            total += costs[hi]
        end
        push!(ranges, lo:hi)
        lo = hi + 1
    end
    return ranges
end

function define_monomial_unit!(range, parent, leaf)
    statements = Expr[]
    for m in range
        j = leaf[m]
        factor = j > 0 ? :(u[$j]) : :(Base.conj(u[$(-j)]))
        push!(statements, :(v[$m] = v[$(parent[m])] * $factor))
    end
    expr = :(
        (v, u) -> begin
            @inbounds begin
                $(statements...)
            end
            return nothing
        end
    )
    return RuntimeGeneratedFunction(@__MODULE__, @__MODULE__, expr)
end

function define_row_unit!(range, pattern)
    statements = Expr[]
    for i in range
        terms = Expr[]
        for k in pattern.colptr[i]:(pattern.colptr[i + 1] - 1)
            push!(terms, :(p.nzval[$k] * v[$(pattern.rowval[k])]))
        end
        value = isempty(terms) ? :(zero(ComplexF64)) : foldl((a, b) -> :($a + $b), terms)
        push!(statements, :(du[$i] = $value))
    end
    expr = :(
        (du, p, v) -> begin
            @inbounds begin
                $(statements...)
            end
            return nothing
        end
    )
    return RuntimeGeneratedFunction(@__MODULE__, @__MODULE__, expr)
end

function define_runner!(units, args)
    # Keep the units in a Vector{Function} and call through a noinline wrapper.
    # This prevents the dispatcher from folding all units into one giant method.
    expr = args == (:v, :u) ?
        :(
            (units, v, u) -> begin
                for unit in units
                    call_generated_unit(unit, v, u)
            end
                return nothing
            end
        ) :
        :(
            (units, du, p, v) -> begin
                for unit in units
                    call_generated_unit(unit, du, p, v)
            end
                return nothing
            end
        )
    return RuntimeGeneratedFunction(@__MODULE__, @__MODULE__, expr)
end

function generated_setup(eqs, ps, target, stage)
    ir, lowering_time, lowering_bytes = timed(() -> QC._lower_moment_ir(eqs))
    pattern = QC.assemble_pattern(ir)
    parameters, parameter_time, parameter_bytes = timed(
        () -> QC.KernelParameters(ir, pattern, QC.kernel_pdict(ir.params, parameter_map(eqs, ps))),
    )
    compact = QC.MomentKernel(ir; parallel = false)
    monomial_costs = Int[compact.fac_ptr[m + 1] - compact.fac_ptr[m] for m in 1:length(ir.parent)]
    row_costs = Int[pattern.colptr[i + 1] - pattern.colptr[i] for i in 1:ir.nstates]
    monomial_ranges = operation_ranges(monomial_costs[2:end], target)
    monomial_ranges = [(first(r) + 1):(last(r) + 1) for r in monomial_ranges]
    row_ranges = operation_ranges(row_costs, target)
    monomial_units = stage in (:monomial, :both) ?
        Function[define_monomial_unit!(r, ir.parent, ir.leaf) for r in monomial_ranges] :
        Function[]
    row_units = stage in (:rows, :both) ?
        Function[define_row_unit!(r, pattern) for r in row_ranges] :
        Function[]
    monomial_runner = define_runner!(monomial_units, (:v, :u))
    row_runner = define_runner!(row_units, (:du, :p, :v))
    generated = GeneratedKernel(
        monomial_units, row_units, monomial_runner, row_runner,
        QC._make_vbufs(length(ir.parent)), compact, stage
    )
    emit("stage moment_ir seconds=$(lowering_time) bytes=$(lowering_bytes) monomials=$(length(ir.parent)) coo=$(length(ir.coo_i))")
    emit("stage parameter_payload seconds=$(parameter_time) bytes=$(parameter_bytes)")
    emit("stage generated_units target=$(target) monomial_units=$(length(monomial_units)) row_units=$(length(row_units)) monomial_ops=$(sum(monomial_costs)) row_ops=$(sum(row_costs))")
    return (; ir, pattern, parameters, generated, lowering_time, parameter_time)
end

function validate_generated(compact, generated_setup)
    cdu = similar(nonzero_state(0, compact.ir.nstates))
    gdu = similar(cdu)
    states = [nonzero_state(seed, compact.ir.nstates) for seed in 1:3]
    errors = Float64[]
    for u in states
        compact.direct_f(cdu, u, compact.parameters, 0.0)
        generated_setup.generated(gdu, u, generated_setup.parameters, 0.0)
        push!(errors, maximum(abs.(cdu .- gdu)))
    end
    update = Dict(first(keys(compact.parameters.values)) => 0.57)
    compact_updated = copy(compact.parameters)
    generated_updated = copy(generated_setup.parameters)
    QC.update_parameters!(compact_updated, update)
    QC.update_parameters!(generated_updated, update)
    updated_error = let u = states[2]
        compact.direct_f(cdu, u, compact_updated, 0.0)
        generated_setup.generated(gdu, u, generated_updated, 0.0)
        maximum(abs.(cdu .- gdu))
    end
    u0 = states[3]
    compact_problem = ODEProblem(compact.direct_f, u0, SHORT_TSPAN, compact.parameters)
    generated_problem = ODEProblem(generated_setup.generated, u0, SHORT_TSPAN, generated_setup.parameters)
    compact_solution = solve(compact_problem, Tsit5(); saveat = range(0.0, SHORT_TSPAN[2]; length = 11))
    generated_solution = solve(generated_problem, Tsit5(); saveat = range(0.0, SHORT_TSPAN[2]; length = 11))
    trajectory_error = maximum(maximum(abs.(a .- b)) for (a, b) in zip(compact_solution.u, generated_solution.u))
    emit("validation rhs_states=$(join(errors, ','))")
    emit("validation rhs_updated_parameter=$(updated_error)")
    emit("validation trajectory_max_error=$(trajectory_error) samples=$(length(compact_solution.t))")
    return (; errors, updated_error, trajectory_error)
end

function compact_run(order)
    stage_header(order, "compact")
    eqs, ps = common_pipeline(order)
    setup = direct_setup(eqs, ps)
    problem, _ = compact_problem(eqs, setup)
    compact_profile(setup)
    return compact_long_integration(problem)
end

function generated_run(order, target, stage)
    stage_header(order, "generated")
    eqs, ps = common_pipeline(order)
    generated, generated_setup_time, generated_setup_bytes = timed(
        () -> generated_setup(eqs, ps, target, stage),
    )
    # Build compact only after generated cold construction. It is needed for
    # validation, but must not warm the generated lowering path.
    compact = direct_setup(eqs, ps)
    u0 = nonzero_state(0, generated.ir.nstates)
    generated_f, function_time, function_bytes = timed(() -> QC._ode_function(generated.generated))
    problem, problem_time, problem_bytes = timed(
        () -> ODEProblem(generated_f, u0, SHORT_TSPAN, generated.parameters),
    )
    _, compile_time, compile_bytes = timed(() -> generated_f(similar(u0), u0, problem.p, 0.0))
    _, solve_time, solve_bytes = timed(() -> solve(problem, Tsit5(); saveat = SHORT_TSPAN[2]))
    emit("stage generated_ode_function seconds=$(function_time) bytes=$(function_bytes)")
    emit("stage generated_problem seconds=$(problem_time) bytes=$(problem_bytes)")
    emit("stage generated_first_call_compile seconds=$(compile_time) bytes=$(compile_bytes)")
    emit("stage generated_first_solution seconds=$(solve_time) bytes=$(solve_bytes)")
    emit("stage generated_setup_total seconds=$(generated_setup_time) bytes=$(generated_setup_bytes)")
    warm = benchmark_stats(() -> generated_f(similar(u0), u0, problem.p, 0.0))
    emit("profile generated_rhs median_seconds=$(warm.median_seconds) allocations=$(warm.allocations) bytes=$(warm.allocated_bytes)")
    validation = validate_generated(compact, generated)
    emit("generated_total_after_complete seconds=$(generated_setup_time + function_time + problem_time + compile_time + solve_time)")
    return validation
end

function mtk_run(order)
    stage_header(order, "mtk")
    eqs, ps = common_pipeline(order)
    u0 = nonzero_state(0, length(eqs.states))
    sys, system_time, system_bytes = timed(() -> System(eqs; name = Symbol(:followup_, order)))
    compiled, compile_time, compile_bytes = timed(() -> mtkcompile(sys))
    values, parameter_time, parameter_bytes = timed(
        () -> parameter_map(compiled, merge(initial_values(eqs, u0), ps)),
    )
    problem, problem_time, problem_bytes = timed(() -> ODEProblem(compiled, values, SHORT_TSPAN))
    _, solve_time, solve_bytes = timed(() -> solve(problem, Tsit5(); saveat = SHORT_TSPAN[2]))
    emit("stage mtk_system seconds=$(system_time) bytes=$(system_bytes)")
    emit("stage mtk_compile seconds=$(compile_time) bytes=$(compile_bytes)")
    emit("stage mtk_parameter_payload seconds=$(parameter_time) bytes=$(parameter_bytes)")
    emit("stage mtk_problem seconds=$(problem_time) bytes=$(problem_bytes)")
    emit("stage mtk_first_solution seconds=$(solve_time) bytes=$(solve_bytes)")
    warm = benchmark_stats(() -> problem.f(similar(u0), problem.u0, problem.p, 0.0))
    emit("profile mtk_rhs median_seconds=$(warm.median_seconds) allocations=$(warm.allocations) bytes=$(warm.allocated_bytes)")
    return emit("mtk_total_after_complete seconds=$(system_time + compile_time + parameter_time + problem_time + solve_time)")
end

mode = get(ENV, "QC_FOLLOWUP_MODE", "compact")
order = parse(Int, get(ENV, "QC_FOLLOWUP_ORDER", "3"))
target = parse(Int, get(ENV, "QC_GEN_TARGET", "256"))
stage = Symbol(get(ENV, "QC_GEN_STAGE", "both"))
mode in ("compact", "generated", "mtk") || throw(ArgumentError("QC_FOLLOWUP_MODE must be compact, generated, or mtk"))
stage in (:monomial, :rows, :both) || throw(ArgumentError("QC_GEN_STAGE must be monomial, rows, or both"))
mode == "compact" && compact_run(order)
mode == "generated" && generated_run(order, target, stage)
mode == "mtk" && mtk_run(order)
