# The RGF gate (spec plan item 2): per-chunk RuntimeGeneratedFunctions instead of eval'd
# top-level functions. Checks (a) constructing per-chunk RGFs at top level does not hit the
# ShardedForm-style codegen deadlock, (b) after a parallel warm phase the first DIRECT call
# (no invokelatest) is warm, (c) the solve is correct, (d) honest timings including the
# System+mtkcompile step shard_manual.jl left untimed.
#
# Fresh session per run (same discipline as shard_measure.jl). Set before include:
#   WARM  = :call | :precompile   (parallel warm strategy; needs `julia -t <n>`)
#   CHUNK = equations per chunk   (default 10)
#   ORDER = 1 | 2 | 3             (default 3)

using QuantumCumulants
using ModelingToolkitBase
using ModelingToolkitBase: unknowns, equations, parameters, get_iv
using Symbolics
using OrdinaryDiffEqLowOrderRK
using QuantumOptics
using RuntimeGeneratedFunctions
RuntimeGeneratedFunctions.init(@__MODULE__)

isdefined(Main, :ORDER) || (ORDER = 3)
isdefined(Main, :WARM) || (WARM = :call)
isdefined(Main, :CHUNK) || (CHUNK = 10)
const N = 6

h = ⊗([PauliSpace(Symbol(:spin, i)) for i in 1:N]...)
σx(i) = Pauli(h, :σ, 1, i); σy(i) = Pauli(h, :σ, 2, i); σz(i) = Pauli(h, :σ, 3, i)
σm(i) = (σx(i) - 1im * σy(i)) / 2
@variables J hx γ
H = -J * sum(σz(i) * σz(i + 1) for i in 1:(N - 1)) - hx * sum(σx(i) for i in 1:N)
eqs = meanfield([σz(i) for i in 1:N], H, [σm(i) for i in 1:N]; rates = [γ for i in 1:N], order = ORDER)
complete!(eqs)
t_sys = @elapsed sys = mtkcompile(System(eqs; name = :sys))
dvs = unknowns(sys); ps = parameters(sys); tiv = get_iv(sys)
rhss = [eq.rhs for eq in equations(sys)]
neq = length(rhss)

ψ0 = tensor([spinup(SpinBasis(1 // 2)) for _ in 1:N]...)
u0 = initial_values(eqs, ψ0)
pvec = Float64[Dict(J => 1.0, hx => 1.0, γ => 0.2)[x] for x in ps]
tspan = (0.0, 5.0)

ranges = [r[1]:r[end] for r in Iterators.partition(1:neq, CHUNK)]
t_build = @elapsed fexprs = [
    Symbolics.build_function(rhss[r], dvs, ps, tiv; expression = Val{true}, cse = true)[2]
        for r in ranges
]
# the deadlock check: one RGF per chunk, constructed at top level (NOT inside codegen)
t_rgf = @elapsed fs = [@RuntimeGeneratedFunction(fe) for fe in fexprs]

Ts = (
    SubArray{ComplexF64, 1, Vector{ComplexF64}, Tuple{UnitRange{Int64}}, true},
    Vector{ComplexF64}, Vector{Float64}, Float64,
)
t_warmphase = @elapsed if WARM === :precompile
    Threads.@threads for f in fs
        precompile(f, Ts)
    end
else
    Threads.@threads for i in eachindex(fs)
        buf = zeros(ComplexF64, neq)
        fs[i](view(buf, ranges[i]), u0, pvec, 0.0)
    end
end

function driver!(du, u, p, t)
    for (f, r) in zip(fs, ranges)
        f(view(du, r), u, p, t)   # direct call; RGFs need no invokelatest
    end
    return nothing
end

# warmness of the first direct call after the parallel warm phase
du = zeros(ComplexF64, neq)
t_call1 = @elapsed driver!(du, u0, pvec, 0.0)
t_call2 = @elapsed driver!(du, u0, pvec, 0.0)

t_prob = @elapsed prob = ODEProblem(ODEFunction(driver!), u0, tspan, pvec)
t_cold = @elapsed sol = solve(prob, RK4(); saveat = 0.5)
t_warm = @elapsed solve(prob, RK4(); saveat = 0.5)
total = t_build + t_rgf + t_warmphase + t_prob + t_cold
finalmax = maximum(abs, sol.u[end])
println(
    "RESULT order=$ORDER neq=$neq WARM=$WARM chunk=$CHUNK nthreads=$(Threads.nthreads()) " *
        "sys=$(round(t_sys, digits = 2))s build=$(round(t_build, digits = 2))s " *
        "rgf=$(round(t_rgf, digits = 2))s warmphase=$(round(t_warmphase, digits = 2))s " *
        "call1=$(round(t_call1 * 1000, digits = 2))ms call2=$(round(t_call2 * 1000, digits = 2))ms " *
        "prob=$(round(t_prob, digits = 2))s first_solve=$(round(t_cold, digits = 2))s " *
        "total_cold=$(round(total, digits = 2))s warm=$(round(t_warm, digits = 3))s " *
        "finalmax=$(round(finalmax, digits = 8))",
)
