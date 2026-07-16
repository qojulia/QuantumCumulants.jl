# Manual sharding: K separate TOP-LEVEL chunk functions (one `build_function` per chunk of
# equations, each `eval`'d individually) instead of ShardedForm's nested RGF tree. This
# unlocks the one lever ShardedForm cannot reach: compiling the chunks on multiple threads
# via `precompile`, since each chunk is an independent top-level function.
#
# Fresh session per run (same discipline as shard_measure.jl). Set before include:
#   MODE      = :serial | :threaded   (threaded needs `julia -t <n>`)
#   CHUNK     = equations per chunk   (default 10, matching ShardedForm leaf size)
#   ORDER     = 1 | 2 | 3             (default 3)

using QuantumCumulants
using ModelingToolkitBase
using ModelingToolkitBase: unknowns, equations, parameters, get_iv
using Symbolics
using OrdinaryDiffEqLowOrderRK
using QuantumOptics

isdefined(Main, :ORDER) || (ORDER = 3)
isdefined(Main, :MODE) || (MODE = :serial)
isdefined(Main, :CHUNK) || (CHUNK = 10)
const N = 6

h = ⊗([PauliSpace(Symbol(:spin, i)) for i in 1:N]...)
σx(i) = Pauli(h, :σ, 1, i); σy(i) = Pauli(h, :σ, 2, i); σz(i) = Pauli(h, :σ, 3, i)
σm(i) = (σx(i) - 1im * σy(i)) / 2
@variables J hx γ
H = -J * sum(σz(i) * σz(i + 1) for i in 1:(N - 1)) - hx * sum(σx(i) for i in 1:N)
eqs = meanfield([σz(i) for i in 1:N], H, [σm(i) for i in 1:N]; rates = [γ for i in 1:N], order = ORDER)
complete!(eqs)
sys = mtkcompile(System(eqs; name = :sys))
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
t_eval = @elapsed fs = [eval(fe) for fe in fexprs]

t_pc = 0.0
if MODE === :threaded
    Ts = (
        SubArray{ComplexF64, 1, Vector{ComplexF64}, Tuple{UnitRange{Int64}}, true},
        Vector{ComplexF64}, Vector{Float64}, Float64,
    )
    t_pc = @elapsed Threads.@threads for f in fs
        precompile(f, Ts)
    end
end

function driver!(du, u, p, t)
    for (f, r) in zip(fs, ranges)
        Base.invokelatest(f, view(du, r), u, p, t)
    end
    return nothing
end

t_prob = @elapsed prob = ODEProblem(ODEFunction(driver!), u0, tspan, pvec)
t_cold = @elapsed solve(prob, RK4(); saveat = 0.5)
t_warm = @elapsed solve(prob, RK4(); saveat = 0.5)
total = t_build + t_eval + t_pc + t_prob + t_cold
println(
    "RESULT order=$ORDER neq=$neq MODE=$MODE chunk=$CHUNK nthreads=$(Threads.nthreads()) " *
        "build=$(round(t_build, digits = 2))s eval=$(round(t_eval, digits = 2))s " *
        "precompile=$(round(t_pc, digits = 2))s prob=$(round(t_prob, digits = 2))s " *
        "first_solve=$(round(t_cold, digits = 2))s total_cold=$(round(total, digits = 2))s " *
        "warm=$(round(t_warm, digits = 3))s",
)
