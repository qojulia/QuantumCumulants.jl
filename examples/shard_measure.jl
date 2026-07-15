# Measure ONE codegen path end-to-end in isolation. Set `PATH` (a Symbol) before include.
# Must be run in a FRESH session per path: compile time is shared across builds within a
# session (and can be perturbed by Revise), so back-to-back paths are not comparable.
#
# Metric = cold first-`solve` wall time and warm second-`solve`; compile ≈ first - warm.
# This is what the user actually experiences and what issue #294's table reports. It does not
# depend on guessing whether codegen happens in mtkcompile, ODEProblem, or the first call.
#
#   PATH  = :ode_default | :ode_eval | :legacy_serial | :legacy_sharded | :legacy_support | :legacy_rcm
#   ORDER = 1 | 2 | 3   (optional, default 3)

using QuantumCumulants
using ModelingToolkitBase
using ModelingToolkitBase: unknowns, equations, parameters, get_iv
using Symbolics
using Symbolics: ShardedForm
using OrdinaryDiffEqLowOrderRK   # exports ODEProblem, ODEFunction, RK4
using QuantumOptics
include(joinpath(@__DIR__, "shard_partition.jl"))

isdefined(Main, :ORDER) || (ORDER = 3)
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
neq = length(dvs)

ψ0 = tensor([spinup(SpinBasis(1 // 2)) for _ in 1:N]...)
u0 = initial_values(eqs, ψ0)
pdict = Dict(J => 1.0, hx => 1.0, γ => 0.2)
tspan = (0.0, 5.0)

# Time from the mtkcompiled sys through the first solution: problem construction (codegen may
# happen here) + first solve (codegen may happen here instead). Their sum captures all
# codegen+compile regardless of where in the path it lands.
t_prob = @elapsed prob = if PATH in (:ode_default, :ode_eval)
    dict = QuantumCumulants.parameter_map(sys, merge(Dict(unknowns(sys) .=> u0), pdict))
    ODEProblem(sys, dict, tspan; eval_expression = (PATH === :ode_eval), build_initializeprob = false)
else
    rhss = [eq.rhs for eq in equations(sys)]
    perm = PATH === :legacy_support ? support_ordering(eqs.graph) :
        PATH === :legacy_rcm ? rcm_ordering(eqs.graph) : collect(1:neq)
    parallel = PATH === :legacy_serial ? nothing : ShardedForm(10, 4)
    kw = parallel === nothing ? (; cse = true) : (; cse = true, parallel)
    # legacy path reorders equations; permute u0 to match so the solve is valid
    fexpr = build_function(rhss[perm], dvs[perm], ps, tiv; expression = Val{true}, kw...)[2]
    f! = eval(fexpr)
    pvec = Float64[pdict[x] for x in ps]
    ODEProblem(ODEFunction((du, u, p, t) -> Base.invokelatest(f!, du, u, p, t)), u0[perm], tspan, pvec)
end

t_cold = @elapsed solve(prob, RK4(); saveat = 0.5)
t_warm = @elapsed solve(prob, RK4(); saveat = 0.5)
total = t_prob + t_cold
println("RESULT order=$ORDER neq=$neq PATH=$PATH prob=$(round(t_prob, digits = 2))s first_solve=$(round(t_cold, digits = 2))s total_cold=$(round(total, digits = 2))s warm=$(round(t_warm, digits = 3))s compile≈$(round(total - t_warm, digits = 2))s")
