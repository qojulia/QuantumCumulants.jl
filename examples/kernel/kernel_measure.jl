# Cold end-to-end measurement of the M·v kernel path (issue #294). Fresh session per run,
# same discipline and same model/metrics as shard_measure.jl. The kernel path starts from
# the completed equations and never touches System/mtkcompile; the RESULT line reports both
# the from-completed-eqs total (comparable to shard_measure's prob+first_solve) and the
# grand total including model construction.
#
#   ORDER = 1 | 2 | 3   (optional, default 3)

using QuantumCumulants
using Symbolics
using SymbolicUtils
using SparseArrays
using OrdinaryDiffEqLowOrderRK
using QuantumOptics

isdefined(Main, :ORDER) || (ORDER = 3)
isdefined(Main, :NSPIN) || (NSPIN = 6)
const N = NSPIN
include(joinpath(@__DIR__, "kernel_lower.jl"))
include(joinpath(@__DIR__, "kernel_eval.jl"))

h = ⊗([PauliSpace(Symbol(:spin, i)) for i in 1:N]...)
σx(i) = Pauli(h, :σ, 1, i); σy(i) = Pauli(h, :σ, 2, i); σz(i) = Pauli(h, :σ, 3, i)
σm(i) = (σx(i) - 1im * σy(i)) / 2
@variables J hx γ
t_model = @elapsed begin
    eqs = meanfield(
        [σz(i) for i in 1:N], -J * sum(σz(i) * σz(i + 1) for i in 1:(N - 1)) -
            hx * sum(σx(i) for i in 1:N), [σm(i) for i in 1:N];
        rates = [γ for i in 1:N], order = ORDER,
    )
    complete!(eqs)
end
neq = length(eqs.graph.nodes)

t_lower = @elapsed ir = lower(eqs)
t_coeff = @elapsed cvals = coefficient_values(ir, Dict(J => 1.0, hx => 1.0, γ => 0.2))
t_mat = @elapsed k = MomentKernel(ir, cvals)

ψ0 = tensor([spinup(SpinBasis(1 // 2)) for _ in 1:N]...)
u0 = initial_values(eqs, ψ0)
tspan = (0.0, 5.0)
t_prob = @elapsed prob = ODEProblem(ODEFunction{true}(k), u0, tspan, cvals)
t_cold = @elapsed solve(prob, RK4(); saveat = 0.5)
t_warm = @elapsed solve(prob, RK4(); saveat = 0.5)

from_eqs = t_lower + t_coeff + t_mat + t_prob + t_cold
println(
    "RESULT order=$ORDER neq=$neq PATH=kernel nmono=$(length(ir.parent)) " *
        "nnz=$(nnz(k.M)) ncoeff=$(length(ir.coeffs)) model=$(round(t_model, digits = 2))s " *
        "lower=$(round(t_lower, digits = 2))s coeff=$(round(t_coeff, digits = 3))s " *
        "mat=$(round(t_mat, digits = 3))s prob=$(round(t_prob, digits = 2))s " *
        "first_solve=$(round(t_cold, digits = 2))s from_eqs=$(round(from_eqs, digits = 2))s " *
        "grand_total=$(round(t_model + from_eqs, digits = 2))s warm=$(round(t_warm, digits = 3))s",
)
