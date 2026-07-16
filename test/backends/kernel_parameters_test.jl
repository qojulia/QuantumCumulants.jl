using QuantumCumulants
using Symbolics: Symbolics, @variables
using SciMLBase: SciMLBase, ODEProblem
using OrdinaryDiffEqLowOrderRK: RK4, solve
using Test

# Transverse-field Ising chain of 3 Pauli spins at order 2 (same model as
# kernel_backend_test.jl; test files run standalone under ParallelTestRunner).
Np = 3
h = ⊗([PauliSpace(Symbol(:spin, i)) for i in 1:Np]...)
σx(i) = Pauli(h, :σ, 1, i)
σy(i) = Pauli(h, :σ, 2, i)
σz(i) = Pauli(h, :σ, 3, i)
σm(i) = (σx(i) - 1im * σy(i)) / 2
@variables J hx γ
H = -J * sum(σz(i) * σz(i + 1) for i in 1:(Np - 1)) - hx * sum(σx(i) for i in 1:Np)
eqs = meanfield(
    [σz(i) for i in 1:Np], H, [σm(i) for i in 1:Np];
    rates = [γ for i in 1:Np], order = 2,
)
complete!(eqs)
ps = Dict(J => 1.0, hx => 1.0, γ => 0.2)
pd2 = Dict(J => 1.7, hx => 0.4, γ => 0.31)
nst = length(eqs.states)
u0 = ComplexF64[k == 1 ? 1.0 : 0.0 for k in 1:nst]
u = ComplexF64[0.1cos(3.7i) + 0.05im * sin(1.3i) for i in 1:nst]

@testset "update_parameters! matches a fresh construction bit-exactly" begin
    prob = ODEProblem(eqs, u0, (0.0, 5.0), ps; backend = KernelBackend())
    update_parameters!(prob, pd2)
    prob_fresh = ODEProblem(eqs, u0, (0.0, 5.0), pd2; backend = KernelBackend())
    du_a, du_b = similar(u0), similar(u0)
    prob.f(du_a, u, prob.p, 0.0)
    prob_fresh.f(du_b, u, prob_fresh.p, 0.0)
    @test du_a == du_b
end

@testset "copy isolation for ensembles" begin
    prob = ODEProblem(eqs, u0, (0.0, 5.0), ps; backend = KernelBackend())
    update_parameters!(prob, pd2)
    du_a, du_b = similar(u0), similar(u0)
    prob.f(du_a, u, prob.p, 0.0)

    f2 = copy(prob.f)
    prob2 = SciMLBase.remake(prob; f = f2, p = f2.f.kp)
    update_parameters!(prob2, Dict(J => 0.3))
    # the original is untouched by the copy's update
    prob.f(du_b, u, prob.p, 0.0)
    @test du_a == du_b
    # and the two problems genuinely evolve differently
    sol1 = solve(prob, RK4(); saveat = 0.5)
    sol2 = solve(prob2, RK4(); saveat = 0.5)
    @test maximum(maximum(abs, a .- b) for (a, b) in zip(sol1.u, sol2.u)) > 1.0e-3
end
