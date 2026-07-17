using QuantumCumulants
using Symbolics: Symbolics, @variables
using SciMLBase: SciMLBase, ODEProblem
using OrdinaryDiffEqLowOrderRK: RK4, solve
using OrdinaryDiffEqRosenbrock: Rodas5P
using QuantumOpticsBase: SpinBasis, spinup, tensor
using Test

# ShardedBackend Jacobian for stiff/implicit solves (issue #294). Both modes compute the
# real-directional linearization the integrator's own finite differencing would, so they
# agree with each other and with an explicit-solver reference, and (unlike the kernel's
# holomorphic M·v Jacobian) work for conjugate-folded closures too. A finite-difference
# time-gradient is attached, so plain `Rodas5P()` runs without an `autodiff` keyword.
#
# Note: solving a complex-state problem with a sparse `jac_prototype` triggers a benign
# "incompatible with sparse automatic differentiation" warning from the integrator; the
# supplied Jacobian is still used and the linear-solve stays sparse.

# Transverse-field Ising chain of 3 Pauli spins at order 2 (36 equations, unfolded).
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
nst = length(eqs.states)
u0 = zeros(ComplexF64, nst)
u = ComplexF64[0.1cos(3.7i) + 0.05im * sin(1.3i) for i in 1:nst]
ψ0 = tensor([spinup(SpinBasis(1 // 2)) for _ in 1:Np]...)

# reference trajectory from an explicit solver (needs no Jacobian)
ref = solve(
    ODEProblem(eqs, ψ0, (0.0, 5.0), ps; backend = ShardedBackend()), RK4();
    abstol = 1.0e-10, reltol = 1.0e-10, saveat = 0.5,
)

@testset "stiff solve with $mode Jacobian matches explicit reference" for mode in (:fd, true, :analytic)
    prob = ODEProblem(eqs, ψ0, (0.0, 5.0), ps; backend = ShardedBackend(), jac = mode)
    @test prob.f.jac !== nothing
    @test prob.f.tgrad !== nothing                        # so plain Rodas5P() needs no autodiff kw
    @test length(prob.f.jac_prototype.nzval) < nst^2      # stored sparsely, not dense
    sol = solve(prob, Rodas5P(); abstol = 1.0e-9, reltol = 1.0e-9, saveat = 0.5)
    @test sol.retcode == SciMLBase.ReturnCode.Success
    @test maximum(abs, sol.u[end] .- ref.u[end]) < 1.0e-5
end

@testset "analytic and colored-FD Jacobians agree" begin
    pfd = ODEProblem(eqs, u0, (0.0, 1.0), ps; backend = ShardedBackend(), jac = :fd)
    pan = ODEProblem(eqs, u0, (0.0, 1.0), ps; backend = ShardedBackend(), jac = :analytic)
    Jf = copy(pfd.f.jac_prototype)
    Ja = copy(pan.f.jac_prototype)
    pfd.f.jac(Jf, u, pfd.p, 0.0)
    pan.f.jac(Ja, u, pan.p, 0.0)
    @test maximum(abs, Array(Ja) .- Array(Jf)) < 1.0e-6   # both are the real-directional J
end

@testset "conjugate-folded closure is supported ($mode)" for mode in (:fd, :analytic)
    hc = FockSpace(:cavity)
    a = Destroy(hc, :a)
    @variables Δ Ω κ U
    Hk = Δ * a' * a + U * a' * a' * a * a + Ω * (a + a')
    eqs_k = meanfield([a], Hk, [a]; rates = [κ], order = 2)
    complete!(eqs_k; get_adjoints = false)          # folded: drift references conj(states)
    ps_k = Dict(Δ => -1.0, Ω => 1.3, κ => 1.0, U => 0.1)
    u0k = zeros(ComplexF64, length(eqs_k.states))
    # kernel rejects this (HolomorphicJacobianError); sharded handles it
    prob = ODEProblem(eqs_k, u0k, (0.0, 20.0), ps_k; backend = ShardedBackend(), jac = mode)
    sol = solve(prob, Rodas5P(); abstol = 1.0e-9, reltol = 1.0e-9, saveat = 5.0)
    ref_k = solve(
        ODEProblem(eqs_k, u0k, (0.0, 20.0), ps_k; backend = ShardedBackend()), RK4();
        abstol = 1.0e-10, reltol = 1.0e-10, saveat = 5.0,
    )
    @test sol.retcode == SciMLBase.ReturnCode.Success
    @test maximum(abs, sol.u[end] .- ref_k.u[end]) < 1.0e-5
end

@testset "Jacobian tracks update_parameters! ($mode)" for mode in (:fd, :analytic)
    prob = ODEProblem(eqs, u0, (0.0, 1.0), ps; backend = ShardedBackend(), jac = mode)
    pd2 = Dict(J => 1.7, hx => 0.4, γ => 0.31)
    update_parameters!(prob, pd2)
    fresh = ODEProblem(eqs, u0, (0.0, 1.0), pd2; backend = ShardedBackend(), jac = mode)
    Ja = copy(prob.f.jac_prototype)
    Jb = copy(fresh.f.jac_prototype)
    prob.f.jac(Ja, u, prob.p, 0.0)
    fresh.f.jac(Jb, u, fresh.p, 0.0)
    @test maximum(abs, Array(Ja) .- Array(Jb)) < 1.0e-10
end

@testset "invalid jac mode errors" begin
    @test_throws ArgumentError ODEProblem(
        eqs, u0, (0.0, 1.0), ps; backend = ShardedBackend(), jac = :bogus
    )
end
