using QuantumCumulants
using SymbolicUtils: SymbolicUtils
using Symbolics: Symbolics, @variables
using SciMLBase: SciMLBase, ODEProblem
using OrdinaryDiffEqLowOrderRK: RK4, solve
using OrdinaryDiffEqRosenbrock: OrdinaryDiffEqRosenbrock, Rodas5P
using QuantumOpticsBase: SpinBasis, spinup, tensor
using Random: MersenneTwister
using Test

# Substitution dict at a numeric point (params + states + conjugate partners + the
# algebra's symbolic imaginary unit), on the public surface only.
function build_subs(eqs, pdict, u)
    subs = Dict{Any, Any}(Symbolics.unwrap(k) => v for (k, v) in pdict)
    for (i, s) in enumerate(eqs.states)
        subs[Symbolics.unwrap(s)] = u[i]
    end
    for (i, s) in enumerate(eqs.states)
        sa = Symbolics.unwrap(average(adjoint(undo_average(s))))
        haskey(subs, sa) || (subs[sa] = conj(u[i]))
    end
    for eq in eqs.equations, v in Symbolics.get_variables(Symbolics.unwrap(eq.rhs))
        uu = SymbolicUtils.unwrap(v)
        SymbolicUtils.issym(uu) && SymbolicUtils.nameof(uu) === :im && (subs[uu] = im)
    end
    return subs
end

# Transverse-field Ising chain of 3 Pauli spins at order 2 (unfolded closure, so the
# drift is holomorphic in the states and the analytic Jacobian exists).
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

@testset "analytic Jacobian vs Symbolics.derivative" begin
    prob = ODEProblem(eqs, u0, (0.0, 1.0), ps; backend = KernelBackend(), jac = true)
    Js = copy(prob.f.jac_prototype)
    prob.f.jac(Js, u, prob.p, 0.0)
    subs = build_subs(eqs, ps, u)
    rng = MersenneTwister(1)
    maxrel = 0.0
    for _ in 1:25
        i, j = rand(rng, 1:nst), rand(rng, 1:nst)
        d = Symbolics.derivative(
            eqs.equations[i].rhs, Symbolics.wrap(Symbolics.unwrap(eqs.states[j]))
        )
        ref = ComplexF64(SymbolicUtils.unwrap_const(Symbolics.substitute(Symbolics.unwrap(d), subs)))
        maxrel = max(maxrel, abs(Js[i, j] - ref) / max(abs(ref), 1.0e-12))
    end
    @test maxrel < 1.0e-10
end

ψ0 = tensor([spinup(SpinBasis(1 // 2)) for _ in 1:Np]...)

@testset "implicit solve with analytic J" begin
    prob = ODEProblem(eqs, ψ0, (0.0, 5.0), ps; jac = true)
    # finite-difference fallback for the W factorization internals; the supplied
    # analytic jac is what actually fills J (ForwardDiff rejects complex state)
    sol = solve(
        prob, Rodas5P(autodiff = OrdinaryDiffEqRosenbrock.AutoFiniteDiff());
        saveat = 0.5,
    )
    @test sol.retcode == SciMLBase.ReturnCode.Success
    sol_rk = solve(
        ODEProblem(eqs, ψ0, (0.0, 5.0), ps), RK4();
        abstol = 1.0e-10, reltol = 1.0e-10, saveat = 0.5,
    )
    @test maximum(abs, sol.u[end] .- sol_rk.u[end]) < 1.0e-4
end

@testset "jac = true survives update_parameters!" begin
    prob = ODEProblem(eqs, u0, (0.0, 1.0), ps; backend = KernelBackend(), jac = true)
    pd2 = Dict(J => 1.7, hx => 0.4, γ => 0.31)
    update_parameters!(prob, pd2)
    fresh = ODEProblem(eqs, u0, (0.0, 1.0), pd2; backend = KernelBackend(), jac = true)
    J_a = copy(prob.f.jac_prototype)
    J_b = copy(fresh.f.jac_prototype)
    prob.f.jac(J_a, u, prob.p, 0.0)
    fresh.f.jac(J_b, u, fresh.p, 0.0)
    @test J_a == J_b
    du_a, du_b = similar(u), similar(u)
    prob.f(du_a, u, prob.p, 0.0)
    fresh.f(du_b, u, fresh.p, 0.0)
    @test du_a == du_b
end

@testset "typed errors" begin
    # conj-folded Kerr: the drift references conj(states), not holomorphic
    hc = FockSpace(:cavity)
    a = Destroy(hc, :a)
    @variables Δ Ω κ U
    Hk = Δ * a' * a + U * a' * a' * a * a + Ω * (a + a')
    eqs_k = meanfield([a], Hk, [a]; rates = [κ], order = 2)
    complete!(eqs_k; get_adjoints = false)
    ps_k = Dict(Δ => -1.0, Ω => 1.3, κ => 1.0, U => 0.1)
    @test_throws HolomorphicJacobianError ODEProblem(
        eqs_k, zeros(ComplexF64, length(eqs_k.states)), (0.0, 1.0), ps_k;
        backend = KernelBackend(), jac = true,
    )
    # sharded has no analytic Jacobian
    @test_throws ArgumentError ODEProblem(
        eqs, u0, (0.0, 1.0), ps; backend = ShardedBackend(), jac = true
    )
end
