using QuantumCumulants
using SymbolicUtils: SymbolicUtils
using Symbolics: Symbolics, @variables
using SciMLBase: SciMLBase, ODEProblem
using OrdinaryDiffEqLowOrderRK: RK4, solve
using Test

# Reference RHS by direct substitution (same helper as kernel_backend_test.jl; test
# files run standalone under ParallelTestRunner).
function reference_du(eqs, pdict, u)
    subs = Dict{Any, Any}(Symbolics.unwrap(k) => v for (k, v) in pdict)
    for (i, s) in enumerate(eqs.states)
        subs[Symbolics.unwrap(s)] = u[i]
    end
    for (i, s) in enumerate(eqs.states)
        sa = Symbolics.unwrap(average(adjoint(undo_average(s))))
        haskey(subs, sa) || (subs[sa] = conj(u[i]))
    end
    ref = Vector{ComplexF64}(undef, length(eqs.states))
    for (i, eq) in enumerate(eqs.equations)
        rhs = Symbolics.unwrap(eq.rhs)
        for v in Symbolics.get_variables(rhs)
            uu = SymbolicUtils.unwrap(v)
            SymbolicUtils.issym(uu) && SymbolicUtils.nameof(uu) === :im && (subs[uu] = im)
        end
        ref[i] = ComplexF64(SymbolicUtils.unwrap_const(Symbolics.substitute(rhs, subs)))
    end
    return ref
end

# Transverse-field Ising chain of 3 Pauli spins at order 2 (36 equations).
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
u0 = zeros(ComplexF64, nst)
u = ComplexF64[0.1cos(3.7i) + 0.05im * sin(1.3i) for i in 1:nst]

@testset "sharded du and trajectory" begin
    pk = ODEProblem(eqs, u0, (0.0, 5.0), ps; backend = KernelBackend())
    psh = ODEProblem(eqs, u0, (0.0, 5.0), ps; backend = ShardedBackend(chunk = 5))
    du_k, du_s = similar(u0), similar(u0)
    pk.f(du_k, u, pk.p, 0.0)
    psh.f(du_s, u, psh.p, 0.0)
    # two independent implementations of the same resolved drift agree
    @test maximum(abs.(du_k .- du_s)) < 1.0e-13
    ref = reference_du(eqs, ps, u)
    @test maximum(abs.(du_s .- ref) ./ max.(abs.(du_s), 1.0e-12)) < 1.0e-12
    sol_k = solve(pk, RK4(); abstol = 1.0e-10, reltol = 1.0e-10, saveat = 0.5)
    sol_s = solve(psh, RK4(); abstol = 1.0e-10, reltol = 1.0e-10, saveat = 0.5)
    @test maximum(maximum(abs, a .- b) for (a, b) in zip(sol_k.u, sol_s.u)) < 1.0e-6
end

# Two-spin order-1 fixture for the t-dependent drift (the kernel rejects it).
Ns = 2
hs = ⊗([PauliSpace(Symbol(:s, i)) for i in 1:Ns]...)
sz(i) = Pauli(hs, :σ, 3, i)
sx(i) = Pauli(hs, :σ, 1, i)
sy(i) = Pauli(hs, :σ, 2, i)
sm(i) = (sx(i) - 1im * sy(i)) / 2
Hs = -J * sz(1) * sz(2) - hx * (sx(1) + sx(2))
eqs_s = meanfield([sz(i) for i in 1:Ns], Hs, [sm(i) for i in 1:Ns]; rates = [γ, γ], order = 1)
complete!(eqs_s)
u0s = zeros(ComplexF64, length(eqs_s.states))

@testset "sharded accepts what the kernel rejects" begin
    eqs_t = modify_equations(eqs_s, (op, d) -> cos(eqs_s.iv) * d)      # t-dependent drift
    prob = ODEProblem(eqs_t, u0s, (0.0, 1.0), ps; backend = ShardedBackend())
    @test solve(prob, RK4()).retcode == SciMLBase.ReturnCode.Success
end

@testset "sharded update_parameters! and dual guard" begin
    psh = ODEProblem(eqs, u0, (0.0, 5.0), ps; backend = ShardedBackend())
    du_k, du_s = similar(u0), similar(u0)
    update_parameters!(psh, pd2)
    fresh = ODEProblem(eqs, u0, (0.0, 5.0), pd2; backend = ShardedBackend())
    psh.f(du_s, u, psh.p, 0.0)
    fresh.f(du_k, u, fresh.p, 0.0)
    @test du_s == du_k
    @test_throws ArgumentError psh.f(du_s, zeros(Float64, length(u0)), psh.p, 0.0)
end
