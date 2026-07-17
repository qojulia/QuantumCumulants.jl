using QuantumCumulants
using SymbolicUtils: SymbolicUtils
using Symbolics: Symbolics, @variables
using SciMLBase: SciMLBase, ODEProblem, ODEFunction
using ModelingToolkitBase: mtkcompile, unknowns
using OrdinaryDiffEqLowOrderRK: RK4, solve
using QuantumOpticsBase: SpinBasis, spinup, tensor
using Test

# Reference RHS by direct substitution into the completed equations, built from the
# public surface only (equations, states, average, adjoint, undo_average).
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

# Transverse-field Ising chain of 3 Pauli spins at order 2: 36 coupled moment equations.
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

@testset "du vs substitution reference (Pauli chain)" begin
    prob = ODEProblem(eqs, u0, (0.0, 1.0), ps; backend = KernelBackend())
    du = similar(u)
    prob.f(du, u, prob.p, 0.0)
    ref = reference_du(eqs, ps, u)
    @test maximum(abs.(du .- ref) ./ max.(abs.(ref), 1.0e-12)) < 1.0e-12
end

@testset "du vs substitution reference (Kerr cavity, conj-folded)" begin
    hc = FockSpace(:cavity)
    a = Destroy(hc, :a)
    @variables Δ Ω κ U
    Hk = Δ * a' * a + U * a' * a' * a * a + Ω * (a + a')
    eqs_k = meanfield([a], Hk, [a]; rates = [κ], order = 2)
    # folded closure: one representative per conjugate pair, so RHS leaves reference
    # conj(state) and the reference substitution below stays well defined off-trajectory
    complete!(eqs_k; get_adjoints = false)
    ps_k = Dict(Δ => -1.0, Ω => 1.3, κ => 1.0, U => 0.1)
    nk = length(eqs_k.states)
    uk = ComplexF64[0.3cos(2.1i) + 0.2im * sin(0.7i) for i in 1:nk]
    prob = ODEProblem(eqs_k, zeros(ComplexF64, nk), (0.0, 1.0), ps_k; backend = KernelBackend())
    du = similar(uk)
    prob.f(du, uk, prob.p, 0.0)
    ref = reference_du(eqs_k, ps_k, uk)
    @test maximum(abs.(du .- ref) ./ max.(abs.(ref), 1.0e-12)) < 1.0e-12
end

ψ0 = tensor([spinup(SpinBasis(1 // 2)) for _ in 1:Np]...)

@testset "trajectory vs the MTK path" begin
    prob = ODEProblem(eqs, ψ0, (0.0, 5.0), ps; backend = KernelBackend())
    sol = solve(prob, RK4(); abstol = 1.0e-10, reltol = 1.0e-10, saveat = 0.5)
    sys = mtkcompile(System(eqs; name = :sys))
    u0_mtk = Dict(unknowns(sys) .=> initial_values(eqs, ψ0))
    prob_mtk = ODEProblem(sys, merge(u0_mtk, Dict(collect(ps))), (0.0, 5.0))
    sol_mtk = solve(prob_mtk, RK4(); abstol = 1.0e-10, reltol = 1.0e-10, saveat = 0.5)
    @test maximum(maximum(abs, a .- b) for (a, b) in zip(sol.u, sol_mtk.u)) < 1.0e-6
end

@testset "construction is deterministic (bit-exact du across builds)" begin
    prob_a = ODEProblem(eqs, u0, (0.0, 1.0), ps; backend = KernelBackend())
    prob_b = ODEProblem(eqs, u0, (0.0, 1.0), ps; backend = KernelBackend())
    du_a, du_b = similar(u), similar(u)
    prob_a.f(du_a, u, prob_a.p, 0.0)
    prob_b.f(du_b, u, prob_b.p, 0.0)
    @test du_a == du_b
end

# The RHS threads both passes (the flat monomial update and the SpMV) through Polyester's
# `@batch` when `parallel` is set. CI is single-threaded, so forcing `parallel = true`
# exercises the threaded code path (Polyester runs serially on one thread, bit-identical to
# the prefix update and CSC-order gather) without depending on real concurrency.
@testset "parallel RHS matches serial (bit-exact)" begin
    ser = ODEProblem(eqs, u0, (0.0, 1.0), ps; backend = KernelBackend(parallel = false))
    par = ODEProblem(eqs, u0, (0.0, 1.0), ps; backend = KernelBackend(parallel = true))
    @test par.f.f.kernel.parallel
    @test !ser.f.f.kernel.parallel
    du_ser, du_par = similar(u), similar(u)
    ser.f(du_ser, u, ser.p, 0.0)
    par.f(du_par, u, par.p, 0.0)
    @test du_par == du_ser
    # small system (36 eqs < KERNEL_PARALLEL_MIN): :auto stays serial
    @test !ODEProblem(eqs, u0, (0.0, 1.0), ps; backend = KernelBackend()).f.f.kernel.parallel
end

@testset "u0 input forms agree" begin
    prob_state = ODEProblem(eqs, ψ0, (0.0, 1.0), ps; backend = KernelBackend())
    vec0 = initial_values(eqs, ψ0)
    prob_vec = ODEProblem(eqs, vec0, (0.0, 1.0), ps; backend = KernelBackend())
    dict0 = Dict(eqs.states[k] => vec0[k] for k in eachindex(vec0))
    prob_dict = ODEProblem(eqs, dict0, (0.0, 1.0), ps; backend = KernelBackend())
    @test prob_state.u0 == prob_vec.u0
    @test prob_dict.u0 == prob_vec.u0
end
