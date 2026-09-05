using QuantumCumulants
using SymbolicUtils: SymbolicUtils
using Symbolics: Symbolics, @variables
using SciMLBase: SciMLBase, ODEProblem
using ModelingToolkitBase: mtkcompile, unknowns
using OrdinaryDiffEqLowOrderRK: RK4, solve
using Test

"""Max deviation over all states and save points between a backend and MTK trajectory."""
function traj_vs_mtk(eqs, ps, backend, tspan; saveat)
    n = length(eqs.states)
    u0 = zeros(ComplexF64, n)
    prob = ODEProblem(eqs, u0, tspan, ps; backend)
    sol = solve(prob, RK4(); saveat, abstol = 1.0e-10, reltol = 1.0e-10)
    sys = mtkcompile(System(eqs; name = :sys))
    prob_m =
        ODEProblem(sys, merge(initial_values(eqs, u0), Dict(parameter_map(eqs, ps))), tspan)
    sol_m = solve(prob_m, RK4(); saveat, abstol = 1.0e-10, reltol = 1.0e-10)
    maxdev = 0.0
    for (i, st) in enumerate(eqs.states)
        ref = get_solution(sol_m, SymbolicUtils.unwrap(st), eqs).(sol.t)
        mine = [sol.u[j][i] for j in eachindex(sol.t)]
        maxdev = max(maxdev, maximum(abs.(ref .- mine)))
    end
    return maxdev
end

# ---- superradiant laser (order 2, phase-invariant completion) --------------------------

hc = FockSpace(:cavity)
ha = NLevelSpace(:atom, 2)
h = hc ⊗ ha
@qnumbers a::Destroy(h)
σ(α, β, k) = IndexedOperator(Transition(h, :σ, α, β), k)
@variables N Δ κ Γ R ν
g(k) = IndexedVariable(:g, k)
i = Index(h, :i, N, ha)
j = Index(h, :j, N, ha)
H = -Δ * a'a + Σ(g(i) * (a' * σ(1, 2, i) + a * σ(2, 1, i)), i)
Jops = [a, σ(1, 2, i), σ(2, 1, i), σ(2, 2, i)]
eqs_sr = meanfield([a' * a, σ(2, 2, j)], H, Jops; rates = [κ, Γ, R, ν], order = 2)
eqs_src = complete(eqs_sr; filter_func = phase_invariant)
pvals = Dict(Δ => 2.0, κ => 8.0, Γ => 1.0, R => 2.0, ν => 1.0)

Nev = 4
eqs_ev = evaluate(eqs_src; limits = (N => Nev))
gvals = [1.5 - 0.1k for k = 1:Nev]
ps_ev = merge(pvals, Dict(g(i) => gvals))

@testset "evaluated superradiant laser (1D array parameter)" begin
    @test traj_vs_mtk(eqs_ev, ps_ev, KernelBackend(), (0.0, 2.0); saveat = 0.2) < 1.0e-6
end

@testset "array-parameter update_parameters!" begin
    g2vals = [0.9 + 0.05k for k = 1:Nev]
    ps_new = merge(pvals, Dict(g(i) => g2vals))
    n = length(eqs_ev.states)
    u0 = zeros(ComplexF64, n)
    u = ComplexF64[0.1cos(3.7k) + 0.05im * sin(1.3k) for k = 1:n]
    du_a, du_b = similar(u0), similar(u0)
    prob = ODEProblem(eqs_ev, u0, (0.0, 1.0), ps_ev; backend = KernelBackend())
    update_parameters!(prob, Dict(g(i) => g2vals))
    fresh = ODEProblem(eqs_ev, u0, (0.0, 1.0), ps_new; backend = KernelBackend())
    prob.f(du_a, u, prob.p, 0.0)
    fresh.f(du_b, u, fresh.p, 0.0)
    @test du_a == du_b
end

@testset "scaled superradiant laser (scalar g)" begin
    eqs_sc = scale(eqs_src)
    ps_sc = merge(pvals, Dict(N => 20.0, g(i) => 1.5))
    @test traj_vs_mtk(eqs_sc, ps_sc, KernelBackend(), (0.0, 2.0); saveat = 0.2) < 1.0e-6
end

# ---- evaluated cavity antiresonance (2D array parameters) ------------------------------

@testset "evaluated cavity antiresonance (2D array parameters)" begin
    hc2 = FockSpace(:cavity)
    ha2 = NLevelSpace(:atom2, 2)
    h2 = hc2 ⊗ ha2
    @variables N2 Δc η Δa κ2
    g2(k) = IndexedVariable(:g2, k)
    Γ2(k, l) = DoubleIndexedVariable(:Γ2, k, l)
    Ω2(k, l) = DoubleIndexedVariable(:Ω2, k, l; identical = false)
    i2 = Index(h2, :i, N2, ha2)
    j2 = Index(h2, :j, N2, ha2)
    @qnumbers b::Destroy(h2)
    s(x, y, k) = IndexedOperator(Transition(h2, :σ, x, y), k)
    Hc2 = Δc * b'b + η * (b' + b)
    Ha2 =
        Δa * Σ(s(2, 2, i2), i2) + Σ(Σ(Ω2(i2, j2) * s(2, 1, i2) * s(1, 2, j2), j2, [i2]), i2)
    Hi2 = Σ(g2(i2) * (b' * s(1, 2, i2) + b * s(2, 1, i2)), i2)
    eqs_a =
        meanfield(b, Hc2 + Ha2 + Hi2, [b, s(1, 2, i2)]; rates = [κ2, Γ2(i2, j2)], order = 1)
    complete!(eqs_a)
    Na = 2
    eqs_ae = evaluate(eqs_a; limits = (N2 => Na))
    Γm = [k == l ? 1.0 : 0.4 for k = 1:Na, l = 1:Na]
    Ωm = [k == l ? 0.0 : 1.1 for k = 1:Na, l = 1:Na]
    ps_a = Dict(
        Δc => 1.0,
        η => 0.2,
        Δa => 0.5,
        κ2 => 20.0,
        g2(i2) => [2.0, -2.0],
        Γ2(i2, j2) => Γm,
        Ω2(i2, j2) => Ωm,
    )
    @test traj_vs_mtk(eqs_ae, ps_a, KernelBackend(), (0.0, 5.0); saveat = 0.5) < 1.0e-6

    # 2D array update: rescale Γ2, compare against a fresh construction, bit-exact
    Γm2 = 1.5 .* Γm
    ps_a2 = merge(ps_a, Dict(Γ2(i2, j2) => Γm2))
    n = length(eqs_ae.states)
    u0 = zeros(ComplexF64, n)
    u = ComplexF64[0.1cos(3.7k) + 0.05im * sin(1.3k) for k = 1:n]
    du_a, du_b = similar(u0), similar(u0)
    prob = ODEProblem(eqs_ae, u0, (0.0, 1.0), ps_a; backend = KernelBackend())
    update_parameters!(prob, Dict(Γ2(i2, j2) => Γm2))
    fresh = ODEProblem(eqs_ae, u0, (0.0, 1.0), ps_a2; backend = KernelBackend())
    prob.f(du_a, u, prob.p, 0.0)
    fresh.f(du_b, u, fresh.p, 0.0)
    @test du_a == du_b
end
