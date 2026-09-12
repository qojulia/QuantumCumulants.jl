using QuantumCumulants
using ModelingToolkitBase: mtkcompile
using OrdinaryDiffEqLowOrderRK: RK4, solve
using SciMLBase: ODEProblem
using Symbolics: @variables
using Test

function array_traj_vs_mtk(eqs, ps, tspan; saveat)
    u0 = zeros(ComplexF64, length(eqs.states))
    direct = solve(
        ODEProblem(eqs, u0, tspan, ps; backend = KernelBackend()),
        RK4();
        saveat,
        abstol = 1.0e-10,
        reltol = 1.0e-10,
    )

    sys = mtkcompile(System(eqs; name = :array_reference))
    values = merge(initial_values(eqs, u0), Dict(parameter_map(eqs, ps)))
    reference = solve(
        ODEProblem(sys, parameter_map(sys, values), tspan),
        RK4();
        saveat,
        abstol = 1.0e-10,
        reltol = 1.0e-10,
    )

    maxdev = 0.0
    for state in eqs.states
        got = get_solution(direct, state, eqs).(direct.t)
        ref = get_solution(reference, state, eqs).(reference.t)
        maxdev = max(maxdev, maximum(abs.(got .- ref)))
    end
    return maxdev
end

@testset "evaluated 1D indexed parameter" begin
    hc = FockSpace(:array_cavity)
    ha = NLevelSpace(:array_atom, 2)
    h = hc ⊗ ha
    @qnumbers a::Destroy(h)
    σ(α, β, k) = IndexedOperator(Transition(h, :σ, α, β), k)
    @variables N Δ κ Γ R ν
    g(k) = IndexedVariable(:g, k)
    i = Index(h, :i, N, ha)
    j = Index(h, :j, N, ha)

    H = -Δ * a'a + Σ(g(i) * (a' * σ(1, 2, i) + a * σ(2, 1, i)), i)
    jumps = [a, σ(1, 2, i), σ(2, 1, i), σ(2, 2, i)]
    eqs = meanfield([a' * a, σ(2, 2, j)], H, jumps; rates = [κ, Γ, R, ν], order = 2)
    closed = complete(eqs; filter_func = phase_invariant)
    evaluated = evaluate(closed; limits = (N => 4))

    gvals = [1.5 - 0.1k for k in 1:4]
    ps = Dict(
        Δ => 2.0,
        κ => 8.0,
        Γ => 1.0,
        R => 2.0,
        ν => 1.0,
        g(i) => gvals,
    )
    @test array_traj_vs_mtk(evaluated, ps, (0.0, 0.5); saveat = 0.1) < 1.0e-6

    changed = [0.9 + 0.05k for k in 1:4]
    u0 = zeros(ComplexF64, length(evaluated.states))
    u = ComplexF64[
        0.1cos(3.7k) + 0.05im * sin(1.3k) for k in eachindex(evaluated.states)
    ]
    prob = ODEProblem(evaluated, u0, (0.0, 1.0), ps; backend = KernelBackend())
    update_parameters!(prob, Dict(g(i) => changed))
    fresh = ODEProblem(
        evaluated,
        u0,
        (0.0, 1.0),
        merge(ps, Dict(g(i) => changed));
        backend = KernelBackend(),
    )
    du_a = similar(u)
    du_b = similar(u)
    prob.f(du_a, u, prob.p, 0.0)
    fresh.f(du_b, u, fresh.p, 0.0)
    @test du_a == du_b
end

@testset "evaluated 2D indexed parameters" begin
    hc = FockSpace(:array_cavity_2d)
    ha = NLevelSpace(:array_atom_2d, 2)
    h = hc ⊗ ha
    @qnumbers a::Destroy(h)
    σ(x, y, k) = IndexedOperator(Transition(h, :σ, x, y), k)
    @variables N Δc η Δa κ
    g(k) = IndexedVariable(:g2, k)
    Γ(k, l) = DoubleIndexedVariable(:Γ2, k, l)
    Ω(k, l) = DoubleIndexedVariable(:Ω2, k, l; identical = false)
    i = Index(h, :i, N, ha)
    j = Index(h, :j, N, ha)

    Hc = Δc * a'a + η * (a' + a)
    Ha =
        Δa * Σ(σ(2, 2, i), i) +
        Σ(Σ(Ω(i, j) * σ(2, 1, i) * σ(1, 2, j), j, [i]), i)
    Hi = Σ(g(i) * (a' * σ(1, 2, i) + a * σ(2, 1, i)), i)
    eqs = meanfield(a, Hc + Ha + Hi, [a, σ(1, 2, i)]; rates = [κ, Γ(i, j)], order = 1)
    complete!(eqs)
    evaluated = evaluate(eqs; limits = (N => 2))

    Γm = [k == l ? 1.0 : 0.4 for k in 1:2, l in 1:2]
    Ωm = [k == l ? 0.0 : 1.1 for k in 1:2, l in 1:2]
    ps = Dict(
        Δc => 1.0,
        η => 0.2,
        Δa => 0.5,
        κ => 20.0,
        g(i) => [2.0, -2.0],
        Γ(i, j) => Γm,
        Ω(i, j) => Ωm,
    )

    u0 = zeros(ComplexF64, length(evaluated.states))
    u = ComplexF64[
        0.1cos(3.7k) + 0.05im * sin(1.3k) for k in eachindex(evaluated.states)
    ]
    changed = 1.5 .* Γm
    prob = ODEProblem(evaluated, u0, (0.0, 1.0), ps; backend = KernelBackend())
    update_parameters!(prob, Dict(Γ(i, j) => changed))
    fresh = ODEProblem(
        evaluated,
        u0,
        (0.0, 1.0),
        merge(ps, Dict(Γ(i, j) => changed));
        backend = KernelBackend(),
    )

    du_a = similar(u)
    du_b = similar(u)
    prob.f(du_a, u, prob.p, 0.0)
    fresh.f(du_b, u, fresh.p, 0.0)
    @test du_a == du_b
end
