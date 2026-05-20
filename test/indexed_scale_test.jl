using QuantumCumulants
using Symbolics: Symbolics, @variables
using ModelingToolkitBase: @named, mtkcompile, ODEProblem, unknowns
using OrdinaryDiffEq: Tsit5, solve, ReturnCode
using Test

# v1 surface: indexed `scale!` over a permutation-symmetric atom subspace.
# Master's `scale(eqs; h = [k])` per-Hilbert-space variant is not yet
# implemented in v1; see test/pending/README.md and TODO.md. The full-
# system scale (no `h` kwarg) is exercised below.

@testset "indexed_scale: Tavis-Cummings closure + scale + ODE" begin
    @variables N::Real Δ::Real g::Real κ::Real Γ::Real R::Real ν::Real

    hc = FockSpace(:cavity)
    ha = NLevelSpace(:atom, 2)
    h = hc ⊗ ha

    k = Index(h, :k, N, ha)
    @qnumbers a::Destroy(h)
    σ(i, j, k) = IndexedOperator(Transition(h, :σ, i, j), k)

    H = -Δ * a' * a + g * (Σ(a' * σ(1, 2, k), k) + Σ(a * σ(2, 1, k), k))
    J = [a, σ(1, 2, k), σ(2, 1, k), σ(2, 2, k)]
    rates = [κ, Γ, R, ν]
    ops = [a' * a, σ(2, 2, k)]

    eqs = meanfield(ops, H, J; rates = rates, order = 1)
    eqs_c = complete(eqs)
    @test length(eqs_c.equations) >= 2

    # Scale collapses the permutation-symmetric atom subspace.
    eqs_sc = scale(eqs_c)
    @test length(eqs_sc.equations) >= 1
    @test isempty(find_missing(eqs_sc; get_adjoints = false))

    @named sys = to_system(eqs_sc)
    sys_c = mtkcompile(sys)
    @test length(unknowns(sys_c)) >= 1
end

@testset "indexed_scale: per-Hilbert-space scale order independence" begin
    # Two indexed subspaces, scaling each separately in either order must
    # produce systems of equal size. From the pending port of master's
    # `test_indexed_scale.jl`.
    @variables N::Real N2::Real Δ::Real g::Real κ::Real Γ::Real R::Real ν::Real

    hc = FockSpace(:cavity)
    ha = NLevelSpace(:atom, 2)
    h = hc ⊗ ha

    k = Index(h, :k, N, ha)
    l = Index(h, :l, N, ha)
    m = Index(h, :m, N2, hc)
    n = Index(h, :n, N2, hc)

    σ(i, j, k) = IndexedOperator(Transition(h, :σ, i, j), k)
    ai(k) = IndexedOperator(Destroy(h, :a), k)

    H_2 = -Δ * ∑(ai(m)' * ai(m), m) +
        g * (∑(Σ(ai(m)' * σ(1, 2, k), k), m) + ∑(Σ(ai(m) * σ(2, 1, k), k), m))
    J_2 = [ai(m), σ(1, 2, k), σ(2, 1, k), σ(2, 2, k)]
    rates_2 = [κ, Γ, R, ν]
    ops_2 = [ai(n)' * ai(n), σ(2, 2, l)]
    eqs_com = complete(meanfield(ops_2, H_2, J_2; rates = rates_2, order = 2))

    s_full_a = scale(scale(eqs_com; h = [1]); h = [2])
    s_full_b = scale(scale(eqs_com; h = [2]); h = [1])
    @test length(s_full_a.equations) == length(s_full_b.equations)
end

@testset "indexed_scale: Ising-XX scale vs evaluate ODE numerical equality" begin
    h2 = NLevelSpace(:spin, 2)
    @variables V::Real Ω::Real N3::Real
    i1 = Index(h2, :i1, N3, h2)
    i2 = Index(h2, :i2, N3, h2)
    i = Index(h2, :i, N3, h2)
    s(α, β, idx) = IndexedOperator(Transition(h2, :S, α, β), idx)
    sp(idx) = s(2, 1, idx); sm(idx) = s(1, 2, idx)
    int_sum = Σ(sp(i1) * sm(i2) + sm(i1) * sp(i2), i1, i2) -
        Σ(sp(i1) * sm(i1) + sm(i1) * sp(i1), i1)
    Hint = V * int_sum + Ω * Σ(sp(i1) + sm(i1), i1)
    eqs_c = complete(meanfield([s(1, 2, i)], Hint; order = 1))

    eqs_sc = scale(eqs_c)
    eqs_ev = evaluate(eqs_c; limits = (N3 => 3))
    @test isempty(find_missing(eqs_sc; get_adjoints = false))
    @test isempty(find_missing(eqs_ev; get_adjoints = false))

    @named sys_sc = to_system(eqs_sc); sys_sc_c = mtkcompile(sys_sc)
    @named sys_ev = to_system(eqs_ev); sys_ev_c = mtkcompile(sys_ev)

    N_, V_, Ω_ = 3, 4.79 / 2, 1.0
    u0_sc = Dict(unknowns(sys_sc_c) .=> zeros(ComplexF64, length(unknowns(sys_sc_c))))
    u0_ev = Dict(unknowns(sys_ev_c) .=> zeros(ComplexF64, length(unknowns(sys_ev_c))))

    prob_sc = ODEProblem(sys_sc_c, merge(u0_sc, Dict([N3, V, Ω] .=> [N_, V_, Ω_])), (0.0, 2.0))
    # post-evaluate, N3 baked in; only V, Ω remain symbolic.
    prob_ev = ODEProblem(sys_ev_c, merge(u0_ev, Dict([V, Ω] .=> [V_, Ω_])), (0.0, 2.0))
    sol_sc = solve(prob_sc, Tsit5(); abstol = 1.0e-10, reltol = 1.0e-10)
    sol_ev = solve(prob_ev, Tsit5(); abstol = 1.0e-10, reltol = 1.0e-10)
    @test sol_sc.retcode == ReturnCode.Success
    @test sol_ev.retcode == ReturnCode.Success

    # Master's actual assertion: the two derivation paths must produce
    # the same physics. Compare the observable ⟨s_12⟩ at the same time
    # points; the canonical state lives at atom position 1 in both cases.
    SQA = QuantumCumulants.SecondQuantizedAlgebra
    # `eqs_sc` has the LHS state ⟨s_(canon)₁₂⟩; the underlying op uses the
    # canonical first-declared atom index.
    obs_op = SQA.undo_average(eqs_sc.states[1])
    val_sc_end = get_solution(sol_sc, obs_op, eqs_sc)(sol_sc.t[end])
    # Evaluate produces ⟨s_(i_1)₁₂⟩, ⟨s_(i_2)₁₂⟩, ⟨s_(i_3)₁₂⟩; by symmetry
    # all three should equal val_sc_end up to numerical noise.
    obs_op_ev_1 = SQA.undo_average(eqs_ev.states[1])
    obs_op_ev_2 = SQA.undo_average(eqs_ev.states[2])
    obs_op_ev_3 = SQA.undo_average(eqs_ev.states[3])
    val_ev_1 = get_solution(sol_ev, obs_op_ev_1, eqs_ev)(sol_ev.t[end])
    val_ev_2 = get_solution(sol_ev, obs_op_ev_2, eqs_ev)(sol_ev.t[end])
    val_ev_3 = get_solution(sol_ev, obs_op_ev_3, eqs_ev)(sol_ev.t[end])
    @test isapprox(val_ev_1, val_sc_end; atol = 1.0e-6)
    @test isapprox(val_ev_2, val_sc_end; atol = 1.0e-6)
    @test isapprox(val_ev_3, val_sc_end; atol = 1.0e-6)
    # Permutation symmetry: the three atoms must give the same value.
    @test isapprox(val_ev_1, val_ev_2; atol = 1.0e-8)
    @test isapprox(val_ev_2, val_ev_3; atol = 1.0e-8)
end
