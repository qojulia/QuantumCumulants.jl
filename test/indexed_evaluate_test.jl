using QuantumCumulants
using Symbolics: Symbolics, @variables
using ModelingToolkitBase: @named, mtkcompile, ODEProblem, unknowns
using OrdinaryDiffEq: Tsit5, solve, ReturnCode
using Test

@testset "evaluate: indexed Ising-XX (N=3) → mtkcompile" begin
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
    eqs = meanfield([s(1, 2, i)], Hint; order = 1)
    eqs_c = complete(eqs)

    @test length(eqs_c.equations) == 2

    # At N=3, each symbolic state with one free atom-space index enumerates
    # into three concrete-position equations.
    eqs_ev = evaluate(eqs_c; limits = (N3 => 3))
    @test length(eqs_ev.equations) == 6
    @test isempty(find_missing(eqs_ev; get_adjoints = false))

    @named sys = System(eqs_ev)
    sys_c = mtkcompile(sys)
    @test length(unknowns(sys_c)) == 6

    N_, V_, Ω_ = 3, 4.79 / 2, 1.0
    u0 = Dict(unknowns(sys_c) .=> zeros(ComplexF64, length(unknowns(sys_c))))
    prob = ODEProblem(
        sys_c, merge(u0, Dict([V, Ω] .=> [V_, Ω_])), (0.0, 2.0)
    )
    sol = solve(prob, Tsit5())
    @test sol.retcode == ReturnCode.Success

    # Physicality: single-atom coherence magnitudes are bounded by 1
    # (|⟨σ_12⟩|² ≤ ⟨σ_21 σ_12⟩ = ⟨σ_22⟩ ≤ 1 by Cauchy-Schwarz).
    SQA = QuantumCumulants.SecondQuantizedAlgebra
    for op_avg in eqs_ev.states
        op = SQA.undo_average(op_avg)
        mag = [abs(get_solution(sol, op, eqs_ev)(t)) for t in sol.t]
        @test maximum(mag) <= 1.0 + 1.0e-9
    end
end

@testset "evaluate: scale agreement at the steady-state limit" begin
    # Scale-first and evaluate-first give the same dynamics on a
    # permutation-symmetric system.
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

    N_, V_, Ω_ = 3, 4.79 / 2, 1.0
    @named sys_sc = System(eqs_sc); sys_sc_c = mtkcompile(sys_sc)
    @named sys_ev = System(eqs_ev); sys_ev_c = mtkcompile(sys_ev)
    u0_sc = Dict(unknowns(sys_sc_c) .=> zeros(ComplexF64, length(unknowns(sys_sc_c))))
    u0_ev = Dict(unknowns(sys_ev_c) .=> zeros(ComplexF64, length(unknowns(sys_ev_c))))

    prob_sc = ODEProblem(sys_sc_c, merge(u0_sc, Dict([N3, V, Ω] .=> [N_, V_, Ω_])), (0.0, 2.0))
    prob_ev = ODEProblem(sys_ev_c, merge(u0_ev, Dict([V, Ω] .=> [V_, Ω_])), (0.0, 2.0))
    sol_sc = solve(prob_sc, Tsit5(); reltol = 1.0e-10, abstol = 1.0e-12)
    sol_ev = solve(prob_ev, Tsit5(); reltol = 1.0e-10, abstol = 1.0e-12)
    @test sol_sc.retcode == ReturnCode.Success
    @test sol_ev.retcode == ReturnCode.Success

    # Coherence magnitudes are bounded by 1 in both representations.
    SQA = QuantumCumulants.SecondQuantizedAlgebra
    for op_avg in eqs_sc.states
        op = SQA.undo_average(op_avg)
        mag = [abs(get_solution(sol_sc, op, eqs_sc)(t)) for t in sol_sc.t]
        @test maximum(mag) <= 1.0 + 1.0e-9
    end
    for op_avg in eqs_ev.states
        op = SQA.undo_average(op_avg)
        mag = [abs(get_solution(sol_ev, op, eqs_ev)(t)) for t in sol_ev.t]
        @test maximum(mag) <= 1.0 + 1.0e-9
    end

    # Both state lists are seeded by `s(1,2,i)` and materialise canonical slot
    # 1 to the same rep, so `states[1]` is the same physical coherence ⟨S^{12}⟩
    # in each; the two solves must agree there.
    @test isequal(eqs_sc.states[1], eqs_ev.states[1])
    sc_op = SQA.undo_average(eqs_sc.states[1])
    ev_op = SQA.undo_average(eqs_ev.states[1])
    for τ in (0.25, 0.5, 1.0, 1.5, 2.0)
        v_sc = get_solution(sol_sc, sc_op, eqs_sc)(τ)
        v_ev = get_solution(sol_ev, ev_op, eqs_ev)(τ)
        @test isapprox(v_sc, v_ev; atol = 1.0e-8)
    end
end

@testset "evaluate: literal-key lookup beats alpha-canonical collapse" begin
    # Per-atom states ⟨σ_{j_1}^{12}⟩ and ⟨σ_{j_2}^{12}⟩ are physically distinct
    # after `evaluate(...; limits = N => 2)` populates the index.
    h2 = NLevelSpace(:spin, 2)
    @variables Ω::Real N3::Real
    i = Index(h2, :i, N3, h2)
    s(α, β, idx) = IndexedOperator(Transition(h2, :S, α, β), idx)
    H = Ω * Σ(s(1, 1, i) + s(2, 2, i), i)
    eqs_c = complete(meanfield([s(1, 2, i)], H; order = 1))
    eqs_ev = evaluate(eqs_c; limits = (N3 => 2))

    # Two atoms produce two distinct concrete states.
    @test length(eqs_ev.states) == 2
    @test !isequal(eqs_ev.states[1], eqs_ev.states[2])

    # The states' RHSes reference distinct averaged atoms.
    SQA = QuantumCumulants.SecondQuantizedAlgebra
    op1 = SQA.undo_average(eqs_ev.states[1])
    op2 = SQA.undo_average(eqs_ev.states[2])
    @test !isequal(op1, op2)
end

@testset "evaluate + mtkcompile: NE-stripped fallback path closes" begin
    # A 2-Hilbert-space indexed system (filter + 2-level atom) at order 2.
    # `scale` over the atom space leaves filter-side per-atom NE pairs; the
    # result must still close through `mtkcompile`.
    @variables κ::Real κf::Real g::Real gf::Real R::Real Γ::Real
    @variables Δ::Real ν::Real N::Real M::Real
    δ(k) = IndexedVariable(:δ, k)
    hc = FockSpace(:cavity); hf = FockSpace(:filter)
    ha = NLevelSpace(:atom, 2)
    h = hc ⊗ hf ⊗ ha
    i = Index(h, :i, M, hf)
    j = Index(h, :j, N, ha)
    @qnumbers a::Destroy(h, 1)
    b(k) = IndexedOperator(Destroy(h, :b, 2), k)
    σ(α, β, k) = IndexedOperator(Transition(h, :σ, α, β, 3), k)
    H = Δ * Σ(σ(2, 2, j), j) + Σ(δ(i) * b(i)' * b(i), i) +
        gf * Σ(a' * b(i) + a * b(i)', i) +
        g * Σ(a' * σ(1, 2, j) + a * σ(2, 1, j), j)
    J = [a, b(i), σ(1, 2, j), σ(2, 1, j), σ(2, 2, j)]
    rates = [κ, κf, Γ, R, ν]
    eqs_c = complete(meanfield([a' * a], H, J; rates = rates, order = 2))
    eqs_sc = scale(eqs_c; h = [3])
    eqs_ev = evaluate(eqs_sc; limits = Dict(M => 2))
    @named sys_ev = System(eqs_ev)
    sys_c = mtkcompile(sys_ev)
    @test length(unknowns(sys_c)) == length(eqs_ev.states)
end
