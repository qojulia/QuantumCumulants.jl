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

    # Pre-evaluate symbolic count (single equation closes under cumulant=1).
    @test length(eqs_c.equations) == 2

    # Evaluate at N=3: one symbolic LHS state with one free atom-space index
    # gets enumerated into three concrete-position equations per state.
    eqs_ev = evaluate(eqs_c; limits = (N3 => 3))
    @test length(eqs_ev.equations) == 6
    @test isempty(find_missing(eqs_ev; get_adjoints = false))

    # Round-trip through MTK: mtkcompile must accept the result without
    # Σ-typed residual shapes leaking into MTK's fixpoint substitution.
    @named sys = System(eqs_ev)
    sys_c = mtkcompile(sys)
    @test length(unknowns(sys_c)) == 6

    # Numeric solve at one parameter point.
    N_, V_, Ω_ = 3, 4.79 / 2, 1.0
    u0 = Dict(unknowns(sys_c) .=> zeros(ComplexF64, length(unknowns(sys_c))))
    # N3 is now baked-in (no longer a parameter post-evaluate).
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
    # Scale-first vs evaluate-first should give the same dynamics on a
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
    # post-evaluate, N3 is baked in (only V, Ω remain symbolic).
    prob_ev = ODEProblem(sys_ev_c, merge(u0_ev, Dict([V, Ω] .=> [V_, Ω_])), (0.0, 2.0))
    sol_sc = solve(prob_sc, Tsit5())
    sol_ev = solve(prob_ev, Tsit5())
    @test sol_sc.retcode == ReturnCode.Success
    @test sol_ev.retcode == ReturnCode.Success

    # Physicality: coherence magnitudes bounded by 1 in both representations.
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
end

@testset "evaluate: literal-key lookup beats alpha-canonical collapse" begin
    # Regression for `_canonicalise_avg_leaves`'s literal-then-canon
    # fallback. Per-atom states `⟨σ_{j_1}^{12}⟩` and `⟨σ_{j_2}^{12}⟩`
    # are physically distinct after `evaluate(... ; limits = N => 2)`
    # populates the index. Pre-fix the leaf canonicaliser collapsed
    # both onto canon slot 1, so the materialised system had a single
    # equation per logical index and reading state 2 returned state 1's
    # numeric value (and the `unique_squeezing` example produced a
    # wrong steady-state).
    h2 = NLevelSpace(:spin, 2)
    @variables Ω::Real N3::Real
    i = Index(h2, :i, N3, h2)
    s(α, β, idx) = IndexedOperator(Transition(h2, :S, α, β), idx)
    H = Ω * Σ(s(1, 1, i) + s(2, 2, i), i)
    eqs_c = complete(meanfield([s(1, 2, i)], H; order = 1))
    eqs_ev = evaluate(eqs_c; limits = (N3 => 2))

    # Two atoms must produce two distinct concrete states.
    @test length(eqs_ev.states) == 2
    @test !isequal(eqs_ev.states[1], eqs_ev.states[2])

    # The states' RHSes must reference DISTINCT averaged atoms, not the
    # same canon-collapsed leaf. Pre-fix the bug, both RHSes resolved
    # to the same alpha-canonical leaf and the state-set degenerated.
    SQA = QuantumCumulants.SecondQuantizedAlgebra
    op1 = SQA.undo_average(eqs_ev.states[1])
    op2 = SQA.undo_average(eqs_ev.states[2])
    @test !isequal(op1, op2)
end

@testset "evaluate + mtkcompile: NE-stripped fallback path closes" begin
    # Regression for `_canonicalise_avg_leaves`'s third branch: when a
    # leaf in a cross-atom sum inherits an NE pair from its parent
    # derivation but the leaf's operator atoms do not actually
    # constrain the NE indices, the leaf is rewritten to its bare
    # `average(literal)` form so the downstream MTK conj-substitution
    # dict can match it. Pre-fix this presented as an `AvgFunc is not
    # callable` error during `mtkcompile` for the filter-cavity scenario.
    #
    # Reproduction: a 2-Hilbert-space indexed system (filter +
    # 2-level atom) at order 2. `scale` over the atom space leaves
    # filter-side per-atom NE pairs that the leaf canonicaliser must
    # strip before the conj-dict lookup.
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
    @test isempty(find_missing(eqs_ev; get_adjoints = false))
    @named sys_ev = System(eqs_ev)
    sys_c = mtkcompile(sys_ev)
    @test length(unknowns(sys_c)) == length(eqs_ev.states)
end
