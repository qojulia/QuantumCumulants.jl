using QuantumCumulants
using Symbolics: Symbolics, @variables
using SymbolicUtils: SymbolicUtils
using ModelingToolkitBase: ModelingToolkitBase, @named, mtkcompile, ODEProblem, unknowns, t_nounits as t
using OrdinaryDiffEqTsit5: Tsit5, solve, ReturnCode
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
    @test isempty(find_missing(eqs_ev))

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

@testset "evaluate: #277 sum over identical=false double-indexed variable" begin
    # H = Σ_{i,j} Ω_ij σ11^i σ22^j with Ω_ii = 0. evaluate must unroll the sum
    # into concrete Ω[i,j] entries, leaving no symbolic index behind, and must
    # honour identical=false so no diagonal Ω[k,k] survives.
    h = NLevelSpace(:atom, 2)
    @variables N::Real
    Ω(a, b) = DoubleIndexedVariable(:Ω, a, b; identical = false)
    σ(x, y, idx) = IndexedOperator(Transition(h, :σ, x, y), idx)
    i = Index(h, :i, N, h); j = Index(h, :j, N, h); k = Index(h, :k, N, h)
    H = Σ(Ω(i, j) * σ(1, 1, i) * σ(2, 2, j), j, i)
    eqs_c = complete(meanfield([σ(1, 2, k)], H; order = 1))

    eqs_ev = evaluate(eqs_c; limits = (N => 2))
    @test length(eqs_ev.equations) == 4

    # A materialised index slot is a `getindex` arg: a numeric constant once
    # resolved, but wrapped as `BasicSymbolic` (so `arg isa Integer` is false and
    # `Int(arg)` throws). `Symbolics.value` collapses the wrap to the integer.
    slots = Set{Tuple{Int, Int}}()
    for eq in eqs_ev.equations, v in Symbolics.get_variables(eq.rhs)
        u = SymbolicUtils.unwrap(v)
        if SymbolicUtils.iscall(u) && SymbolicUtils.operation(u) === getindex
            a = SymbolicUtils.arguments(u)
            occursin("Ω", string(a[1])) &&
                push!(slots, (Int(Symbolics.value(a[2])), Int(Symbolics.value(a[3]))))
        end
    end
    # identical=false ⇒ no diagonal Ω[k,k] survives; only off-diagonal entries do.
    @test slots == Set([(1, 2), (2, 1)])

    # The populations ⟨σ22^k⟩ are conserved: a diagonal-free Ω makes H commute
    # with each σ22^k, so exactly the two ⟨σ22⟩ equations vanish. A zero rhs is a
    # `BasicSymbolic` that is not `isequal` to `0`, but carries no free variables.
    @test count(eq -> isempty(Symbolics.get_variables(eq.rhs)), eqs_ev.equations) == 2
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

@testset "evaluate: #288/#198 double-indexed coupling leaves no free sum index" begin
    # Regression for #198/#288: an order-2 term inside a sum whose coefficient is a
    # double-indexed variable J(l,k) factorises into a product of averages. The bound
    # index `k` lives in the c-number coefficient AND in several averaged leaves at
    # once, so stamping the scope onto a single leaf left `k` dangling on the
    # coefficient (`UndefVarError: δjk`). The sum node must wrap the whole
    # index-dependent product, so after `evaluate` no RHS carries a free index.
    ha = NLevelSpace(:atoms, 2); hb = FockSpace(:bus); h = hb ⊗ ha
    b = Destroy(h, :b)
    @variables L::Real Δ_b::Real κ::Real γ::Real ν::Real
    j = Index(h, :j, L, ha); k = Index(h, :k, L, ha); l = Index(h, :l, L, ha)
    χ(x) = IndexedVariable(:χ, x)
    Δ(x) = IndexedVariable(:Δ, x)
    J(x, y) = DoubleIndexedVariable(:J, x, y)
    σz(i) = IndexedOperator(Transition(h, :σ, 2, 2), i)
    σp(i) = IndexedOperator(Transition(h, :σ, 1, 2), i)
    σm(i) = IndexedOperator(Transition(h, :σ, 2, 1), i)
    H = Δ_b * b' * b + Σ(Δ(j) * σz(j), j) +
        Σ(χ(j) * (b' * σm(j) + b * σp(j)), j) +
        Σ(J(j, k) * (σp(j) * σm(k)), j, k)
    ops = [b' * b, σz(l)]
    Jops = [b, σp(j), σm(j)]
    eqs_c = complete(meanfield(ops, H, Jops; rates = [κ, γ, ν], order = 2))
    eqs_ev = evaluate(eqs_c; limits = (L => 2))

    # The system closes: every moment on an RHS has its own equation.
    @test isempty(find_missing(eqs_ev))

    # The summation indices `j`, `k` are fully unrolled; only concrete atom labels
    # (the observable index `l` materialised to `l_1`/`l_2`) may remain. A surviving
    # `j`/`k` is the #198 dangling-index bug.
    SQA = QuantumCumulants.SecondQuantizedAlgebra
    for eq in eqs_ev.equations
        names = Set(SQA.index_name(i) for i in SQA.get_indices(SymbolicUtils.unwrap(eq.rhs)))
        @test :j ∉ names && :k ∉ names
    end

    # The concrete system builds and solves (the #198 acceptance criterion).
    sys_c = mtkcompile(System(eqs_ev; name = :sys288))
    init = initial_values(eqs_ev, zeros(ComplexF64, length(eqs_ev.equations)))
    pmap = Dict{Any, Any}()
    for p in ModelingToolkitBase.parameters(sys_c)
        nm = string(p)
        nm == "Δ_b" && (pmap[p] = 0.0)
        nm == "κ" && (pmap[p] = 0.5)
        nm == "γ" && (pmap[p] = 0.1)
        nm == "ν" && (pmap[p] = 0.05)
        nm == "χ" && (pmap[p] = [0.7, 0.9])
        nm == "Δ" && (pmap[p] = [-2.5, -1.9])
        nm == "J" && (pmap[p] = [0.0 0.1; 0.12 0.0])
    end
    prob = ODEProblem(sys_c, merge(init, pmap), (0.0, 5.0))
    sol = solve(prob, Tsit5())
    @test sol.retcode == ReturnCode.Success
end

@testset "evaluate: coefficient-only summation index matches single-atom exact (#198)" begin
    # The numerical half of #198. Each atom `l` is driven at Rabi frequency
    # Ω_l = Σ_k u(l,k) built from a DoubleIndexedVariable, so the summation index `k`
    # appears ONLY in the coefficient (no operator carries it). The atoms are otherwise
    # uncoupled, so a second-order cumulant expansion is exact and every atom must
    # reproduce a single two-level atom driven at its own Ω_l. Before the fix the
    # diagonal `k=l` slice of the coefficient sum was over-counted (the off-diagonal body
    # picked up `u(l,k)+u(l,l)` and the diagonal was applied a second time), which
    # inflated the populations (atom 1 read 0.447 against the exact 0.325) while the system
    # still built and solved, so the existing structural #288/#198 test could not catch it.
    # Needs SecondQuantizedAlgebra ≥ 0.9.3 for the companion diagonal-split fix.
    ha = NLevelSpace(:atom, 2); h = ha
    @variables L::Real Δ::Real γ::Real
    u(a, b) = DoubleIndexedVariable(:u, a, b)
    j = Index(h, :j, L, ha); k = Index(h, :k, L, ha); l = Index(h, :l, L, ha)
    σ(a, b, idx) = IndexedOperator(Transition(h, :σ, a, b), idx)
    H = Δ * Σ(σ(2, 2, j), j) + Σ(Σ(u(j, k) * (σ(2, 1, j) + σ(1, 2, j)), j), k)
    eqs_ev = evaluate(
        complete(meanfield(σ(2, 2, l), H, [σ(1, 2, j)]; rates = [γ], order = 2));
        limits = (L => 2),
    )
    @test isempty(find_missing(eqs_ev))

    umat = [0.3 0.5; 0.4 0.2]; Δv, γv = 0.7, 1.0
    sys = mtkcompile(System(eqs_ev; name = :cc198))
    init = initial_values(eqs_ev, zeros(ComplexF64, length(eqs_ev.equations)))
    pmap = Dict{Any, Any}()
    for p in ModelingToolkitBase.parameters(sys)
        nm = string(p)
        nm == "Δ" && (pmap[p] = Δv)
        nm == "γ" && (pmap[p] = γv)
        nm == "u" && (pmap[p] = umat)
    end
    sol = solve(ODEProblem(sys, merge(init, pmap), (0.0, 5.0)), Tsit5(); reltol = 1.0e-10, abstol = 1.0e-12)
    @test sol.retcode == ReturnCode.Success

    # Reference: a single two-level atom driven at a scalar Ω (order 2 is exact for one atom).
    haR = NLevelSpace(:atomR, 2)
    @variables Ω::Real
    σR(a, b) = Transition(haR, :σ, a, b)
    HR = Δ * σR(2, 2) + Ω * (σR(2, 1) + σR(1, 2))
    eqsR = complete(meanfield(σR(2, 2), HR, [σR(1, 2)]; rates = [γ], order = 2))
    sysR = mtkcompile(System(eqsR; name = :ref198))
    function ref_solve(Ωv)
        initR = initial_values(eqsR, zeros(ComplexF64, length(eqsR.equations)))
        pmapR = Dict{Any, Any}()
        for p in ModelingToolkitBase.parameters(sysR)
            nm = string(p)
            nm == "Δ" && (pmapR[p] = Δv)
            nm == "γ" && (pmapR[p] = γv)
            nm == "Ω" && (pmapR[p] = Ωv)
        end
        solve(ODEProblem(sysR, merge(initR, pmapR), (0.0, 5.0)), Tsit5(); reltol = 1.0e-10, abstol = 1.0e-12)
    end
    sol1 = ref_solve(sum(umat[1, :]))   # Ω_1 = u(1,1) + u(1,2)
    sol2 = ref_solve(sum(umat[2, :]))   # Ω_2 = u(2,1) + u(2,2)
    @test sol1.retcode == ReturnCode.Success
    @test sol2.retcode == ReturnCode.Success

    # Each atom's excited population tracks the single-atom solution at its own Ω_l.
    for τ in (0.5, 1.0, 2.0, 5.0)
        @test isapprox(
            real(get_solution(sol, σ(2, 2, l(1)), eqs_ev)(τ)),
            real(get_solution(sol1, σR(2, 2), eqsR)(τ)); atol = 1.0e-6,
        )
        @test isapprox(
            real(get_solution(sol, σ(2, 2, l(2)), eqs_ev)(τ)),
            real(get_solution(sol2, σR(2, 2), eqsR)(τ)); atol = 1.0e-6,
        )
    end
end

@testset "evaluate: #288 time-dependent multi-factor indexed coefficient keeps scope" begin
    # The issue's exact failure mode: a sum coefficient `sin(δ(k)·t)·u(k)` is a *product* of a
    # time+index-dependent factor and an index-dependent factor. Canonicalising such a
    # multi-factor coefficient is what stripped the old leaf-level scope; the sum node must keep
    # the coefficient scoped so `evaluate` leaves no free summation index and preserves `t`.
    ha = NLevelSpace(:atoms, 2); hb = FockSpace(:bus); h = hb ⊗ ha
    b = Destroy(h, :b)
    @variables L::Real κ::Real γ::Real
    k = Index(h, :k, L, ha); l = Index(h, :l, L, ha)
    u(x) = IndexedVariable(:u, x)
    δ(x) = IndexedVariable(:δ, x)
    σz(idx) = IndexedOperator(Transition(h, :σ, 2, 2), idx)
    σp(idx) = IndexedOperator(Transition(h, :σ, 1, 2), idx)
    σm(idx) = IndexedOperator(Transition(h, :σ, 2, 1), idx)
    H = Σ(sin(δ(k) * t) * u(k) * (b' * σm(k) + b * σp(k)), k)
    eqs_c = complete(meanfield([b' * b, σz(l)], H, [b, σm(k)]; rates = [κ, γ], order = 2))
    eqs_ev = evaluate(eqs_c; limits = (L => 2))

    @test isempty(find_missing(eqs_ev))
    SQA = QuantumCumulants.SecondQuantizedAlgebra
    for eq in eqs_ev.equations
        nms = Set(SQA.index_name(x) for x in SQA.get_indices(SymbolicUtils.unwrap(eq.rhs)))
        @test :k ∉ nms
    end
    # The independent variable survives the unrolling (time-dependence preserved).
    t_uw = SymbolicUtils.unwrap(t)
    @test any(
        eq -> any(v -> isequal(v, t_uw), Symbolics.get_variables(eq.rhs)),
        eqs_ev.equations,
    )
end

@testset "evaluate: invalid limits type throws" begin
    h2 = NLevelSpace(:spin, 2)
    @variables Ω::Real N3::Real
    i = Index(h2, :i, N3, h2)
    s(α, β, idx) = IndexedOperator(Transition(h2, :S, α, β), idx)
    H = Ω * Σ(s(1, 1, i) + s(2, 2, i), i)
    eqs_c = complete(meanfield([s(1, 2, i)], H; order = 1))
    # `limits` must be a Pair, Tuple of Pairs, or Dict; a bare scalar is rejected.
    @test_throws ArgumentError evaluate(eqs_c; limits = 5)
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
