using QCNew
using Symbolics: Symbolics, @variables
using ModelingToolkitBase: @named, mtkcompile, ODEProblem, unknowns
using OrdinaryDiffEq: Tsit5, solve, ReturnCode
using Test

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

    eqs_sc = scale(eqs_c)
    @test length(eqs_sc.equations) >= 1
    @test isempty(find_missing(eqs_sc; get_adjoints = false))

    @named sys = System(eqs_sc)
    sys_c = mtkcompile(sys)
    @test length(unknowns(sys_c)) >= 1
end

@testset "indexed_scale: per-Hilbert-space scale order independence" begin
    # Scaling two indexed subspaces in either order must produce systems of equal size.
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

    @named sys_sc = System(eqs_sc); sys_sc_c = mtkcompile(sys_sc)
    @named sys_ev = System(eqs_ev); sys_ev_c = mtkcompile(sys_ev)

    N_, V_, Ω_ = 3, 4.79 / 2, 1.0
    u0_sc = Dict(unknowns(sys_sc_c) .=> zeros(ComplexF64, length(unknowns(sys_sc_c))))
    u0_ev = Dict(unknowns(sys_ev_c) .=> zeros(ComplexF64, length(unknowns(sys_ev_c))))

    prob_sc = ODEProblem(sys_sc_c, merge(u0_sc, Dict([N3, V, Ω] .=> [N_, V_, Ω_])), (0.0, 2.0))
    prob_ev = ODEProblem(sys_ev_c, merge(u0_ev, Dict([V, Ω] .=> [V_, Ω_])), (0.0, 2.0))
    sol_sc = solve(prob_sc, Tsit5(); abstol = 1.0e-10, reltol = 1.0e-10)
    sol_ev = solve(prob_ev, Tsit5(); abstol = 1.0e-10, reltol = 1.0e-10)
    @test sol_sc.retcode == ReturnCode.Success
    @test sol_ev.retcode == ReturnCode.Success

    # The scale and evaluate paths must produce the same physics. Compare ⟨s_12⟩
    # at the same time points.
    SQA = QCNew.SecondQuantizedAlgebra
    obs_op = SQA.undo_average(eqs_sc.states[1])
    val_sc_end = get_solution(sol_sc, obs_op, eqs_sc)(sol_sc.t[end])
    obs_op_ev_1 = SQA.undo_average(eqs_ev.states[1])
    obs_op_ev_2 = SQA.undo_average(eqs_ev.states[2])
    obs_op_ev_3 = SQA.undo_average(eqs_ev.states[3])
    val_ev_1 = get_solution(sol_ev, obs_op_ev_1, eqs_ev)(sol_ev.t[end])
    val_ev_2 = get_solution(sol_ev, obs_op_ev_2, eqs_ev)(sol_ev.t[end])
    val_ev_3 = get_solution(sol_ev, obs_op_ev_3, eqs_ev)(sol_ev.t[end])
    @test isapprox(val_ev_1, val_sc_end; atol = 1.0e-6)
    @test isapprox(val_ev_2, val_sc_end; atol = 1.0e-6)
    @test isapprox(val_ev_3, val_sc_end; atol = 1.0e-6)
    # Permutation symmetry: the three atoms give the same value.
    @test isapprox(val_ev_1, val_ev_2; atol = 1.0e-8)
    @test isapprox(val_ev_2, val_ev_3; atol = 1.0e-8)
end

@testset "indexed vs explicit: Tavis-Cummings physics agreement" begin
    # Build the same 2-atom driven Tavis-Cummings system explicitly (two atom
    # Hilbert spaces, no `Index`) and indexed (one atom space + `Index`, reduced
    # via `scale` and via `evaluate(N => 2)`), then require the cavity observables
    # ⟨a'a⟩ and ⟨a⟩ to agree across all three.
    @variables Δ::Real g::Real κ::Real γ::Real η::Real N::Real

    # (a) explicit two-atom build
    hc = FockSpace(:cavity); ha1 = NLevelSpace(:atom1, 2); ha2 = NLevelSpace(:atom2, 2)
    hx = hc ⊗ ha1 ⊗ ha2
    ax = Destroy(hx, :a, 1)
    sx1(α, β) = Transition(hx, :σ1, α, β, 2)
    sx2(α, β) = Transition(hx, :σ2, α, β, 3)
    Hx = -Δ * ax' * ax + η * (ax' + ax) +
        g * (ax' * sx1(1, 2) + ax * sx1(2, 1)) +
        g * (ax' * sx2(1, 2) + ax * sx2(2, 1))
    eqs_ex = complete(
        meanfield(
            [ax' * ax, sx1(2, 2), sx2(2, 2)], Hx,
            [ax, sx1(1, 2), sx2(1, 2)]; rates = [κ, γ, γ], order = 2
        )
    )

    # (b) indexed build
    hci = FockSpace(:cavity); hai = NLevelSpace(:atom, 2); hi = hci ⊗ hai
    ii = Index(hi, :i, N, hai)
    ai = Destroy(hi, :a, 1)
    σi(α, β, k) = IndexedOperator(Transition(hi, :σ, α, β, 2), k)
    Hi = -Δ * ai' * ai + η * (ai' + ai) +
        g * Σ(ai' * σi(1, 2, ii) + ai * σi(2, 1, ii), ii)
    eqs_i = complete(
        meanfield(
            [ai' * ai, σi(2, 2, ii)], Hi,
            [ai, σi(1, 2, ii)]; rates = [κ, γ], order = 2
        )
    )
    eqs_sc = scale(eqs_i)
    eqs_ev = evaluate(eqs_i; limits = (N => 2))

    function _solve(eqs, pvals)
        sys = mtkcompile(System(eqs; name = :s))
        u0 = Dict(unknowns(sys) .=> zeros(ComplexF64, length(unknowns(sys))))
        prob = ODEProblem(sys, merge(u0, pvals), (0.0, 3.0))
        return solve(prob, Tsit5(); reltol = 1.0e-10, abstol = 1.0e-12)
    end

    base = Dict(Δ => 0.0, g => 1.5, κ => 0.7, γ => 0.3, η => 1.2)
    sol_ex = _solve(eqs_ex, base)
    sol_ev = _solve(eqs_ev, base)
    sol_sc = _solve(eqs_sc, merge(base, Dict(N => 2.0)))
    @test sol_ex.retcode == ReturnCode.Success
    @test sol_ev.retcode == ReturnCode.Success
    @test sol_sc.retcode == ReturnCode.Success

    # Cavity observables must agree across all three representations.
    for τ in (0.5, 1.0, 2.0, 3.0)
        n_ex = get_solution(sol_ex, ax' * ax, eqs_ex)(τ)
        a_ex = get_solution(sol_ex, ax, eqs_ex)(τ)
        @test isapprox(get_solution(sol_ev, ai' * ai, eqs_ev)(τ), n_ex; atol = 1.0e-8)
        @test isapprox(get_solution(sol_sc, ai' * ai, eqs_sc)(τ), n_ex; atol = 1.0e-8)
        @test isapprox(get_solution(sol_ev, ai, eqs_ev)(τ), a_ex; atol = 1.0e-8)
        @test isapprox(get_solution(sol_sc, ai, eqs_sc)(τ), a_ex; atol = 1.0e-8)
    end
end
