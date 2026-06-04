using QCNew
using ModelingToolkitBase
using ModelingToolkitBase: t_nounits as t
using Symbolics: @variables, Num
using OrdinaryDiffEqTsit5: ODEProblem, solve, Tsit5
using Test

# Regression tests for every example in `examples/`.

import QCNew.SecondQuantizedAlgebra as SQA

_phase_inv(x) = 0
_phase_inv(::Destroy) = -1
_phase_inv(::Create) = 1
function _phase_inv(τ::Transition)
    return (τ.i == 2 && τ.j == 1) ? 1 :
        (τ.i == 1 && τ.j == 2) ? -1 : 0
end
function _phase_inv(q::SQA.QAdd)
    for (term, _) in q.arguments
        p = 0
        for op in term.ops; p += _phase_inv(op); end
        return p
    end
    return 0
end
function _phase_inv(avg)
    SQA.is_average(avg) || return 0
    return _phase_inv(SQA.undo_average(avg))
end
phase_invariant(x) = iszero(_phase_inv(x))

@testset "mollow" begin
    @variables Δ Ω γ
    h = NLevelSpace(:atom, (:g, :e))
    σ(α, β) = Transition(h, :σ, α, β)
    H = Δ * σ(:e, :e) + Ω * (σ(:e, :g) + σ(:g, :e))
    eqs = meanfield([σ(:e, :g), σ(:e, :e)], H, [σ(:g, :e)]; rates = [γ])
    complete!(eqs)
    @test length(eqs.equations) == 2
    sys_c = mtkcompile(System(eqs; name = :mollow_sys))
    @test length(unknowns(sys_c)) == 2
    u0 = zeros(ComplexF64, length(eqs.equations))
    init = initial_values(eqs, u0)
    ps = Dict(Δ => 1.0, Ω => 0.5, γ => 0.1)
    prob = ODEProblem(sys_c, merge(init, ps), (0.0, 80.0))
    sol = solve(prob, Tsit5(); reltol = 1.0e-10, abstol = 1.0e-10)
    @test isapprox(
        real(get_solution(sol, σ(:e, :e), eqs).(sol.t[end])),
        0.16491193049237177; rtol = 1.0e-5,
    )
    @test isapprox(
        get_solution(sol, σ(:e, :g), eqs).(sol.t[end]),
        -0.33088938061453155 + 0.016559911920571824im; rtol = 1.0e-5,
    )
    c = CorrelationFunction(σ(:e, :g), σ(:g, :e), eqs; steady_state = true)
    @test length(c.eqs.equations) == 3
    ps_tup = (Δ, Ω, γ)
    p0 = (0.0, 2.0, 1.0)
    sys = mtkcompile(System(eqs; name = :mollow_steady))
    u0_ss = zeros(ComplexF64, length(unknowns(sys)))
    dict_ss = merge(Dict(unknowns(sys) .=> u0_ss), Dict(ps_tup .=> p0))
    prob_ss = ODEProblem(sys, dict_ss, (0.0, 20.0))
    sol_ss = solve(prob_ss, Tsit5(); reltol = 1.0e-10, abstol = 1.0e-10)
    S = Spectrum(c, ps_tup)
    s_vals = S([0.0, 2.0, -2.0], sol_ss.u[end], p0)
    @test isapprox(s_vals[1], 1.0257952485186328; rtol = 1.0e-5)
    @test isapprox(s_vals[2], 0.14359486828829415; rtol = 1.0e-5)
    # spectrum symmetry S(ω) = S(-ω)
    @test isapprox(s_vals[2], s_vals[3]; rtol = 1.0e-10)
end

@testset "single-atom-laser-spectrum" begin
    @variables Δ g γ κ ν
    hf = FockSpace(:cavity); ha = NLevelSpace(:atom, (:g, :e))
    h = hf ⊗ ha
    a = Destroy(h, :a); s = Transition(h, :σ, :g, :e)
    H = Δ * a' * a + g * (a' * s + a * s')
    eq_n = meanfield(a' * a, H, [a, s, s']; rates = [κ, γ, ν], order = 2)
    eqs = complete(eq_n; filter_func = phase_invariant)
    @test length(eqs.equations) == 4
    sys_c = mtkcompile(System(eqs; name = :sal_sys))
    @test length(unknowns(sys_c)) == 4
    u0 = zeros(ComplexF64, length(eqs.equations))
    init = initial_values(eqs, u0)
    ps = Dict(Δ => 0.0, g => 1.5, γ => 0.25, κ => 1.0, ν => 4.0)
    prob = ODEProblem(sys_c, merge(init, ps), (0.0, 20.0))
    sol = solve(prob, Tsit5(); reltol = 1.0e-10, abstol = 1.0e-10)
    @test isapprox(
        real(get_solution(sol, a' * a, eqs).(sol.t[end])),
        1.4835417514686384; rtol = 1.0e-5,
    )
    c = CorrelationFunction(
        a', a, eqs; steady_state = true, filter_func = phase_invariant,
    )
    @test length(c.eqs.equations) == 2
    ps_tup = (Δ, g, γ, κ, ν)
    p0 = (0.0, 1.5, 0.25, 1.0, 4.0)
    sys_steady = mtkcompile(System(eqs; name = :sal_steady))
    u0_ss = zeros(ComplexF64, length(unknowns(sys_steady)))
    dict_ss = merge(Dict(unknowns(sys_steady) .=> u0_ss), Dict(ps_tup .=> p0))
    prob_ss = ODEProblem(sys_steady, dict_ss, (0.0, 20.0))
    sol_ss = solve(prob_ss, Tsit5(); reltol = 1.0e-10, abstol = 1.0e-10)
    S = Spectrum(c, ps_tup)
    s_vals = S([-0.5, 0.0, 0.5], sol_ss.u[end], p0)
    @test isapprox(s_vals[2], 12.01916023864571; rtol = 1.0e-4)
    @test isapprox(s_vals[1], 2.683110742979243; rtol = 1.0e-4)
    # spectrum symmetry S(ω) = S(-ω)
    @test isapprox(s_vals[1], s_vals[3]; rtol = 1.0e-10)
end

@testset "ramsey_spectroscopy" begin
    @variables Δ Ω Γ ν
    @register_symbolic _ramsey_f(tt)
    h = NLevelSpace(:atom, 2)
    σ(i, j) = Transition(h, :σ, i, j)
    eqs_seed = meanfield(
        [σ(2, 2), σ(1, 2)], -Δ * σ(2, 2),
        [σ(1, 2), σ(2, 2)]; rates = [Γ, ν],
    )
    tv = eqs_seed.iv
    H = -Δ * σ(2, 2) + _ramsey_f(tv) * Ω / 2 * (σ(1, 2) + σ(2, 1))
    eqs = meanfield(
        [σ(2, 2), σ(1, 2)], H, [σ(1, 2), σ(2, 2)];
        rates = [Γ, ν], iv = tv,
    )
    complete!(eqs)
    @test length(eqs.equations) == 2
    sys_c = mtkcompile(System(eqs; name = :ramsey_sys))
    @test length(unknowns(sys_c)) == 2
    _ramsey_f(_) = 1.0
    u0 = zeros(ComplexF64, length(eqs.equations))
    init = initial_values(eqs, u0)
    ps = Dict(Δ => 0.0, Ω => 1.0, Γ => 0.1, ν => 0.05)
    prob = ODEProblem(sys_c, merge(init, ps), (0.0, 20.0))
    sol = solve(prob, Tsit5(); reltol = 1.0e-10, abstol = 1.0e-10)
    @test isapprox(
        real(get_solution(sol, σ(2, 2), eqs).(sol.t[end])),
        0.45407705143934646; rtol = 1.0e-5,
    )
    @test isapprox(
        get_solution(sol, σ(1, 2), eqs).(sol.t[end]),
        0.0 - 0.12468144623895208im; rtol = 1.0e-5,
    )
end

@testset "many-atom-laser" begin
    N = 2
    @variables κ g Γ23 Γ13 Γ12 Ω Δc Δ3
    hf = FockSpace(:cavity)
    ha = ⊗([NLevelSpace(Symbol(:atom, i), 3) for i in 1:N]...)
    h = hf ⊗ ha
    a = Destroy(h, :a)
    σ(i, j, k) = Transition(h, Symbol("σ_{$k}"), i, j, k + 1)
    H = -Δc * a' * a +
        sum(g * (a' * σ(1, 2, i) + a * σ(2, 1, i)) for i in 1:N) +
        sum(Ω * (σ(3, 1, i) + σ(1, 3, i)) for i in 1:N) -
        sum(Δ3 * σ(3, 3, i) for i in 1:N)
    J = [a; [σ(1, 2, i) for i in 1:N]; [σ(1, 3, i) for i in 1:N]; [σ(2, 3, i) for i in 1:N]]
    rates = [κ; [Γ12 for _ in 1:N]; [Γ13 for _ in 1:N]; [Γ23 for _ in 1:N]]
    ops = [a' * a, σ(2, 2, 1), σ(3, 3, 1)]
    eqs = meanfield(ops, H, J; rates = rates, order = 2)
    complete!(eqs)
    @test length(eqs.equations) == 117
    sys_c = mtkcompile(System(eqs; name = :mal_sys))
    @test length(unknowns(sys_c)) == 117
    u0 = zeros(ComplexF64, length(eqs.equations))
    init = initial_values(eqs, u0)
    ps = Dict(
        κ => 0.5, g => 2.0, Γ23 => 20.0, Γ13 => 2.0, Γ12 => 1.0,
        Ω => 10.0, Δc => 0.0, Δ3 => 0.0,
    )
    prob = ODEProblem(sys_c, merge(init, ps), (0.0, 10.0))
    sol = solve(prob, Tsit5(); reltol = 1.0e-9, abstol = 1.0e-9)
    @test isapprox(
        real(get_solution(sol, a' * a, eqs).(sol.t[end])),
        12.908870121878037; rtol = 1.0e-4,
    )
    @test isapprox(
        real(get_solution(sol, σ(2, 2, 1), eqs).(sol.t[end])),
        0.4479907519944381; rtol = 1.0e-4,
    )
end

@testset "optomechanical-cooling" begin
    hc = FockSpace(:cavity); hm = FockSpace(:motion); h = hc ⊗ hm
    @qnumbers a::Destroy(h, 1) b::Destroy(h, 2)
    @variables Δ ωm E G κ
    H = -Δ * a' * a + ωm * b' * b + G * a' * a * (b + b') + E * (a + a')
    eqs = meanfield([a' * a, b' * b], H, [a]; rates = [κ], order = 2)
    eqs_c = complete!(deepcopy(eqs))
    @test length(eqs_c.equations) == 14
    sys_c = mtkcompile(System(eqs_c; name = :om_sys))
    @test length(unknowns(sys_c)) == 14
    u0 = zeros(ComplexF64, length(eqs_c.equations))
    init = initial_values(
        eqs_c, u0; defaults = Dict(average(b' * b) => 100.0 + 0im),
    )
    ps = Dict(Δ => -1.0, ωm => 1.0, E => 0.5, G => 0.05, κ => 0.1)
    prob = ODEProblem(sys_c, merge(init, ps), (0.0, 50.0))
    sol = solve(prob, Tsit5(); reltol = 1.0e-9, abstol = 1.0e-9)
    @test isapprox(
        real(get_solution(sol, a' * a, eqs_c).(sol.t[end])),
        0.28986377256026735; rtol = 1.0e-4,
    )
    @test isapprox(
        real(get_solution(sol, b' * b, eqs_c).(sol.t[end])),
        0.03417729886018334; rtol = 1.0e-3,
    )
end

@testset "waveguide" begin
    M_p = 2; M_np = 2; M = M_p + M_np
    h_spin(i) = SpinSpace(Symbol("spin_$(i)"))
    h = tensor([h_spin(i) for i in 1:M]...)
    Sx(i) = Spin(h, Symbol("S$(i)"), 1, i)
    Sy(i) = Spin(h, Symbol("S$(i)"), 2, i)
    Sz(i) = Spin(h, Symbol("S$(i)"), 3, i)
    Sm(i) = (Sx(i) - 1im * Sy(i))
    Sp(i) = (Sx(i) + 1im * Sy(i))
    _Ω_cache = Dict()
    _Γ_cache = Dict()
    Ωp(i, j) = (k = i > j ? (j, i) : (i, j);
        get!(_Ω_cache, k, (ModelingToolkitBase.@variables $(Symbol("Ω_$(k[1])_$(k[2])")))[1]))
    Γp(i, j) = (k = i > j ? (j, i) : (i, j);
        get!(_Γ_cache, k, (ModelingToolkitBase.@variables $(Symbol("Γ_$(k[1])_$(k[2])")))[1]))
    H = sum((i ≠ j) * Ωp(i, j) * Sp(i) * Sm(j) for i in 1:M for j in 1:M)
    J = [Sm(c1) for c1 in 1:M]
    rates = [Γp(c1, c2) for c1 in 1:M, c2 in 1:M]
    S(i) = [Sx(i), Sy(i), Sz(i)]
    SiSi(i) = [Sx(i)Sx(i), Sx(i)Sy(i), Sx(i)Sz(i), Sy(i)Sy(i), Sy(i)Sz(i), Sz(i)Sz(i)]
    ops = []
    for i in 1:M; push!(ops, S(i)...); end
    for i in 1:M; push!(ops, SiSi(i)...); end
    for i in 1:M, j in i:M
        if i ≠ j
            for α in 1:3, β in 1:3
                push!(ops, S(i)[α] * S(j)[β])
            end
        end
    end
    @test length(ops) == 90
    eqs = meanfield(ops, H, J; rates = rates, order = 2)
    @test length(eqs.equations) == 90
    sys_c = mtkcompile(System(eqs; name = :wg_sys))
    @test length(unknowns(sys_c)) == 90
    u0 = zeros(ComplexF64, length(eqs.equations))
    for k in 1:M
        for (idx, s) in enumerate(eqs.states)
            if string(s) == string(average(Sz(k)))
                u0[idx] = 0.5 + 0im
            elseif string(s) == string(average(Sz(k) * Sz(k)))
                u0[idx] = 0.25 + 0im
            end
        end
    end
    init = initial_values(eqs, u0)
    pmap_dict = Dict{Any, Any}()
    for i in 1:M, j in (i + 1):M
        pmap_dict[Ωp(i, j)] = 0.1 * (i + j)
    end
    for i in 1:M, j in 1:M
        pmap_dict[Γp(i, j)] = (i == j ? 1.0 : 0.2)
    end
    prob = ODEProblem(sys_c, merge(init, pmap_dict), (0.0, 2.0))
    sol = solve(prob, Tsit5(); reltol = 1.0e-9, abstol = 1.0e-9)
    @test isapprox(
        real(get_solution(sol, Sz(1), eqs).(sol.t[end])),
        -0.029879167913446643; rtol = 1.0e-4,
    )
    @test isapprox(
        real(get_solution(sol, Sz(M), eqs).(sol.t[end])),
        -0.030352188453010635; rtol = 1.0e-4,
    )
end

@testset "excitation-transport-chain (N=10 structural)" begin
    N = 10
    h_chain = ⊗([NLevelSpace(Symbol(:atom, i), (:g, :e)) for i in 1:N]...)
    σ_ch(i, j, k) = Transition(h_chain, Symbol(:σ_, k), i, j, k)
    @variables Ω γ Δ J0
    x_ = [first(@variables $(Symbol("x_$i"))) for i in 1:N]
    Jc(xi, xj) = J0 / abs(xi - xj)^3
    H = -Δ * sum(σ_ch(:e, :e, k) for k in 1:N) +
        Ω * (σ_ch(:e, :g, 1) + σ_ch(:g, :e, 1)) +
        sum(
        Jc(x_[k], x_[k + 1]) * (
            σ_ch(:e, :g, k) * σ_ch(:g, :e, k + 1) +
                σ_ch(:g, :e, k) * σ_ch(:e, :g, k + 1)
        )
            for k in 1:(N - 1)
    )
    c_ops = [σ_ch(:g, :e, k) for k in 1:N]
    eqs = meanfield(σ_ch(:g, :e, 1), H, c_ops; rates = [γ for _ in 1:N], order = 2)
    complete!(eqs)
    # The full order-2 closure has 434 distinct moments. complete!'s iteration cap
    # must be high enough to reach the fixpoint (the N=10 chain needs > 200 BFS
    # node-expansions; the cap was raised in the moment-class rebuild so closure
    # finishes instead of truncating). Coordinate-consistent find_missing (Task 2)
    # confirms the closed system has no missing leaves.
    @test isempty(find_missing(eqs; get_adjoints = false))
    @test length(eqs.equations) == 434
    sys_c = mtkcompile(System(eqs; name = :chain_sys))
    @test length(unknowns(sys_c)) == 434
end

# N=4 chain integrated as a mean-field ODE, pinning the end-of-chain population.
@testset "excitation-transport-chain (N=4 numerical)" begin
    N = 4
    h_chain = ⊗([NLevelSpace(Symbol(:atom, i), (:g, :e)) for i in 1:N]...)
    σ_ch(i, j, k) = Transition(h_chain, Symbol(:σ_, k), i, j, k)
    @variables Ω γ Δ J0
    x_ = [first(@variables $(Symbol("x4_$i"))) for i in 1:N]
    Jc(xi, xj) = J0 / abs(xi - xj)^3
    H = -Δ * sum(σ_ch(:e, :e, k) for k in 1:N) +
        Ω * (σ_ch(:e, :g, 1) + σ_ch(:g, :e, 1)) +
        sum(
        Jc(x_[k], x_[k + 1]) * (
            σ_ch(:e, :g, k) * σ_ch(:g, :e, k + 1) +
                σ_ch(:g, :e, k) * σ_ch(:e, :g, k + 1)
        )
            for k in 1:(N - 1)
    )
    c_ops = [σ_ch(:g, :e, k) for k in 1:N]
    eqs = meanfield(σ_ch(:g, :e, 1), H, c_ops; rates = [γ for _ in 1:N], order = 2)
    complete!(eqs)
    @test length(eqs.equations) == 65
    sys_c = mtkcompile(System(eqs; name = :chain4_sys))
    @test length(unknowns(sys_c)) == 65
    u0 = zeros(ComplexF64, length(eqs.equations))
    init = initial_values(eqs, u0)
    ps = Dict(
        Ω => 2.0, γ => 1.0, Δ => 0.0, J0 => 1.0,
        x_[1] => 0.0, x_[2] => 1.0, x_[3] => 2.0, x_[4] => 3.0,
    )
    prob = ODEProblem(sys_c, merge(init, ps), (0.0, 15.0))
    sol = solve(prob, Tsit5(); reltol = 1.0e-9, abstol = 1.0e-9)
    @test isapprox(
        real(get_solution(sol, σ_ch(:e, :e, N), eqs).(sol.t[end])),
        0.09034045248945875; rtol = 1.0e-4,
    )
end

@testset "superradiant_laser_indexed" begin
    hc = FockSpace(:cavity); ha = NLevelSpace(:atom, 2); h = hc ⊗ ha
    @qnumbers a::Destroy(h)
    σ(α, β, idx) = IndexedOperator(Transition(h, :σ, α, β), idx)
    @variables N Δ κ Γ R ν
    g_v(idx) = IndexedVariable(:g, idx)
    i = Index(h, :i, N, ha); j = Index(h, :j, N, ha)
    H = -Δ * a' * a + Σ(g_v(i) * (a' * σ(1, 2, i) + a * σ(2, 1, i)), i)
    J = [a, σ(1, 2, i), σ(2, 1, i), σ(2, 2, i)]
    rates = [κ, Γ, R, ν]
    ops = [a' * a, σ(2, 2, j)]
    eqs = meanfield(ops, H, J; rates = rates, order = 2)
    eqs_c = complete(eqs; filter_func = phase_invariant)
    @test length(eqs_c.equations) == 6
    eqs_sc = scale(eqs_c)
    @test length(eqs_sc.equations) == 5
    sys_c = mtkcompile(System(eqs_sc; name = :sr_sys))
    @test length(unknowns(sys_c)) == 5
    u0 = zeros(ComplexF64, length(eqs_sc.equations))
    init = initial_values(eqs_sc, u0)
    pmap = parameter_map(eqs_sc, Dict(
        N => 50.0, Δ => 0.0, g_v(i) => 1.0, κ => 1.0, Γ => 0.25,
        R => 4.0, ν => 0.0,
    ))
    prob = ODEProblem(sys_c, merge(init, pmap), (0.0, 30.0))
    sol = solve(prob, Tsit5(); reltol = 1.0e-8, abstol = 1.0e-8)
    @test isapprox(
        real(get_solution(sol, a' * a, eqs_sc).(sol.t[end])),
        91.96090096151146; rtol = 1.0e-4,
    )
    corr = CorrelationFunction(
        a', a, eqs_c; steady_state = true, filter_func = phase_invariant,
    )
    corr_sc = scale(corr)
    @test length(corr_sc.eqs.equations) >= 1
    ps_tup = (Δ, κ, Γ, R, ν, N, g_v(i))
    p0 = (0.0, 1.0, 0.25, 4.0, 0.0, 50.0, 1.0)
    S = Spectrum(corr_sc, ps_tup)
    s_vals = S([-0.5, 0.0, 0.5], sol.u[end], p0)
    @test all(isfinite, s_vals)
    # spectrum symmetry S(ω) = S(-ω)
    @test isapprox(s_vals[1], s_vals[3]; rtol = 1.0e-6)
end

@testset "unique_squeezing" begin
    hf = FockSpace(:harmonic); ha = NLevelSpace(:spin, 2); h = hf ⊗ ha
    @variables ω Ω ωd η κ g γ N ξ
    @qnumbers a::Destroy(h)
    σ(x, y, idx) = IndexedOperator(Transition(h, :σ, x, y), idx)
    b = a * cosh(ξ) + a' * sinh(ξ)
    i = Index(h, :i, N, ha); j = Index(h, :j, N, ha)
    Hf = ω * a' * a + η * (b' * exp(-1im * ωd * t) + b * exp(1im * ωd * t))
    Ha = Ω * Σ(σ(2, 2, i) - σ(1, 1, i), i) / 2
    Hi = g * Σ((σ(1, 2, i) + σ(2, 1, i)) * (a + a'), i) / 2
    H = Hf + Ha + Hi
    J = [b, σ(1, 2, i)]
    rates = [κ, γ]
    eqs = meanfield([a, a' * a, σ(2, 2, j)], H, J; rates = rates, order = 2)
    eqs_c = complete(eqs)
    @test length(eqs_c.equations) == 22
    eqs_sc = scale(eqs_c)
    @test length(eqs_sc.equations) == 19
    sys_c = mtkcompile(System(eqs_sc; name = :us_sys))
    @test length(unknowns(sys_c)) == 19
    ω_ = 1.0; Ω_ = 2.0e3; η_ = 4.0; κ_ = 1.0; γ_ = 1.0
    N_ = 1
    gc_ = sqrt(Ω_ * ω_ / N_)
    g_ = 0.9 * gc_
    ωd_ = sqrt(1 - g_^2 / gc_^2) * ω_
    ξ_ = (1 / 4) * log(1 - N_ * g_^2 / (ω_ * Ω_))
    u0 = zeros(ComplexF64, length(eqs_sc.equations))
    init = initial_values(eqs_sc, u0)
    ps = Dict(
        ω => ω_, Ω => Ω_, ωd => ωd_, g => g_, η => η_,
        κ => κ_, γ => γ_, N => Float64(N_), ξ => ξ_,
    )
    tend = Float64(2π / ωd_)
    prob = ODEProblem(sys_c, merge(init, ps), (0.0, tend))
    sol = solve(prob, Tsit5(); reltol = 1.0e-8, abstol = 1.0e-8)
    @test isapprox(
        real(get_solution(sol, a' * a, eqs_sc).(sol.t[end])),
        30.149750050640534; rtol = 1.0e-3,
    )
end

@testset "unique_squeezing (free-j shape, N=100 plateau)" begin
    hf = FockSpace(:harmonic); ha = NLevelSpace(:spin, 2); h = hf ⊗ ha
    @variables ω Ω ωd η κ g γ N ξ
    @qnumbers a::Destroy(h)
    σ(x, y, idx) = IndexedOperator(Transition(h, :σ, x, y), idx)
    b = a * cosh(ξ) + a' * sinh(ξ)
    i = Index(h, :i, N, ha); j = Index(h, :j, N, ha)
    Hf = ω * a' * a + η * (b' * exp(-1im * ωd * t) + b * exp(1im * ωd * t))
    Ha = Ω * Σ(σ(2, 2, i) - σ(1, 1, i), i) / 2
    Hi = g * Σ((σ(1, 2, i) + σ(2, 1, i)) * (a + a'), i) / 2
    H = Hf + Ha + Hi
    J = [b, σ(1, 2, i)]
    rates = [κ, γ]
    eqs = meanfield([a, a' * a, σ(2, 2, j)], H, J; rates = rates, order = 2)
    eqs_c = complete(eqs)
    @test length(eqs_c.equations) == 22
    eqs_sc = scale(eqs_c)
    @test length(eqs_sc.equations) == 19
    sys_c = mtkcompile(System(eqs_sc; name = :us_freej_sys))
    @test length(unknowns(sys_c)) == 19
    ω_ = 1.0; Ω_ = 2.0e3; η_ = 4.0; κ_ = 1.0; γ_ = 1.0
    N_ = 100
    gc_ = sqrt(Ω_ * ω_ / N_)
    g_ = 0.9 * gc_
    ωd_ = sqrt(1 - g_^2 / gc_^2) * ω_
    ξ_ = (1 / 4) * log(1 - N_ * g_^2 / (ω_ * Ω_))
    u0 = zeros(ComplexF64, length(eqs_sc.equations))
    u0d = Dict{Any, Any}(unknowns(sys_c) .=> u0)
    ps = Dict{Any, Any}(
        ω => ω_, Ω => Ω_, ωd => ωd_, g => g_, η => η_,
        κ => κ_, γ => γ_, N => Float64(N_), ξ => ξ_,
    )
    tend = Float64(4π / ωd_)
    prob = ODEProblem(sys_c, merge(u0d, ps), (0.0, tend))
    sol = solve(
        prob, Tsit5(); saveat = π / 30ωd_,
        reltol = 1.0e-10, abstol = 1.0e-10,
    )
    t_ = sol.t
    adag_a = get_solution(sol, a' * a, eqs_sc).(t_)
    aa = get_solution(sol, a * a, eqs_sc).(t_)
    adag_adag = get_solution(sol, a' * a', eqs_sc).(t_)
    a_ = get_solution(sol, a, eqs_sc).(t_)
    adag = get_solution(sol, a', eqs_sc).(t_)
    sqx = real.(adag_adag + aa + 2 * adag_a .+ 1 - (adag + a_) .^ 2)
    sqy = real.(adag_adag + aa - 2 * adag_a .- 1 - (adag - a_) .^ 2)
    # N=100 squeezing plateau, Gietka et al., PRL 131, 223604: X ≈ 2.29, P ≈ 0.44.
    @test isapprox(sqx[end], 2.292622635966603; rtol = 1.0e-3)
    @test isapprox(-sqy[end], 0.43703195626758884; rtol = 1.0e-3)
end

@testset "heterodyne_detection" begin
    @variables N ωa γ η χ ωc κ g ξ ωl
    @register_symbolic _het_pulse(tt)
    hc = FockSpace(:resonator); ha = NLevelSpace(:atom, 2); h = hc ⊗ ha
    j = Index(h, :j, N, ha); k = Index(h, :k, N, ha)
    @qnumbers a::Destroy(h, 1)
    σ(α, β, idx) = IndexedOperator(Transition(h, :σ, α, β, 2), idx)
    eqs_seed = meanfield(a, ωc * a' * a, [a]; rates = [κ])
    tv = eqs_seed.iv
    H = ωc * a' * a + ωa * Σ(σ(2, 2, j), j) +
        g * a' * Σ(σ(1, 2, j), j) +
        g * a * Σ(σ(2, 1, j), j)
    J = [a * exp(1.0im * ωl * tv), σ(1, 2, j), σ(2, 1, j), σ(2, 2, j)]
    rates = [κ, γ, η * _het_pulse(tv), 2 * χ]
    efficiencies = [ξ, 0, 0, 0]
    ops = [a', a' * a, σ(2, 2, k), σ(1, 2, k), a * a]
    eqs = meanfield(
        ops, H, J;
        rates = rates, efficiencies = efficiencies,
        direction = Forward(), order = 2, iv = tv,
    )
    eqs_c = complete(eqs; get_adjoints = false)
    @test length(eqs_c.equations) == 13
    eqs_sc = scale(eqs_c)
    @test length(eqs_sc.equations) == 12
    sys_c = mtkcompile(System(eqs_sc; name = :het_sys))
    @test length(unknowns(sys_c)) == 12
end

@testset "heterodyne_detection (single-trajectory SDE bounded)" begin
    using StochasticDiffEq: SDEProblem, EM, RealWienerProcess
    using StochasticDiffEq.SciMLBase: ReturnCode
    import Random
    @variables N ωa γ η χ ωc κ g ξ ωl
    @register_symbolic _het_pulse2(tt)
    hc = FockSpace(:resonator); ha = NLevelSpace(:atom, 2); h = hc ⊗ ha
    j = Index(h, :j, N, ha); k = Index(h, :k, N, ha)
    @qnumbers a::Destroy(h, 1)
    σ(α, β, idx) = IndexedOperator(Transition(h, :σ, α, β, 2), idx)
    eqs_seed = meanfield(a, ωc * a' * a, [a]; rates = [κ])
    tv = eqs_seed.iv
    H = ωc * a' * a + ωa * Σ(σ(2, 2, j), j) +
        g * a' * Σ(σ(1, 2, j), j) +
        g * a * Σ(σ(2, 1, j), j)
    J = [a * exp(1.0im * ωl * tv), σ(1, 2, j), σ(2, 1, j), σ(2, 2, j)]
    rates = [κ, γ, η * _het_pulse2(tv), 2 * χ]
    efficiencies = [ξ, 0, 0, 0]
    ops = [a', a' * a, σ(2, 2, k), σ(1, 2, k), a * a]
    eqs = meanfield(
        ops, H, J;
        rates = rates, efficiencies = efficiencies,
        direction = Forward(), order = 2, iv = tv,
    )
    eqs_c = complete(eqs; get_adjoints = false)
    eqs_sc = scale(eqs_c)
    sys_st = mtkcompile(System(eqs_sc; name = :hetst_sys))
    ωc_ = 0.0; κ_ = 2π * 1.13e3; ξ_ = 0.12; N_ = 5.0e4
    ωa_ = 0.0; γ_ = 2π * 0.375; η_ = 2π * 20; χ_ = 0.016
    g_ = 2π * 0.73; ωl_ = 2π * 1.0e3
    t0 = 0.0; t1 = 20.0e-3
    _het_pulse2(tt) = (tt > t0 && tt < t0 + t1) * 1.0
    T_end = 0.1
    p_pairs = Dict(
        N => N_, ωa => ωa_, γ => γ_, η => η_, χ => χ_,
        ωc => ωc_, κ => κ_, g => g_, ξ => ξ_, ωl => ωl_,
    )
    u0_vec = zeros(ComplexF64, length(eqs_sc))
    dict_st = parameter_map(
        sys_st,
        merge(Dict{Any, Any}(unknowns(sys_st) .=> u0_vec),
            Dict{Any, Any}(p_pairs)),
    )
    Random.seed!(2)
    noise = RealWienerProcess(0.0, 0.0)
    prob_st = SDEProblem(sys_st, dict_st, (0.0, T_end); noise = noise)
    sol = solve(prob_st, EM(); dt = T_end / 2.0e5)
    # The SDE must stay bounded; example trajectory peaks are ~600.
    adag_a_traj = real.(get_solution(sol, a' * a, eqs_sc).(sol.t))
    @test sol.retcode == ReturnCode.Success
    @test all(isfinite, adag_a_traj)
    @test maximum(abs, adag_a_traj) < 1.0e8
end

@testset "heterodyne_detection (drift-only numerical)" begin
    @variables N ωa γ η χ ωc κ g ωl
    @register_symbolic _het_pulse(tt)
    hc = FockSpace(:resonator); ha = NLevelSpace(:atom, 2); h = hc ⊗ ha
    j = Index(h, :j, N, ha); k = Index(h, :k, N, ha)
    @qnumbers a::Destroy(h, 1)
    σ(α, β, idx) = IndexedOperator(Transition(h, :σ, α, β, 2), idx)
    eqs_seed = meanfield(a, ωc * a' * a, [a]; rates = [κ])
    tv = eqs_seed.iv
    H = ωc * a' * a + ωa * Σ(σ(2, 2, j), j) +
        g * a' * Σ(σ(1, 2, j), j) +
        g * a * Σ(σ(2, 1, j), j)
    J = [a * exp(1.0im * ωl * tv), σ(1, 2, j), σ(2, 1, j), σ(2, 2, j)]
    rates = [κ, γ, η * _het_pulse(tv), 2 * χ]
    ops = [a', a' * a, σ(2, 2, k), σ(1, 2, k), a * a]
    eqs = meanfield(
        ops, H, J;
        rates = rates, direction = Forward(), order = 2, iv = tv,
    )
    eqs_c = complete(eqs; get_adjoints = false)
    eqs_sc = scale(eqs_c)
    sys_c = mtkcompile(System(eqs_sc; name = :het_det_sys))
    _het_pulse(_) = 1.0
    u0 = zeros(ComplexF64, length(eqs_sc.equations))
    init = initial_values(eqs_sc, u0)
    pmap = parameter_map(eqs_sc, Dict(
        ωc => 0.0, ωa => 0.0, ωl => 0.0,
        κ => 1.0, γ => 0.1, η => 1.0, χ => 0.1, g => 1.0, N => 4.0,
    ))
    prob = ODEProblem(sys_c, merge(init, pmap), (0.0, 5.0))
    sol = solve(prob, Tsit5(); reltol = 1.0e-8, abstol = 1.0e-8)
    @test isapprox(
        real(get_solution(sol, a' * a, eqs_sc).(sol.t[end])),
        1.906585813419867; rtol = 1.0e-4,
    )
    @test isapprox(
        real(get_solution(sol, eqs_sc.states[3], eqs_sc).(sol.t[end])),
        0.4702676694506166; rtol = 1.0e-4,
    )
end

@testset "cavity_antiresonance_indexed" begin
    hc = FockSpace(:cavity); ha = NLevelSpace(:atom, 2); h = hc ⊗ ha
    @variables N Δc η Δa κ
    g_v(idx) = IndexedVariable(:g, idx)
    Γ_v(i, j) = DoubleIndexedVariable(:Γ, i, j)
    Ω_v(i, j) = DoubleIndexedVariable(:Ω, i, j; identical = false)
    i = Index(h, :i, N, ha); j = Index(h, :j, N, ha)
    @qnumbers a::Destroy(h)
    σ(x, y, idx) = IndexedOperator(Transition(h, :σ, x, y), idx)
    Hc = Δc * a' * a + η * (a' + a)
    Ha = Δa * Σ(σ(2, 2, i), i) +
        Σ(Σ(Ω_v(i, j) * σ(2, 1, i) * σ(1, 2, j), j, [i]), i)
    Hi = Σ(g_v(i) * (a' * σ(1, 2, i) + a * σ(2, 1, i)), i)
    H = Hc + Ha + Hi
    J = [a, σ(1, 2, i)]
    rates = [κ, Γ_v(i, j)]
    eqs = meanfield(a, H, J; rates = rates, order = 1)
    complete!(eqs)
    @test length(eqs.equations) == 4
    eqs_ = evaluate(eqs; limits = (N => 2))
    @test length(eqs_.equations) == 5
    sys_c = mtkcompile(System(eqs_; name = :ca_sys))
    @test length(unknowns(sys_c)) == 5
    u0 = zeros(ComplexF64, length(eqs_.equations))
    init = initial_values(eqs_, u0)
    pmap_dict = Dict{Any, Any}()
    for p in ModelingToolkitBase.parameters(sys_c)
        nm = string(p)
        if nm == "Δc";    pmap_dict[p] = 0.0
        elseif nm == "Δa"; pmap_dict[p] = 0.0
        elseif nm == "η";  pmap_dict[p] = 0.1
        elseif nm == "κ";  pmap_dict[p] = 1.0
        elseif nm == "g";  pmap_dict[p] = [1.0, 1.0]
        elseif nm == "Γ";  pmap_dict[p] = [1.0 1.0; 1.0 1.0]
        elseif nm == "Ω";  pmap_dict[p] = [0.0 0.5; 0.5 0.0]
        end
    end
    prob = ODEProblem(sys_c, merge(init, pmap_dict), (0.0, 30.0))
    sol = solve(prob, Tsit5(); reltol = 1.0e-9, abstol = 1.0e-9)
    @test isapprox(
        abs2(get_solution(sol, a, eqs_).(sol.t[end])),
        0.0019889934766074693; rtol = 1.0e-4,
    )
end

@testset "filter-cavity_indexed" begin
    @variables κ g gf κf R Γ Δ ν N M
    δ_v(idx) = IndexedVariable(:δ, idx)
    hc = FockSpace(:cavity); hf = FockSpace(:filter); ha = NLevelSpace(:atom, 2)
    h = hc ⊗ hf ⊗ ha
    i = Index(h, :i, M, hf); j = Index(h, :j, N, ha)
    @qnumbers a::Destroy(h, 1)
    b(idx) = IndexedOperator(Destroy(h, :b, 2), idx)
    b(idx::Integer) = IndexedOperator(Destroy(h, :b, 2), i(2)(idx))
    σ(α, β, idx) = IndexedOperator(Transition(h, :σ, α, β, 3), idx)
    σ(α, β, idx::Integer) = IndexedOperator(Transition(h, :σ, α, β, 3), j(idx))
    H = Δ * Σ(σ(2, 2, j), j) +
        Σ(δ_v(i) * b(i)' * b(i), i) +
        gf * (Σ(a' * b(i) + a * b(i)', i)) +
        g * (Σ(a' * σ(1, 2, j) + a * σ(2, 1, j), j))
    J = [a, b(i), σ(1, 2, j), σ(2, 1, j), σ(2, 2, j)]
    rates = [κ, κf, Γ, R, ν]
    eqs = meanfield(a' * a, H, J; rates = rates, order = 2)
    eqs_c = complete!(deepcopy(eqs))
    # closure: 42 equations, 39 after scale, 554 after evaluate at M=20.
    @test length(eqs_c.equations) == 42
    eqs_sc = scale(eqs_c; h = [3])
    @test length(eqs_sc.equations) == 39
    eqs_eval = evaluate(eqs_sc; limits = Dict(M => 20))
    @test length(eqs_eval.equations) == 554
    sys_c = mtkcompile(System(eqs_eval; name = :fc_sys))
    @test length(unknowns(sys_c)) == 554
end

# Reduced M=3 variant integrated as a mean-field ODE.
@testset "filter-cavity_indexed (M=3 numerical)" begin
    @variables κ g gf κf R Γ Δ ν N M
    δ_v(idx) = IndexedVariable(:δ, idx)
    hc = FockSpace(:cavity); hf = FockSpace(:filter); ha = NLevelSpace(:atom, 2)
    h = hc ⊗ hf ⊗ ha
    i = Index(h, :i, M, hf); j = Index(h, :j, N, ha)
    @qnumbers a::Destroy(h, 1)
    b(idx) = IndexedOperator(Destroy(h, :b, 2), idx)
    b(idx::Integer) = IndexedOperator(Destroy(h, :b, 2), i(2)(idx))
    σ(α, β, idx) = IndexedOperator(Transition(h, :σ, α, β, 3), idx)
    σ(α, β, idx::Integer) = IndexedOperator(Transition(h, :σ, α, β, 3), j(idx))
    H = Δ * Σ(σ(2, 2, j), j) +
        Σ(δ_v(i) * b(i)' * b(i), i) +
        gf * (Σ(a' * b(i) + a * b(i)', i)) +
        g * (Σ(a' * σ(1, 2, j) + a * σ(2, 1, j), j))
    J = [a, b(i), σ(1, 2, j), σ(2, 1, j), σ(2, 2, j)]
    rates = [κ, κf, Γ, R, ν]
    eqs = meanfield(a' * a, H, J; rates = rates, order = 2)
    eqs_c = complete!(deepcopy(eqs))
    eqs_sc = scale(eqs_c; h = [3])
    M_ = 3
    eqs_eval = evaluate(eqs_sc; limits = Dict(M => M_))
    @test length(eqs_eval.equations) == 44
    sys_c = mtkcompile(System(eqs_eval; name = :fc3_sys))
    u0 = zeros(ComplexF64, length(eqs_eval.equations))
    init = initial_values(eqs_eval, u0)
    pmap_dict = Dict{Any, Any}(
        κ => 1.0, κf => 0.05, g => 0.5, gf => 0.05,
        R => 1.0, Γ => 0.01, Δ => 0.0, ν => 0.01, N => 10.0,
    )
    δ_vals = [-1.0, 0.0, 1.0]
    for k in 1:M_
        pmap_dict[δ_v(i(2)(k))] = δ_vals[k]
    end
    pmap = parameter_map(eqs_eval, pmap_dict)
    prob = ODEProblem(sys_c, merge(init, pmap), (0.0, 5.0))
    sol = solve(prob, Tsit5(); reltol = 1.0e-8, abstol = 1.0e-8)
    @test isapprox(
        real(get_solution(sol, a' * a, eqs_eval).(sol.t[end])),
        4.406206739122889; rtol = 1.0e-4,
    )
end

@testset "retrodiction_homodyne (forward, noise)" begin
    using StochasticDiffEq: SDEProblem, EM, RealWienerProcess
    using StochasticDiffEq.SciMLBase: ReturnCode
    import Random
    h = PhaseSpace(:motion)
    @qnumbers x::Position(h)
    p = Momentum(h, :p)
    @variables Ω Γ η s
    m = 1
    a = (x + 1im * p) * s
    H = p^2 / (2m) + 0.5m * Ω^2 * x^2
    eqs = meanfield(
        [x, p, x * x, x * p, p * p], H, [a];
        rates = [Γ], efficiencies = [η], order = 2,
    )
    @test length(eqs.equations) == 5
    sys_c = mtkcompile(System(eqs; name = :ret_sys))
    @test length(unknowns(sys_c)) == 5
    u0 = ComplexF64[5.0, 0.0, 5 + 25.0, 0.5im, 5.0]
    init = initial_values(eqs, u0)
    pmap = Dict(Ω => 1.0, Γ => 1 / 6, η => 0.5, s => 1 / √2)
    Random.seed!(11)
    noise = RealWienerProcess(0.0, 0.0)
    prob_st = SDEProblem(sys_c, merge(init, pmap), (0.0, 2.0); noise = noise)
    sol = solve(prob_st, EM(); dt = 1.0e-4)
    x_traj = real.(get_solution(sol, x, eqs).(sol.t))
    @test all(isfinite, x_traj)
    @test maximum(abs, x_traj) < 1.0e6
end

# Forward drift only (no measurement noise), solved as an ODE over the
# PhaseSpace / Position / Momentum operator path.
@testset "retrodiction_homodyne (forward drift only)" begin
    h = PhaseSpace(:motion)
    @qnumbers x::Position(h)
    p = Momentum(h, :p)
    @variables Ω Γ s
    m = 1
    a = (x + 1im * p) * s
    H = p^2 / (2m) + 0.5m * Ω^2 * x^2
    eqs = meanfield(
        [x, p, x * x, x * p, p * p], H, [a];
        rates = [Γ], order = 2,
    )
    sys_c = mtkcompile(System(eqs; name = :ret_drift_sys))
    u0 = ComplexF64[5.0, 0.0, 5 + 25.0, 0.5im, 5.0]
    init = initial_values(eqs, u0)
    pmap = Dict(Ω => 1.0, Γ => 1 / 6, s => 1 / √2)
    prob = ODEProblem(sys_c, merge(init, pmap), (0.0, 5.0))
    sol = solve(prob, Tsit5(); reltol = 1.0e-9, abstol = 1.0e-9)
    @test isapprox(
        real(get_solution(sol, x, eqs).(sol.t[end])),
        0.935008189535323; rtol = 1.0e-5,
    )
    @test isapprox(
        real(get_solution(sol, p, eqs).(sol.t[end])),
        3.1608092157050347; rtol = 1.0e-5,
    )
    @test isapprox(
        real(get_solution(sol, x * x, eqs).(sol.t[end])),
        3.3299322539872325; rtol = 1.0e-5,
    )
end
