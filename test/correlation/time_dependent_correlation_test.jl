using QuantumCumulants
using Symbolics: Symbolics, @variables, @register_symbolic
using SymbolicUtils
using ModelingToolkitBase: ModelingToolkitBase, @named, mtkcompile, ODEProblem, unknowns
using OrdinaryDiffEqTsit5: Tsit5, solve, ReturnCode
using Test

# Registered time-dependent coefficient; constant ≡ 1 so the τ-equations reduce to the
# constant-coefficient damped cavity and we can check against the analytic result.
@register_symbolic _corr_f(tt)
_corr_f(_) = 1.0

# Gaussian pulse centred at t = 5 (issue #93 "pulse not re-applied" test).
@register_symbolic _pulse_f(tt)
_pulse_f(tt) = exp(-(tt - 5.0)^2)

# Five independent drive envelopes for the multi-mode reproduction (issue #171).
@register_symbolic _f171_1(tt)
@register_symbolic _f171_2(tt)
@register_symbolic _f171_3(tt)
@register_symbolic _f171_4(tt)
@register_symbolic _f171_5(tt)
_f171_1(tt) = exp(-(tt - 3.0)^2)
_f171_2(tt) = exp(-(tt - 4.0)^2)
_f171_3(tt) = exp(-(tt - 5.0)^2)
_f171_4(tt) = exp(-(tt - 6.0)^2)
_f171_5(tt) = exp(-(tt - 7.0)^2)

@testset "time-dependent Hamiltonian: correlation ODE solves (issue #171)" begin
    hc = FockSpace(:cavity)
    a = Destroy(hc, :a)
    @variables ωc::Real κ::Real
    @variables t0::Real

    # Seed run to grab the parent independent variable, then build the time-dependent H.
    seed = meanfield(a' * a, ωc * a' * a, [a]; rates = [κ])
    tv = seed.iv
    H = _corr_f(tv) * ωc * a' * a
    he = meanfield(a' * a, H, [a]; rates = [κ])
    complete!(he)

    # Omitting iv0 on a time-dependent system is a user error with a clear message.
    @test_throws ArgumentError CorrelationFunction(a, a', he)

    # t -> t0 + τ substitution; with _corr_f ≡ 1 the dynamics are the constant case.
    c = CorrelationFunction(a, a', he; iv0 = t0)
    @test c isa CorrelationFunction

    ps = (ωc, κ)
    pn = ComplexF64[1.0, 0.5]
    u_end = Dict(SymbolicUtils.unwrap(average(a' * a)) => 0.0 + 0im)
    u0 = correlation_u0(c, u_end)
    # t0 is an ordinary parameter of the τ-system; supply it like any other.
    p0_dict = correlation_p0(c, u_end, [ωc => pn[1], κ => pn[2], t0 => 0.0])

    @named csys = System(c)
    csys_c = mtkcompile(csys)
    prob = ODEProblem(csys_c, merge(u0, p0_dict), (0.0, 30.0))
    sol = solve(prob, Tsit5(); abstol = 1.0e-10, reltol = 1.0e-10)
    @test sol.retcode == ReturnCode.Success

    g_var = first(unknowns(csys_c))
    τs = [0.5, 1.0, 2.0, 4.0]
    g_num = [sol(τ; idxs = g_var) for τ in τs]
    g_an = @. exp(-(im * pn[1] + pn[2] / 2) * τs)
    @test maximum(abs.(g_num .- g_an)) < 1.0e-7
end

@testset "time-dependent Hamiltonian: pulse is not re-applied at τ=0 (issue #93)" begin
    # Gaussian-pulse-driven cavity. The correlation starts at t0 = 10, long after the
    # pulse (centred at 5) has passed. With the correct t -> t0 + τ substitution the drive
    # seen during the delay is ≈ 0, so the field merely decays:
    #   g(τ) = ⟨a†(t0+τ) a(t0)⟩ = ⟨a†a⟩(t0) exp(-κτ/2).
    # A naive t -> τ would re-fire the pulse at τ=0 and break this.
    h = FockSpace(:cavity)
    a = Destroy(h, :a)
    @variables Ω0::Real κ::Real
    @variables t0::Real

    seed = meanfield(a, Ω0 * (a + a'), [a]; rates = [κ])
    tv = seed.iv
    H = Ω0 * _pulse_f(tv) * (a + a')   # no bare frequency: only the (vanishing) drive remains
    he = meanfield([a, a' * a], H, [a]; rates = [κ])
    complete!(he)

    c = CorrelationFunction(a', a, he; iv0 = t0)

    Ω0v, κv, tstop = 2.0, 1.0, 10.0
    sys = mtkcompile(System(he; name = :pulse_sys))
    u0 = initial_values(he, zeros(ComplexF64, length(he)))
    sol = solve(
        ODEProblem(sys, merge(u0, Dict(Ω0 => Ω0v, κ => κv)), (0.0, tstop)),
        Tsit5(); abstol = 1.0e-10, reltol = 1.0e-10,
    )
    @test sol.retcode == ReturnCode.Success

    csys = mtkcompile(System(c; name = :pulse_corr))
    u0_c = correlation_u0(c, sol.u[end])
    p0_c = correlation_p0(c, sol.u[end], [Ω0 => Ω0v, κ => κv, t0 => tstop])
    sol_c = solve(
        ODEProblem(csys, merge(u0_c, Dict(p0_c)), (0.0, 8.0)),
        Tsit5(); abstol = 1.0e-10, reltol = 1.0e-10, save_idxs = 1,
    )
    @test sol_c.retcode == ReturnCode.Success

    # Pure free decay from g(0): confirms the pulse is not re-applied during the delay.
    g0 = sol_c.u[1]
    τs = [0.5, 1.0, 2.0, 4.0]
    g_num = [sol_c(τ)[1] for τ in τs]
    g_decay = @. g0 * exp(-κv / 2 * τs)
    @test maximum(abs.(g_num .- g_decay)) < 1.0e-6
end

@testset "time-dependent Hamiltonian: multi-mode correlation solves (issue #171)" begin
    # Faithful reproduction of the three-cavity, five-drive system from issue #171, which
    # crashed with `UndefVarError: t not defined`. Regression guard: it must now solve.
    h1 = FockSpace(:cavity1); h2 = FockSpace(:cavity2); h3 = FockSpace(:cavity3)
    h = h1 ⊗ h2 ⊗ h3
    a1 = Destroy(h, :a1, 1); a2 = Destroy(h, :a2, 2); a3 = Destroy(h, :a3, 3)
    @variables η_1::Real η_3::Real θ_1::Real θ_2::Real g_2::Real κ::Real
    @variables s0::Real

    seed = meanfield(a1, η_1 * (a1 + a1'), [a1, a2, a3]; rates = [κ, κ, κ])
    t = seed.iv
    H = -_f171_1(t) * η_3 * (a3 + a3') + _f171_2(t) * θ_2 * (a2' * a3 + a2 * a3') +
        _f171_3(t) * g_2 * (a2 * a2 + a2' * a2') + _f171_4(t) * θ_1 * (a1' * a2 + a1 * a2') +
        _f171_5(t) * η_1 * (a1 + a1')
    me = meanfield([a1, a2, a3], H, [a1, a2, a3]; rates = [κ, κ, κ], order = 2)
    complete!(me)

    c = CorrelationFunction(a1', a1, me; iv0 = s0)   # ⟨a₁†(t0+τ) a₁(t0)⟩
    @test c isa CorrelationFunction

    ps = [η_1 => 0.5, η_3 => 0.4, θ_1 => 0.3, θ_2 => 0.2, g_2 => 0.1, κ => 1.0]
    sys = mtkcompile(System(me; name = :n171))
    u0 = initial_values(me, zeros(ComplexF64, length(me)))
    sol = solve(
        ODEProblem(sys, merge(u0, Dict(ps)), (0.0, 10.0)),
        Tsit5(); maxiters = Int(1.0e7),
    )
    @test sol.retcode == ReturnCode.Success

    csys = mtkcompile(System(c; name = :cc171))
    u0_c = correlation_u0(c, sol.u[end])
    p0_c = correlation_p0(c, sol.u[end], vcat(ps, [s0 => 10.0]))
    sol_c = solve(
        ODEProblem(csys, merge(u0_c, Dict(p0_c)), (0.0, 50.0)),   # 5*stop_time, as in the issue
        Tsit5(); save_idxs = 1,
    )
    @test sol_c.retcode == ReturnCode.Success
    @test all(isfinite, abs.(sol_c.u))
end
