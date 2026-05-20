using QuantumCumulants
using Symbolics: Symbolics, @variables, @register_symbolic
using ModelingToolkitBase: mtkcompile, unknowns
using OrdinaryDiffEq: Tsit5, ODEProblem, ReturnCode, solve
using StochasticDiffEq: StochasticDiffEq, SDEProblem, EM, RealWienerProcess
using Statistics: cor
using Random
using Test

# Continuous monitoring of a harmonic oscillator: simulate the forward
# SDE, extract the measurement record dY from the realised Wiener path,
# then build the deterministic Kalman ODE that is driven by dY via
# `modify_equations`. On the *same* noise realisation, the Kalman ODE
# must track the SDE trajectory (perfect agreement modulo integration
# scheme: EM vs RK; we assert correlation >= 0.95 and relative endpoint
# error < 0.2). Split out of measurement_retrodiction_test.jl: this
# testset's SDE forward solve + MTK compile + ODE solve dominates the
# file's wall-clock, isolating it lets ParallelTestRunner schedule it on
# its own worker.
@testset "Kalman: forward SDE tracked by measurement-record-driven ODE" begin
    h = PhaseSpace(:motion)
    @qnumbers x::Position(h)
    p = Momentum(h, :p)
    @variables Ω Γ η s
    m = 1
    a = (x + 1im * p) * s
    H = p^2 / (2m) + 0.5m * Ω^2 * x^2
    J = [a]
    R = [Γ]
    eff = [η]
    ops = [x, p, x * x, x * p, p * p]

    eqs = meanfield(ops, H, J; rates = R, efficiencies = eff, order = 2)
    @test eqs isa NoiseMeanFieldEquations

    Ω_, Γ_, η_ = 1.0, 1 / 6, 0.5
    Tend = 0.5 / Γ_
    dt = Tend / 2.0e3
    T_saveat = collect(0:dt:Tend)
    ps = [Ω, Γ, η, s]
    pn = [Ω_, Γ_, η_, 1 / √2]
    x0, p0 = 5.0, 0.0
    u0 = ComplexF64[x0, p0, 5 + x0 * x0, im / 2 + x0 * p0, 5 + p0 * p0]

    Random.seed!(2)
    sys_fw = mtkcompile(System(eqs; name = :sys_fw))
    dict_fw = merge(Dict(unknowns(sys_fw) .=> u0), Dict(ps .=> pn))
    noise = StochasticDiffEq.RealWienerProcess(0.0, 0.0; save_everystep = true)
    prob_fw = SDEProblem(sys_fw, dict_fw, (0, Tend); noise = noise)
    sol_fw = solve(prob_fw, EM(); dt = dt, saveat = T_saveat)
    @test sol_fw.retcode == ReturnCode.Success
    @test all(isfinite.(real.(get_solution(sol_fw, x, eqs).(T_saveat))))

    t_W = noise.t
    W_W = noise.W
    Nstep = length(t_W) - 1
    sol_x_at_W = real(get_solution(sol_fw, x, eqs).(t_W))
    dY_W = zeros(Nstep)
    for k in 1:Nstep
        dW = W_W[k + 1] - W_W[k]
        cd_p_c = sol_x_at_W[k] * √(2) * √(Γ_)
        dY_W[k] = dW + √η_ * cd_p_c * dt
    end
    Ydot_data = [dY_W; dY_W[end]] / dt
    function Ydot_impl(time::Real)
        k = Int(round(time / dt)) + 1
        k > length(Ydot_data) && return Ydot_data[end]
        k < 1 && return Ydot_data[1]
        return Ydot_data[k]
    end
    @register_symbolic Ydot_meas(time)
    Ydot_meas(time::Real) = Ydot_impl(time)

    eqs_det = meanfield(ops, H, J; rates = R, order = 2)
    iv = eqs_det.iv
    function f_measure(lhs, rhs)
        avg_apad = average(a + a')
        inn = lhs * a + a' * lhs - avg_apad * lhs
        drive = Ydot_meas(iv) - √(η * Γ) * avg_apad
        term_ = √(η * Γ) * inn * drive
        return rhs + cumulant_expansion(average(term_), 2)
    end
    eqs_kal = modify_equations(eqs_det, f_measure)
    @test eqs_kal isa MeanFieldEquations
    @test all(
        !_is_zero(eqs_kal.equations[i].rhs - eqs_det.equations[i].rhs)
            for i in eachindex(eqs_kal.equations)
    )

    sys_kal = mtkcompile(System(eqs_kal; name = :sys_kal))
    dict_kal = merge(Dict(unknowns(sys_kal) .=> u0), Dict(ps .=> pn))
    prob_kal = ODEProblem(sys_kal, dict_kal, (0, Tend))
    sol_kal = solve(prob_kal, Tsit5(); dt = dt, adaptive = false, saveat = T_saveat)
    @test sol_kal.retcode == ReturnCode.Success

    ts_compare = T_saveat[1:max(1, div(length(T_saveat), 50)):end]
    x_fw = real(get_solution(sol_fw, x, eqs).(ts_compare))
    x_kal = real(get_solution(sol_kal, x, eqs_det).(ts_compare))
    @test cor(x_fw, x_kal) >= 0.95
    @test maximum(abs.(x_fw .- x_kal)) / maximum(abs.(x_fw)) < 0.2
end
