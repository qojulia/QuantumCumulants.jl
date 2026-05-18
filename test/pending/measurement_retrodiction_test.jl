# PENDING: port of master test/test_measurement_retrodiction.jl
#
# Status: blocked on several v1 features —
#   1. `meanfield_backward(...)` is renamed to `meanfield(...; direction =
#      Backward())` (CHANGELOG); the explicit `meanfield_backward` symbol is
#      not exported, the existing tests use the new form.
#   2. `modify_equations(eqs, f)` and `modify_equations!(eqs, f)` for
#      Kalman-filter-style measurement-record terms are not implemented.
#   3. `translate_W_to_Y(eqs)` for converting Wiener-driven SDEs into
#      measurement-record-driven ones is not implemented.
#   4. `@rnumbers` (real-valued cnumber) was replaced by `@variables x::Real`.
#   5. SDESystem code-gen for v1's `NoiseMeanFieldEquations` is built (see
#      examples/heterodyne_detection.jl) but the ensemble-averaging
#      retrodiction examples are flagged as follow-up.
#
# To enable: implement modify_equations + translate_W_to_Y + the SDESystem
# integration path used here.

#=
using QuantumCumulants
using Symbolics: Symbolics, @variables, @syms, @register_symbolic
using SymbolicUtils
using ModelingToolkitBase: @named, mtkcompile, ODEProblem, unknowns
using OrdinaryDiffEq: Euler, solve
using StochasticDiffEq: StochasticDiffEq, SDEProblem, EM
using DiffEqNoiseProcess: NoiseGrid
using Random
using Test

@testset "measurement_retrodiction: harmonic oscillator Kalman smoothing" begin
    h = PhaseSpace(:motion)
    @qnumbers x::Position(h, 1)
    p = Momentum(h, :p, 1)
    @variables Ω::Real Γ::Real η::Real s::Real
    m = 1
    @syms t::Real
    @register_symbolic Ydot(t)
    @register_symbolic Ydot_rev(t)
    a = (x + 1im * p) * s

    H = p^2 / (2m) + 0.5m * Ω^2 * x^2
    J = [a]
    R = [Γ]
    eff = [η]
    ops = [x, p, x * x, x * p, p * p]
    eqs = meanfield(ops, H, J; rates = R, efficiencies = eff, order = 2)

    Ω_, Γ_, η_ = 1, 1 / 6, 0.5
    Tend = 0.1 / Γ_
    dt = Tend / 5e3
    T_saveat = collect(0:dt:Tend)
    ps = [Ω, Γ, η, s]
    pn = [Ω_, Γ_, η_, 1 / √2]
    x0, p0 = 5.0, 0.0
    u0 = [x0, p0, 5 + x0 * x0, im / 2 + x0 * p0, 5 + p0 * p0]

    Random.seed!(1)
    @named sys_fw = to_system(eqs)  # NoiseMeanFieldEquations → SDE-capable System
    dict_fw = merge(Dict(unknowns(sys_fw) .=> u0), Dict(ps .=> pn))
    noise = StochasticDiffEq.RealWienerProcess(0.0, 0.0; save_everystep = true)
    prob_fw = SDEProblem(sys_fw, dict_fw, (0, Tend); noise = noise)
    sol_fw = solve(prob_fw, EM(); dt = dt, saveat = T_saveat)

    # Build measurement record dY from the simulated Wiener path
    t_W = noise.t; W_W = noise.W; N = length(t_W) - 1
    dY_W = zeros(N)
    sol_x = real(get_solution(sol_fw, x, eqs).(t_W[1:N]))
    for k in 1:N
        dW = W_W[k + 1] - W_W[k]
        cd_p_c = sol_x[k] * √(2) * √(Γ_)
        dY_W[k] = dW + √η_ * cd_p_c * dt
    end
    dYdt_data = [dY_W; dY_W[end]] / dt
    dYdt(t) = dYdt_data[Int(round(t / dt)) + 1]
    Ydot(t) = dYdt(t)

    # Kalman-style equations: deterministic + measurement-record drive
    eqs_det = meanfield(ops, H, J; rates = R, order = 2)
    function f_measure(lhs, rhs)
        term_ = √(η * Γ) * (lhs * a + a' * lhs - average(a + a') * lhs) *
                (Ydot(t) - √(η * Γ) * average(a + a'))
        return rhs + cumulant_expansion(average(term_), 2)
    end
    eqs_kal = modify_equations(eqs_det, f_measure)
    @test isequal(eqs_kal[1].lhs, eqs_det[1].lhs)

    # Backward propagation with adjusted recycling + backward measurement
    eqs_back = meanfield(ops[1:(end - 1)], H, J;
                         rates = R, efficiencies = eff, order = 2,
                         direction = Backward())
    eqs_back_c = complete(eqs_back)

    # Translate Wiener-driven SDE to measurement-record-driven SDE
    eqs_c_Y = translate_W_to_Y(eqs)
    @test eqs_c_Y isa NoiseMeanFieldEquations
end
=#
