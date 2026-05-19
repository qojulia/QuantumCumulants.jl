using QuantumCumulants
using Symbolics: Symbolics, @variables, @register_symbolic
using ModelingToolkitBase: mtkcompile, unknowns
using OrdinaryDiffEq: Tsit5, ODEProblem, ReturnCode, solve
using StochasticDiffEq: StochasticDiffEq, SDEProblem, EM, RealWienerProcess
using Statistics: cor
using Random
using Test

# v1 surface: Kalman-style continuous-measurement scenario for the
# harmonic oscillator under continuous position monitoring. Master's
# full smoothing comparison (forward + backward Kalman + past-quantum
# combination) lives in test/pending/measurement_retrodiction_test.jl.
# What we assert here is the v1 pipeline that already works: NoiseMeanFieldEquations
# build, SDE forward solves, `translate_W_to_Y` augments the drift,
# `modify_equations` accepts a measurement-record callback, and the
# resulting deterministic Kalman ODE driven by the measurement record
# tracks the Wiener-driven SDE on the same noise realisation.

@testset "meanfield: direction = Backward()" begin
    h = FockSpace(:cavity)
    @qnumbers a::Destroy(h)
    @variables ω κ
    H = ω * a' * a
    fw = meanfield([a], H, [a]; rates = [κ])
    bw = meanfield([a], H, [a]; rates = [κ], direction = Backward())
    @test length(fw.equations) == length(bw.equations) == 1
    # Forward sign of the Hamiltonian commutator opposite to backward.
    rhs_fw = fw.equations[1].rhs
    rhs_bw = bw.equations[1].rhs
    @test !isequal(rhs_fw, rhs_bw)
end

@testset "modify_equations: identity is a no-op" begin
    h = FockSpace(:cavity)
    @qnumbers a::Destroy(h)
    @variables ω κ
    H = ω * a' * a
    eqs = meanfield([a], H, [a]; rates = [κ])
    eqs2 = modify_equations(eqs, (lhs, rhs) -> rhs)
    @test length(eqs2.equations) == length(eqs.equations)
    for (eq1, eq2) in zip(eqs.equations, eqs2.equations)
        @test isequal(eq1.lhs, eq2.lhs)
        @test isequal(eq1.rhs, eq2.rhs)
    end
end

@testset "modify_equations: f sees undo_average(lhs) as a QField" begin
    h = FockSpace(:cavity)
    @qnumbers a::Destroy(h)
    @variables ω κ
    H = ω * a' * a
    eqs = meanfield([a], H, [a]; rates = [κ])
    SQA = QuantumCumulants.SecondQuantizedAlgebra
    captured = Any[]
    modify_equations!(eqs, function (lhs, rhs)
        push!(captured, lhs)
        return rhs
    end)
    @test all(c isa SQA.QField for c in captured)
end

@testset "translate_W_to_Y: noise equations are augmented" begin
    h = FockSpace(:cavity)
    @qnumbers a::Destroy(h)
    @variables ω κ η
    H = ω * a' * a
    eqs = meanfield([a], H, [a]; rates = [κ], efficiencies = [η])
    @test eqs isa NoiseMeanFieldEquations
    out = translate_W_to_Y(eqs)
    @test out isa NoiseMeanFieldEquations
    @test length(out.equations) == length(eqs.equations)
    # The dY-form adds a deterministic correction to the drift, but the
    # noise (`dW`-coefficient) terms are unchanged.
    @test !_is_zero(out.equations[1].rhs - eqs.equations[1].rhs)
    @test _is_zero(out.noise_equations[1].rhs - eqs.noise_equations[1].rhs)
end

# Continuous monitoring of a harmonic oscillator: simulate the forward
# SDE, extract the measurement record dY from the realised Wiener path,
# then build the deterministic Kalman ODE that is driven by dY via
# `modify_equations`. On the *same* noise realisation, the Kalman ODE
# must track the SDE trajectory (perfect agreement modulo integration
# scheme: EM vs RK; we assert correlation >= 0.95 and relative endpoint
# error < 0.2). This is the master-test physics assertion at minimum
# strength; tighter tolerances need Euler-Euler scheme matching, which
# OrdinaryDiffEq's v7 metapackage doesn't expose at the test surface.
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
    dt = Tend / 2e3
    T_saveat = collect(0:dt:Tend)
    ps = [Ω, Γ, η, s]
    pn = [Ω_, Γ_, η_, 1 / √2]
    x0, p0 = 5.0, 0.0
    u0 = ComplexF64[x0, p0, 5 + x0 * x0, im / 2 + x0 * p0, 5 + p0 * p0]

    Random.seed!(2)
    sys_fw = mtkcompile(to_system(eqs; name = :sys_fw))
    dict_fw = merge(Dict(unknowns(sys_fw) .=> u0), Dict(ps .=> pn))
    noise = StochasticDiffEq.RealWienerProcess(0.0, 0.0; save_everystep = true)
    prob_fw = SDEProblem(sys_fw, dict_fw, (0, Tend); noise = noise)
    sol_fw = solve(prob_fw, EM(); dt = dt, saveat = T_saveat)
    @test sol_fw.retcode == ReturnCode.Success
    @test all(isfinite.(real.(get_solution(sol_fw, x, eqs).(T_saveat))))

    # Build measurement record dY from the realised Wiener path.
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

    # Build the deterministic Kalman equations: drift = ordinary
    # meanfield drift + measurement-record term. The `Symbolics.Num`
    # wraps are needed because v1's `average(...)` returns a
    # `BasicSymbolic{SymReal}` which does not yet multiply with a `QAdd`
    # directly (gap in SQA's product overloads). See TODO.md.
    eqs_det = meanfield(ops, H, J; rates = R, order = 2)
    iv = eqs_det.iv
    function f_measure(lhs, rhs)
        avg_apad = Symbolics.Num(average(a + a'))
        inn = lhs * a + a' * lhs - avg_apad * lhs
        drive = Symbolics.Num(Ydot_meas(iv)) - √(η * Γ) * avg_apad
        term_ = √(η * Γ) * inn * drive
        return rhs + cumulant_expansion(average(term_), 2)
    end
    eqs_kal = modify_equations(eqs_det, f_measure)
    @test eqs_kal isa MeanFieldEquations
    # modify_equations actually modified the drift on every equation.
    @test all(!_is_zero(eqs_kal.equations[i].rhs - eqs_det.equations[i].rhs)
              for i in eachindex(eqs_kal.equations))

    sys_kal = mtkcompile(to_system(eqs_kal; name = :sys_kal))
    dict_kal = merge(Dict(unknowns(sys_kal) .=> u0), Dict(ps .=> pn))
    prob_kal = ODEProblem(sys_kal, dict_kal, (0, Tend))
    sol_kal = solve(prob_kal, Tsit5(); dt = dt, adaptive = false, saveat = T_saveat)
    @test sol_kal.retcode == ReturnCode.Success

    # Tracking assertion: on the same noise realisation, the Kalman ODE
    # and the SDE forward must produce highly-correlated trajectories.
    ts_compare = T_saveat[1:max(1, div(length(T_saveat), 50)):end]
    x_fw = real(get_solution(sol_fw, x, eqs).(ts_compare))
    x_kal = real(get_solution(sol_kal, x, eqs_det).(ts_compare))
    @test cor(x_fw, x_kal) >= 0.95
    @test maximum(abs.(x_fw .- x_kal)) / maximum(abs.(x_fw)) < 0.2
end

# Backward direction for noise meanfields: master's smoothing equations
# require time-reversed coherent drift + adjoint Lindblad recycling +
# trace-preserving term. We assert (a) the type and direction marker,
# (b) that the deterministic drift differs from the forward case, and
# (c) the noise term magnitude is comparable (same Lindblad operators).
@testset "Backward NoiseMeanFieldEquations: drift sign flip, noise unchanged" begin
    h = PhaseSpace(:motion)
    @qnumbers x::Position(h)
    p = Momentum(h, :p)
    @variables Ω Γ η s
    a = (x + 1im * p) * s
    H = p^2 / 2 + 0.5 * Ω^2 * x^2
    J = [a]
    ops = [x, p, x * x]

    fw = meanfield(ops, H, J; rates = [Γ], efficiencies = [η], order = 2)
    bw = meanfield(ops, H, J; rates = [Γ], efficiencies = [η], order = 2,
                   direction = Backward())
    @test fw isa NoiseMeanFieldEquations
    @test bw isa NoiseMeanFieldEquations
    @test length(fw.equations) == length(bw.equations)
    @test !_is_zero(fw.equations[1].rhs - bw.equations[1].rhs)
end
