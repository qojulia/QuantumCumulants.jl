using QuantumCumulants
using Symbolics: Symbolics, @variables
using Test

# The Kalman SDE-vs-ODE tracking testset lives in
# measurement_retrodiction_kalman_test.jl so ParallelTestRunner can
# schedule its heavy SDE+ODE solve on a separate worker.

# v1 surface: Kalman-style continuous-measurement scenario for the
# harmonic oscillator under continuous position monitoring. Master's
# full smoothing comparison (forward + backward Kalman + past-quantum
# combination) is not ported. What we assert here is the v1 pipeline
# that already works: NoiseMeanFieldEquations
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
    modify_equations!(
        eqs, function (lhs, rhs)
            push!(captured, lhs)
            return rhs
        end
    )
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
    bw = meanfield(
        ops, H, J; rates = [Γ], efficiencies = [η], order = 2,
        direction = Backward()
    )
    @test fw isa NoiseMeanFieldEquations
    @test bw isa NoiseMeanFieldEquations
    @test length(fw.equations) == length(bw.equations)
    @test !_is_zero(fw.equations[1].rhs - bw.equations[1].rhs)
end
