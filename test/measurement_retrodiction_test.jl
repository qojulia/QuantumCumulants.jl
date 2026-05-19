using QuantumCumulants
using Symbolics: Symbolics, @variables
using Test

# v1 surface only. Master's full Kalman-smoothing scenario lives in
# test/pending/measurement_retrodiction_test.jl and requires
# Wiener-process integration + measurement-record code we'll re-port once
# the SDE bridge stabilises (CHANGELOG follow-up).

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
end
