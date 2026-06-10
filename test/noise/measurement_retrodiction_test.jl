using QuantumCumulants
using Symbolics: Symbolics, @variables
using Test

@testset "meanfield: direction = Backward()" begin
    h = FockSpace(:cavity)
    @qnumbers a::Destroy(h)
    @variables ω κ
    H = ω * a' * a
    fw = meanfield([a], H, [a]; rates = [κ])
    bw = meanfield([a], H, [a]; rates = [κ], direction = Backward())
    @test length(fw.equations) == length(bw.equations) == 1
    rhs_fw = fw.equations[1].rhs
    rhs_bw = bw.equations[1].rhs
    @test !isequal(rhs_fw, rhs_bw)
end

@testset "Backward: deterministic drift includes Lindblad recycling + trace" begin
    # For H = ω·a†a the commutator [H, a†a] vanishes, so a commutator-only backward RHS
    # (the bug where a misplaced `return` dropped the recycling and trace terms) would be
    # exactly 0. The correct backward retrodiction drift of ⟨a†a⟩ under jump a is
    # κ(1 + ⟨a†a⟩); this pins the recycling + trace contributions, which `!isequal(fw, bw)`
    # alone does not (the commutator sign flip makes fw ≠ bw even with the terms dropped).
    h = FockSpace(:cavity)
    @qnumbers a::Destroy(h)
    @variables ω κ
    H = ω * a' * a
    bw = meanfield([a' * a], H, [a]; rates = [κ], direction = Backward())
    rhs = bw.equations[1].rhs
    @test !_is_zero(rhs)                              # would be 0 under the dropped-terms bug
    @test _is_zero(rhs - (κ + κ * average(a' * a)))   # exact backward drift κ(1 + ⟨a†a⟩)
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
    @test eqs isa NoiseMeanfieldEquations
    out = translate_W_to_Y(eqs)
    @test out isa NoiseMeanfieldEquations
    @test length(out.equations) == length(eqs.equations)
    # The dY-form adds a deterministic correction to the drift; the noise terms
    # are unchanged.
    @test !_is_zero(out.equations[1].rhs - eqs.equations[1].rhs)
    @test _is_zero(out.noise_equations[1].rhs - eqs.noise_equations[1].rhs)
end

@testset "Backward NoiseMeanfieldEquations: drift sign flip, noise unchanged" begin
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
    @test fw isa NoiseMeanfieldEquations
    @test bw isa NoiseMeanfieldEquations
    @test length(fw.equations) == length(bw.equations)
    @test !_is_zero(fw.equations[1].rhs - bw.equations[1].rhs)
end
