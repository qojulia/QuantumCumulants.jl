using QuantumCumulants
using Symbolics: Symbolics, @variables
using Test

function _damped_cavity_eqs()
    hc = FockSpace(:cavity)
    a = Destroy(hc, :a)
    @variables ω κ
    H = ω * a' * a
    return meanfield([a], H, [a]; rates = [κ])
end

function _noise_eqs()
    hc = FockSpace(:cavity)
    a = Destroy(hc, :a)
    @variables ω κ η
    H = ω * a' * a
    return meanfield([a], H, [a]; rates = [κ], efficiencies = [η], order = 2)
end

function _correlation_and_spectrum()
    hc = FockSpace(:cavity)
    a = Destroy(hc, :a)
    @variables ω κ
    H = ω * a' * a
    eqs = meanfield([a], H, [a]; rates = [κ])
    eqs_c = complete(eqs; order = 2)
    c = CorrelationFunction(a', a, eqs_c)
    s = Spectrum(c)
    return c, s
end

@testset "show: MeanfieldEquations" begin
    eqs = _damped_cavity_eqs()
    s = repr(eqs)
    @test occursin("∂ₜ(", s)
    @test occursin("⟨a⟩", s)
end

@testset "show: NoiseMeanfieldEquations" begin
    eqs = _noise_eqs()
    s = repr(eqs)
    @test occursin("∂ₜ(", s)
end

@testset "show: CorrelationFunction" begin
    c, _ = _correlation_and_spectrum()
    s = repr(c)
    @test s isa String
    @test !isempty(s)
    @test occursin("a_0", s)
    @test !occursin("⟨a' * a⟩", s)

    eqs_s = repr(c.eqs)
    @test occursin("a_0", eqs_s)
end

@testset "show: Spectrum" begin
    _, spec = _correlation_and_spectrum()
    s = repr(spec)
    @test occursin("ℱ(", s)
    @test occursin(")(ω)", s)
end

@testset "show MIME text/latex: MeanfieldEquations (single row)" begin
    hc = FockSpace(:cavity)
    a = Destroy(hc, :a)
    @variables ω κ
    eqs = meanfield([a], ω * a' * a, [a]; rates = [κ])
    @test repr(MIME("text/latex"), eqs) == "\$\$\n\\begin{aligned}\n" *
        "\\partial_{t} \\langle a \\rangle &= \\langle a \\rangle \\left(  - 0.5 \\kappa - i \\omega \\right)\n" *
        "\\end{aligned}\n\n\$\$"
end

@testset "show MIME text/latex: MeanfieldEquations (row break)" begin
    hc = FockSpace(:cavity)
    a = Destroy(hc, :a)
    @variables ω κ
    # Two rows exercise the tightened `\\[-0.0em]` separator.
    eqs = meanfield([a, a'], ω * a' * a, [a]; rates = [κ])
    @test repr(MIME("text/latex"), eqs) == "\$\$\n\\begin{aligned}\n" *
        "\\partial_{t} \\langle a \\rangle &= \\langle a \\rangle \\left(  - 0.5 \\kappa - i \\omega \\right) \\\\[-0.0em]\n" *
        "\\partial_{t} \\langle a^{\\dagger} \\rangle &= \\langle a^{\\dagger} \\rangle \\left(  - 0.5 \\kappa + i \\omega \\right)\n" *
        "\\end{aligned}\n\n\$\$"
end

@testset "show MIME text/latex: NoiseMeanfieldEquations" begin
    hc = FockSpace(:cavity)
    a = Destroy(hc, :a)
    @variables ω κ η
    # The noise term renders as a factored `dW/dt`, not `\frac{coeff·dW}{dt}`.
    eqs = meanfield([a], ω * a' * a, [a]; rates = [κ], efficiencies = [η], order = 2)
    @test repr(MIME("text/latex"), eqs) == "\$\$\n\\begin{aligned}\n" *
        "\\partial_{t} \\langle a \\rangle &= \\langle a \\rangle \\left(  - 0.5 \\kappa - i \\omega \\right) + " *
        "\\frac{\\mathrm{d}W}{\\mathrm{d}t} \\left( \\langle aa \\rangle \\sqrt{\\eta \\kappa} + " *
        "\\langle a^{\\dagger}a \\rangle \\sqrt{\\eta \\kappa} + " *
        "\\langle a \\rangle \\left(  - \\langle a \\rangle - \\langle a^{\\dagger} \\rangle \\right) \\sqrt{\\eta \\kappa} \\right)\n" *
        "\\end{aligned}\n\n\$\$"
end

@testset "show MIME text/latex: CorrelationFunction" begin
    c, _ = _correlation_and_spectrum()
    s = repr(MIME("text/latex"), c)
    @test !isempty(s)
    # The renamed ancilla `a_0` renders as a subscripted `a_{\mathrm{0}}`.
    @test occursin("a_{\\mathrm{0}}", s)
end

@testset "show MIME text/latex: Spectrum" begin
    _, spec = _correlation_and_spectrum()
    s = repr(MIME("text/latex"), spec)
    @test occursin("mathcal{F}", s)
end
