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

@testset "show MIME text/latex: MeanfieldEquations" begin
    eqs = _damped_cavity_eqs()
    s = repr(MIME("text/latex"), eqs)
    # Equations render as a display-math `aligned` environment so Documenter/Markdown does
    # not mangle the subscripts (see `_latex_display`).
    @test occursin("begin{aligned}", s)
    @test startswith(s, "\$\$") && endswith(rstrip(s), "\$\$")
end

@testset "show MIME text/latex: NoiseMeanfieldEquations" begin
    eqs = _noise_eqs()
    s = repr(MIME("text/latex"), eqs)
    @test !isempty(s)
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
