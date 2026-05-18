# PENDING: port of master test/test_code_quality.jl
#
# Status: superseded by `test/quality/quality_test.jl` (Aqua +
# ExplicitImports + CheckConcreteStructs) and `test/quality/type_stability_test.jl`
# (JET) in the active suite. The master assertions are kept here only as a
# record of the original test thresholds (e.g. `length(get_reports(rep)) <=
# 615`) and the original struct list, which referenced the v0 types
# `MeanfieldEquations`, `IndexedMeanfieldEquations`, `EvaledMeanfieldEquations`,
# `ScaledMeanfieldEquations`, `IndexedMeanfieldNoiseEquations`,
# `MeanfieldNoiseEquations` (all collapsed in v1 into `MeanFieldEquations`
# and `NoiseMeanFieldEquations`).
#
# Verbatim master source (unchanged):

#=
using QuantumCumulants, Test

@testset "best practices" begin
    using Aqua

    Aqua.test_ambiguities([QuantumCumulants]; broken = true)
    Aqua.test_piracies(QuantumCumulants; broken = true)
    Aqua.test_all(
        QuantumCumulants;
        ambiguities = false,
        piracies = false,
        persistent_tasks = false,
    )
end

@testset "ExplicitImports" begin
    using ExplicitImports
    @test check_no_implicit_imports(QuantumCumulants) == nothing
    @test check_all_explicit_imports_via_owners(QuantumCumulants) == nothing
    # @test check_all_explicit_imports_are_public(QuantumCumulants) == nothing
    @test check_no_stale_explicit_imports(QuantumCumulants) == nothing
    @test check_all_qualified_accesses_via_owners(QuantumCumulants) == nothing
    # @test check_all_qualified_accesses_are_public(QuantumCumulants) == nothing
    @test check_no_self_qualified_accesses(QuantumCumulants) == nothing
end

if isempty(VERSION.prerelease)
    @testset "Code linting" begin
        using JET
        # JET.test_package(SecondQuantizedAlgebra; target_defined_modules=true)
        rep = report_package("QuantumCumulants")
        @show rep
        @test length(JET.get_reports(rep)) <= 615 # 317
        @test_broken length(JET.get_reports(rep)) == 0
    end
end

@testset "Concretely typed" begin
    import QuantumCumulants as QC
    using CheckConcreteStructs

    all_concrete(QC.CorrelationFunction)
    all_concrete(QC.Spectrum)

    all_concrete(QC.MeanfieldEquations)
    all_concrete(QC.IndexedMeanfieldEquations)
    all_concrete(QC.EvaledMeanfieldEquations)
    all_concrete(QC.ScaledMeanfieldEquations)

    all_concrete(QC.IndexedMeanfieldNoiseEquations)
    all_concrete(QC.MeanfieldNoiseEquations)
end
=#
