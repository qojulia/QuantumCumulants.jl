using QuantumCumulants, Test

@testset "Code linting" begin
    using JET
    JET.test_package(QuantumCumulants; target_defined_modules=true)
end

@testset "Code quality" begin
    using ExplicitImports, Aqua
    @test check_no_stale_explicit_imports(QuantumCumulants) == nothing
    @test check_all_explicit_imports_via_owners(QuantumCumulants) == nothing
    Aqua.test_ambiguities([QuantumCumulants])
    Aqua.test_all(QuantumCumulants;ambiguities=false)
end
