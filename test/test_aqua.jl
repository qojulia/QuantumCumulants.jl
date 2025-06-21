using QuantumCumulants, Test

@testset "Code Quality" begin
    using Aqua
    Aqua.test_piracies(QuantumCumulants, broken = true)
    Aqua.test_undefined_exports(QuantumCumulants)
    Aqua.test_project_extras(QuantumCumulants)
    Aqua.test_stale_deps(QuantumCumulants)
end
