using QuantumCumulants
using Test

@testset "find_operators: SQA reexport" begin
    @test isdefined(QuantumCumulants, :find_operators)

    hc = FockSpace(:c)
    ops_fock = find_operators(hc, 2)
    # Single Fock mode at order 2: {a, a*a, a'*a}
    @test length(ops_fock) == 3

    ha = NLevelSpace(:a, 2)
    ops_nlevel_1 = find_operators(ha, 1)
    @test length(ops_nlevel_1) == 2

    h = hc ⊗ ha
    ops_prod = find_operators(h, 2)
    # Same coverage works on a product space; closure under the cumulant
    # expansion is checked elsewhere. Here we just sanity-check the call.
    @test !isempty(ops_prod)
end
