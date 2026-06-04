using QCNew
using Test

@testset "find_operators: SQA reexport" begin
    @test isdefined(QCNew, :find_operators)

    hc = FockSpace(:c)
    ops_fock = find_operators(hc, 2)
    # Single Fock mode at order 2: {a, a*a, a'*a}
    @test length(ops_fock) == 3

    ha = NLevelSpace(:a, 2)
    ops_nlevel_1 = find_operators(ha, 1)
    @test length(ops_nlevel_1) == 2

    h = hc ⊗ ha
    ops_prod = find_operators(h, 2)
    @test !isempty(ops_prod)
end
