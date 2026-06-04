using QCNew
using Symbolics: Symbolics, @variables
using Test

@testset "scale! on damped cavity (idempotent)" begin
    hc = FockSpace(:cavity)
    @qnumbers a::Destroy(hc)
    @variables ω κ
    H = ω * a' * a
    eqs = meanfield([a], H, [a]; rates = [κ])
    complete!(eqs)
    n_before = length(eqs.equations)
    scale!(eqs)
    @test length(eqs.equations) == n_before
end

@testset "scale (non-mutating)" begin
    hc = FockSpace(:cavity)
    @qnumbers a::Destroy(hc)
    @variables ω κ
    H = ω * a' * a
    eqs = meanfield([a], H, [a]; rates = [κ])
    complete!(eqs)
    n_before = length(eqs.equations)
    eqs2 = scale(eqs)
    @test length(eqs.equations) == n_before
    @test length(eqs2.equations) == n_before
end
