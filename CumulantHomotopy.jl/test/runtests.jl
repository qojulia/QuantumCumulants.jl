using Test
using CumulantHomotopy

@testset "CumulantHomotopy" begin
    @test isdefined(CumulantHomotopy, :stationary_state)
    @test isdefined(CumulantHomotopy, :stationary_sequence)
end
