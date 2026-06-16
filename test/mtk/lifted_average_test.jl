using QuantumCumulants
using Symbolics: @variables, simplify
using SymbolicUtils: SymbolicUtils
using ModelingToolkitBase: unknowns, System
import SecondQuantizedAlgebra as SQA
import QuantumCumulants: mapleaves, eachleaf
import QuantumCumulants.SecondQuantizedAlgebra: make_time_dependent, is_average
using Test

@testset "lifted average conj is sound (Number symtype)" begin
    h = FockSpace(:c)
    a = Destroy(h, :a)
    @variables t
    u = make_time_dependent(average(a), t)
    @test SymbolicUtils.symtype(u) === Number
    folded = simplify(u - conj(u))
    @test !isequal(folded, 0)            # Real symtype would collapse this to 0
end

@testset "registry: lifted vars, structural identity, conj folding" begin
    h = FockSpace(:c)
    a = Destroy(h, :a)
    @variables Δ κ
    H = Δ * a' * a
    eqs = meanfield([a, a' * a], H, [a]; rates = [κ], order = 2)
    sys = System(eqs; name = :s)
    us = unknowns(sys)
    @test length(us) == 2                                  # ⟨a⟩ and ⟨a†a⟩, conjugates folded
    @test all(u -> SymbolicUtils.symtype(SymbolicUtils.unwrap(u)) === Number, us)
    @test all(u -> is_average(SymbolicUtils.unwrap(u)), us)
end

@testset "structural signature total order" begin
    h = FockSpace(:c)
    a = Destroy(h, :a)
    s1 = SQA.qadd_order_key(a * 1)
    s2 = SQA.qadd_order_key(a' * 1)
    s3 = SQA.qadd_order_key((a' * a) * 1)
    sigs = [s1, s2, s3]
    @test allunique(sigs)                                  # distinct operators -> distinct keys
    @test issorted(sort(sigs))                             # comparable / total order
    @test SQA.qadd_order_key(a * 1) == SQA.qadd_order_key(a * 1)  # deterministic
    @test (s1 < s2) || (s2 < s1)                           # any two are strictly ordered
end

@testset "walkers see lifted averages as leaves" begin
    h = FockSpace(:c)
    a = Destroy(h, :a)
    @variables t
    expr = make_time_dependent(2 * average(a) + average(a' * a), t)
    leaves = eachleaf(expr)
    @test length(leaves) == 2
    @test all(is_average, leaves)
    @test isequal(mapleaves(identity, expr), expr)
end
