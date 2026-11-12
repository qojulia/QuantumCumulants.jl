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
    symtypes = Set(SymbolicUtils.symtype(SymbolicUtils.unwrap(u)) for u in us)
    @test symtypes == Set([Number, Real])
    @test all(u -> is_average(SymbolicUtils.unwrap(u)), us)
end

@testset "structural signature total order" begin
    h = FockSpace(:c)
    a = Destroy(h, :a)
    lt = SQA.canonical_isless
    ops = [a * 1, a' * 1, (a' * a) * 1]
    sorted = sort(ops; lt)
    @test issorted(sorted; lt)                             # comparable / total order
    @test !lt(a * 1, a * 1)                                # irreflexive, deterministic
    for x in ops, y in ops                                 # distinct operators strictly ordered
        isequal(x, y) || @test lt(x, y) != lt(y, x)
    end
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
