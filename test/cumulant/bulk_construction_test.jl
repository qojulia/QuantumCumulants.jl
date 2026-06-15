using QuantumCumulants
using SymbolicUtils: SymbolicUtils
using Symbolics: Symbolics, expand
using Test

# Robust zero-check, copied verbatim from test/cumulant/cumulant_test.jl (test files
# run in separate workers, so the helper is not shared and must be redefined here).
_iz(x) =
    x isa Number ? iszero(x) :
    x isa SymbolicUtils.BasicSymbolic ? SymbolicUtils._iszero(expand(x)) : iszero(x)

# Length-N product over N distinct Fock modes (order N), as an average leaf.
function _prodavg(N)
    h = reduce(⊗, (FockSpace(Symbol(:f, i)) for i in 1:N))
    return average(reduce(*, [Destroy(h, Symbol(:a, j), j) for j in 1:N]))
end

@testset "bulk-construction: output is unchanged (golden)" begin
    h = FockSpace(:f)
    a = Destroy(h, :a)
    # no-op identity
    @test isequal(cumulant_expansion(average(a' * a), 2), average(a' * a))
    # order-1 factorization of three distinct modes
    h3 = reduce(⊗, (FockSpace(Symbol(:f, i)) for i in 1:3))
    a1 = Destroy(h3, :a1, 1)
    a2 = Destroy(h3, :a2, 2)
    a3 = Destroy(h3, :a3, 3)
    @test _iz(
        cumulant_expansion(average(a1 * a2 * a3), 1) -
            average(a1) * average(a2) * average(a3)
    )
    # Round-trip stability: expanding an already-expanded result to the same order is a fixed point.
    for N in 3:5, ord in 2:(N - 1)
        e1 = cumulant_expansion(_prodavg(N), ord)
        @test _iz(cumulant_expansion(e1, ord) - e1)
    end
end

@testset "bulk-construction: allocation budget" begin
    avg = _prodavg(7)
    cumulant_expansion(avg, 2)                       # compile
    # Regression guard for the memoised expansion. A length-7 order-2 product reaches only
    # 99 distinct sub-blocks and 28 distinct moment leaves, but a non-memoised expansion
    # reconstructs them ~98k times and allocates ~250 MB here; the memoised build allocates
    # ~9 MB. 30 MB sits well between, so losing the memo (or reverting to incremental
    # accumulation) trips this while the memoised build stays comfortably under.
    @test (@allocated cumulant_expansion(avg, 2)) < 30_000_000
end

@testset "bulk-construction: no-op returns the identical object" begin
    h = FockSpace(:f)
    a = Destroy(h, :a)
    x = SymbolicUtils.unwrap(average(a' * a) + average(a))   # max order 2
    @test cumulant_expansion(x, 2) === x                      # nothing exceeds order 2
    @test cumulant_expansion(x, [2, 2]) === x
end

@testset "bulk-construction: re-expanding to the same order returns the same eqs" begin
    h = FockSpace(:f)
    a = Destroy(h, :a)
    eqs = meanfield([a' * a], a' * a, [a]; rates = [1.0])
    e2 = cumulant_expansion(eqs, 2)
    @test cumulant_expansion(e2, 2) === e2          # idempotent: already at order [2]
end
