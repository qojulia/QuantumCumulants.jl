using QuantumCumulants
using Symbolics: Symbolics, @variables
using SymbolicUtils: SymbolicUtils
using TermInterface: TermInterface
using Test
import QuantumCumulants as QC

# Regression: a bare `catch` masked real failures in `op(args...)` as a structural rebuild.
@testset "_qc_maketerm: narrow fallback, real errors surface" begin
    @variables x
    xu = SymbolicUtils.unwrap(x)
    T = typeof(xu)

    # `complex(re, im)` rewrites to the factored `re + im*IM` form.
    @test isequal(
        QC._qc_maketerm(T, complex, [xu, 2xu], nothing),
        SymbolicUtils.unwrap(x + 2x * Symbolics.IM),
    )

    # `op` not callable on the rewritten args (here only on `Int`) -> structural rebuild.
    only_int(y::Int) = y + 1
    @test isequal(
        QC._qc_maketerm(T, only_int, [xu], nothing),
        TermInterface.maketerm(T, only_int, [xu], nothing),
    )

    # A genuine error must propagate, not be masked by the fallback.
    @test_throws ErrorException QC._qc_maketerm(T, (_ -> error("boom")), [xu], nothing)
end
