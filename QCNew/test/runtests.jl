using QCNew
using ParallelTestRunner: ParallelTestRunner

# Start with autodiscovered tests
testsuite = ParallelTestRunner.find_tests(@__DIR__)

# Parse arguments
args = ParallelTestRunner.parse_args(ARGS)

if ParallelTestRunner.filter_tests!(testsuite, args)
    delete!(testsuite, "quality/JET")
end

# Shared helpers available to every test file. Defined here (not in a
# separate test/common.jl) because ParallelTestRunner discovers every *.jl
# file as a test; using `init_code` is the documented pattern for sharing
# code across test workers without polluting the test list.
init_code = quote
    using QCNew
    using SymbolicUtils: SymbolicUtils
    using Symbolics: Symbolics
    using Test

    """
        assert_real(sol, op, eqs; atol=1e-6, ts=sol.t)

    Assert observable trajectory is real (within tolerance). For Hermitian
    observables (number ops, populations), the average must be real along
    the trajectory.
    """
    function assert_real(sol, op, eqs; atol = 1.0e-6, ts = sol.t)
        traj = [get_solution(sol, op, eqs)(t) for t in ts]
        return @test maximum(abs.(imag.(traj))) < atol
    end

    """
        assert_bounded(sol, op, eqs, lo, hi; atol=1e-6, ts=sol.t)

    Assert observable trajectory stays within [lo, hi].
    """
    function assert_bounded(sol, op, eqs, lo, hi; atol = 1.0e-6, ts = sol.t)
        traj = [real(get_solution(sol, op, eqs)(t)) for t in ts]
        @test minimum(traj) >= lo - atol
        return @test maximum(traj) <= hi + atol
    end

    """
        assert_nonneg(sol, op, eqs; atol=1e-6, ts=sol.t)

    Assert observable trajectory stays non-negative. Use for photon number,
    populations, and other physical quantities required to be ≥ 0.
    """
    assert_nonneg(sol, op, eqs; atol = 1.0e-6, ts = sol.t) =
        assert_bounded(sol, op, eqs, 0.0, Inf; atol, ts)

    """
        assert_population(sol, op, eqs; atol=1e-6, ts=sol.t)

    Assert observable is a valid population: trajectory ⊂ [0, 1].
    """
    assert_population(sol, op, eqs; atol = 1.0e-6, ts = sol.t) =
        assert_bounded(sol, op, eqs, 0.0, 1.0; atol, ts)

    """
        assert_steady(sol; atol=1e-4, k_last=2)

    Assert the trajectory reached a steady state at the END: the numerical
    time-derivative at sol.t[end] is below atol componentwise.
    """
    function assert_steady(sol; atol = 1.0e-4, k_last::Int = 2)
        n = length(sol.t)
        n >= k_last || return
        t_end = sol.t[end]
        t_prev = sol.t[end - k_last + 1]
        dt = t_end - t_prev
        dt > 0 || return
        for k in eachindex(sol.u[end])
            du = abs(sol.u[end][k] - sol.u[end - k_last + 1][k]) / dt
            @test du < atol
        end
        return
    end

    """
        _is_zero(x)

    Decide whether a symbolic / numeric expression is zero. Use for symbolic
    RHS comparisons where `iszero` alone is too coarse due to representational
    quirks.
    """
    _is_zero(x::Number) = iszero(x)
    function _is_zero(x)
        s = try
            Symbolics.simplify(x; expand = true)
        catch
            x
        end
        s isa Number && return iszero(s)
        if s isa SymbolicUtils.BasicSymbolic
            SymbolicUtils.isconst(s) && return iszero(s.val)
            return isequal(s, 0)
        end
        return isequal(s, 0)
    end
end

ParallelTestRunner.runtests(QCNew, args; testsuite, init_code)
