using QuantumCumulants
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
    using QuantumCumulants
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

    """
        phase_invariant(x) -> Bool

    U(1) phase-invariance filter for laser-type models. A moment is kept when its
    net phase vanishes: an annihilator contributes -1, a creator +1, and a
    transition `σ_{l1,l2}` contributes `l1 - l2`. Pass as `filter_func` to
    `complete`/`complete!`/`CorrelationFunction` to drop the fast-rotating moments.
    """
    phase_invariant(x) = iszero(_net_phase(x))
    _net_phase(_) = 0
    function _net_phase(op::Op)
        is_destroy(op) && return -1
        is_create(op) && return 1
        return is_transition(op) ? op.l1 - op.l2 : 0
    end
    _net_phase(t::QuantumCumulants.SecondQuantizedAlgebra.QTerm) =
        isempty(t.ops) ? 0 : sum(_net_phase(op) for op in t.ops)
    function _net_phase(q::QuantumCumulants.SecondQuantizedAlgebra.QAdd)
        isempty(q.arguments) && return 0
        term, _ = first(q.arguments)
        return _net_phase(term)
    end
    function _net_phase(avg::SymbolicUtils.BasicSymbolic)
        is_average(avg) || return 0
        return _net_phase(undo_average(avg))
    end
end

ParallelTestRunner.runtests(QuantumCumulants, args; testsuite, init_code)
