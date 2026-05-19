using QuantumCumulants
using Symbolics: Symbolics, @variables
using SymbolicUtils
using ModelingToolkitBase: mtkcompile
using OrdinaryDiffEq: ODEProblem, Tsit5, solve
using Test

# v1 surface: indexed `CorrelationFunction` (unified API, no separate
# `IndexedCorrelationFunction` type). Master's full test exercises
# `IndexedCorrelationFunction`, `evaluate(corr, 1, 2; limits=...)`,
# `scale(corr)`, `split_sums`, and several SQA helpers (`DoubleNumberedVariable`,
# `SingleNumberedVariable`, `value_map`) that don't exist in v1; those parts
# stay in test/pending/indexed_correlation_test.jl until the order=2 indexed
# JC closure speeds up (the same `complete!` mixed-order issue blocks
# indexed_mixed_order_test).

# Master also exercises an order=2 JC closure with a `phase_invariant`
# filter and `split_sums` to chop the symmetric atom block; in v1 the
# mixed-order `complete!` does not close that system within reasonable
# iteration bounds (see TODO.md "Open v1 feature gaps"). The order=1
# tests below exercise the phase filter on the available surface.

const _SQA = QuantumCumulants.SecondQuantizedAlgebra

# Phase-invariant filter for the JC laser: counts excess of creators
# over annihilators (cavity) and σ_{2,1} over σ_{1,2} (atom) per term.
# Same construction as in master's `test_indexed_correlation.jl` and
# `higher_order_test.jl`.
_ϕ(::_SQA.Destroy) = -1
_ϕ(::_SQA.Create) = 1
function _ϕ(t::_SQA.Transition)
    t.i == t.j && return 0
    return t.i == 2 ? 1 : -1
end
function _ϕ(q::_SQA.QAdd)
    isempty(q.arguments) && return 0
    qt, _ = first(q.arguments)
    return sum(_ϕ(op) for op in qt.ops; init = 0)
end
function _ϕ(avg)
    avg isa SymbolicUtils.BasicSymbolic || return 0
    _SQA.is_average(avg) || return 0
    return _ϕ(_SQA.undo_average(avg))
end
_phase_invariant(x) = iszero(_ϕ(x))

@testset "indexed CorrelationFunction: JC laser, order=1" begin
    @variables N::Real Δ::Real g::Real κ::Real Γ::Real R::Real ν::Real
    hc = FockSpace(:cavity)
    ha = NLevelSpace(:atom, 2)
    h = hc ⊗ ha
    k = Index(h, :k, N, ha)
    @qnumbers a::Destroy(h)
    σ(i, j, idx) = IndexedOperator(Transition(h, :σ, i, j), idx)

    H = -Δ * a' * a + g * (Σ(a' * σ(1, 2, k), k) + Σ(a * σ(2, 1, k), k))
    J = [a, σ(1, 2, k), σ(2, 1, k), σ(2, 2, k)]
    rates = [κ, Γ, R, ν]
    ops = [a' * a, a]

    eqs = meanfield(ops, H, J; rates = rates, order = 1)
    eqs_c = complete(eqs)
    @test length(eqs_c.equations) >= 1

    # Scale should produce a finite-size closed system.
    eqs_sc = scale(eqs_c)
    @test length(eqs_sc.equations) >= 1
    @test isempty(find_missing(eqs_sc; get_adjoints = false))
end

@testset "indexed CorrelationFunction: phase-invariant filter trims LHS" begin
    @variables N::Real Δ::Real g::Real κ::Real Γ::Real R::Real ν::Real
    hc = FockSpace(:cavity)
    ha = NLevelSpace(:atom, 2)
    h = hc ⊗ ha
    k = Index(h, :k, N, ha)
    @qnumbers a::Destroy(h)
    σ(i, j, idx) = IndexedOperator(Transition(h, :σ, i, j), idx)

    H = -Δ * a' * a + g * (Σ(a' * σ(1, 2, k), k) + Σ(a * σ(2, 1, k), k))
    J = [a, σ(1, 2, k), σ(2, 1, k), σ(2, 2, k)]
    rates = [κ, Γ, R, ν]
    ops = [a' * a, σ(2, 2, k)]   # both phase-invariant

    eqs = meanfield(ops, H, J; rates = rates, order = 1)
    eqs_c_nofilt = complete(eqs)
    eqs_c_filt = complete(eqs; filter_func = _phase_invariant)

    # The full closure picks up the coherences ⟨a⟩, ⟨a'⟩, ⟨σ12⟩, ⟨σ21⟩;
    # the phase-invariant filter keeps only ⟨a'a⟩ and ⟨σ22⟩.
    @test length(eqs_c_nofilt.equations) > length(eqs_c_filt.equations)
    @test length(eqs_c_filt.equations) == 2
    @test isempty(find_missing(eqs_c_filt; filter_func = _phase_invariant))
end

@testset "indexed CorrelationFunction: order=1 laser steady state via evaluate(N=>1)" begin
    # Numerical correctness check at order=1 (rate equations) for a single
    # atom: at order=1 mean-field has no cavity coherence, so ⟨a'a⟩ → 0
    # in steady state, and ⟨σ22⟩ → R/(R+Γ) is the bare two-level pumping
    # steady state. We use `evaluate(eqs_c; limits = (N => 1))` rather than
    # `scale(eqs_c)`; the v1 `scale` path puts an extra `N` factor on the
    # per-atom decay coefficient (see TODO.md "Open v1 feature gaps") that
    # masks this analytic check. `evaluate(N=>1)` collapses the bound atom
    # index to a single concrete atom and reproduces master physics.
    @variables N::Real Δ::Real g::Real κ::Real Γ::Real R::Real ν::Real
    hc = FockSpace(:cavity)
    ha = NLevelSpace(:atom, 2)
    h = hc ⊗ ha
    k = Index(h, :k, N, ha)
    @qnumbers a::Destroy(h)
    σ(i, j, idx) = IndexedOperator(Transition(h, :σ, i, j), idx)

    H = -Δ * a' * a + g * (Σ(a' * σ(1, 2, k), k) + Σ(a * σ(2, 1, k), k))
    J = [a, σ(1, 2, k), σ(2, 1, k), σ(2, 2, k)]
    rates = [κ, Γ, R, ν]
    ops = [a' * a, σ(2, 2, k)]

    eqs = meanfield(ops, H, J; rates = rates, order = 1)
    eqs_c = complete(eqs)
    eqs_ev = evaluate(eqs_c; limits = (N => 1))
    sys = mtkcompile(to_system(eqs_ev; name = :sys))

    Δ_, g_, κ_, Γ_, R_, ν_ = 0.0, 1.0, 1.0, 0.25, 4.0, 1.0
    p_map = parameter_map(
        eqs_ev, [Δ => Δ_, g => g_, κ => κ_, Γ => Γ_, R => R_, ν => ν_]
    )
    u0 = initial_values(eqs_ev, ComplexF64[0.0 for _ in eqs_ev.states])
    prob = ODEProblem(sys, merge(u0, p_map), (0.0, 50.0))
    sol = solve(prob, Tsit5(); reltol = 1e-8, abstol = 1e-10, saveat = 50.0)

    n_ss = real(get_solution(sol, a' * a, eqs_ev)(50.0))
    @test isapprox(n_ss, 0.0; atol = 1e-6)

    # The atom-state observable in eqs_ev is indexed at the concrete atom
    # `k_1` that `evaluate(N=>1)` minted; pull it from the LHS list rather
    # than reconstructing.
    σ22_lhs = first(eq.lhs for eq in eqs_ev.equations
                   if startswith(string(eq.lhs), "⟨σ"))
    σ22_op = _SQA.undo_average(σ22_lhs)
    σ22_ss = real(get_solution(sol, σ22_op, eqs_ev)(50.0))
    @test isapprox(σ22_ss, R_ / (R_ + Γ_); atol = 1e-3)
    @test 0.0 <= σ22_ss <= 1.0
end

@testset "indexed CorrelationFunction: g^(1)(τ) construction" begin
    @variables N::Real Δ::Real g::Real κ::Real Γ::Real R::Real ν::Real
    hc = FockSpace(:cavity)
    ha = NLevelSpace(:atom, 2)
    h = hc ⊗ ha
    k = Index(h, :k, N, ha)
    @qnumbers a::Destroy(h)
    σ(i, j, idx) = IndexedOperator(Transition(h, :σ, i, j), idx)

    H = -Δ * a' * a + g * (Σ(a' * σ(1, 2, k), k) + Σ(a * σ(2, 1, k), k))
    J = [a, σ(1, 2, k)]
    rates = [κ, Γ]
    eqs = meanfield([a' * a, a], H, J; rates = rates, order = 1)
    eqs_c = complete(eqs)

    # First-order correlation of the cavity mode.
    corr = CorrelationFunction(a', a, eqs_c)
    @test corr isa CorrelationFunction
    @test length(corr.eqs.equations) >= 1

    # Scale on the correlation function: same surface as on regular equations.
    corr_sc = scale(corr)
    @test corr_sc isa CorrelationFunction
end
