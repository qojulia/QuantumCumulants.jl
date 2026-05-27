using QuantumCumulants
using Symbolics: Symbolics, @variables
using SymbolicUtils
using ModelingToolkitBase: mtkcompile
using OrdinaryDiffEq: ODEProblem, Tsit5, solve
using Test

# v1 surface: indexed `CorrelationFunction` (unified API, no separate
# `IndexedCorrelationFunction` type). Master's full test additionally
# exercises `split_sums` and SQA helpers (`DoubleNumberedVariable`,
# `SingleNumberedVariable`, `value_map`) that are intentionally absent
# in v1. The order=1 tests below exercise the v1 phase-filter surface.

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
    sys = mtkcompile(System(eqs_ev; name = :sys))

    Δ_, g_, κ_, Γ_, R_, ν_ = 0.0, 1.0, 1.0, 0.25, 4.0, 1.0
    p_map = parameter_map(
        eqs_ev, [Δ => Δ_, g => g_, κ => κ_, Γ => Γ_, R => R_, ν => ν_]
    )
    u0 = initial_values(eqs_ev, ComplexF64[0.0 for _ in eqs_ev.states])
    prob = ODEProblem(sys, merge(u0, p_map), (0.0, 50.0))
    sol = solve(prob, Tsit5(); reltol = 1.0e-8, abstol = 1.0e-10, saveat = 50.0)

    n_ss = real(get_solution(sol, a' * a, eqs_ev)(50.0))
    @test isapprox(n_ss, 0.0; atol = 1.0e-6)

    # The atom-state observable in eqs_ev is indexed at the concrete atom
    # `k_1` that `evaluate(N=>1)` minted; pull it from the LHS list rather
    # than reconstructing.
    σ22_lhs = first(
        eq.lhs for eq in eqs_ev.equations
            if startswith(string(eq.lhs), "⟨σ")
    )
    σ22_op = _SQA.undo_average(σ22_lhs)
    σ22_ss = real(get_solution(sol, σ22_op, eqs_ev)(50.0))
    @test isapprox(σ22_ss, R_ / (R_ + Γ_); atol = 1.0e-3)
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


@testset "indexed CorrelationFunction: ancilla canon threads parent canon" begin
    # Regression for `_complete_ancilla!`'s `parent_canon` kwarg
    # (src/correlation.jl). After `scale`, the parent system's
    # canonical-index pool on the atom subspace is e.g. `[j_2_1, j_2_2]`
    # (materialised via the fallback `Index(:j)(2)` because user `j` is
    # bound by H). When the ancilla derivation for the correlation
    # function builds its own canon from scratch, it must reuse the
    # parent's vocabulary so that steady-state lookups against
    # `eqs0.states` succeed. Pre-fix the ancilla would mint fresh
    # `j_3`-suffixed slots, the lookup missed, and `Spectrum` collapsed
    # to a constant (superradiant_laser_indexed spectrum bug).
    @variables Δ::Real g::Real κ::Real Γ::Real R::Real ν::Real N::Real
    hc = FockSpace(:cavity); ha = NLevelSpace(:atom, 2); h = hc ⊗ ha
    @qnumbers a::Destroy(h, 1)
    σ(α, β, k) = IndexedOperator(Transition(h, :σ, α, β, 2), k)
    j_ind = Index(h, :j, N, ha)
    H = Δ * Σ(σ(2, 2, j_ind), j_ind) +
        g * Σ(a' * σ(1, 2, j_ind) + a * σ(2, 1, j_ind), j_ind)
    J = [a, σ(1, 2, j_ind), σ(2, 1, j_ind), σ(2, 2, j_ind)]
    eqs_c = complete(meanfield([a' * a], H, J; rates = [κ, Γ, R, ν], order = 2))
    eqs_sc = scale(eqs_c)

    parent_canon = QuantumCumulants._build_canonical_indices(eqs_sc)
    corr = CorrelationFunction(a', a, eqs_sc)
    ancilla_canon = QuantumCumulants._build_canonical_indices(corr.eqs)

    # The ancilla atom-canon list must be a subset of the parent's.
    atom_space = 2
    parent_pool = Set(get(parent_canon, atom_space, []))
    ancilla_pool = Set(get(ancilla_canon, atom_space, []))
    @test !isempty(parent_pool)
    @test issubset(ancilla_pool, parent_pool)

    # Stronger: every atom-space index appearing in an ancilla state's
    # operator must come from the parent's pool, not from a freshly
    # minted name like `j_3`.
    SQA = QuantumCumulants.SecondQuantizedAlgebra
    for s in corr.eqs.states
        op = SQA.undo_average(s)
        for idx in QuantumCumulants._all_indices(op)
            idx.space_index == atom_space || continue
            @test idx in parent_pool
        end
    end
end
