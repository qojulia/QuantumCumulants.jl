# Allocation regression tests for the symbolic hot spots.
#
# These passes allocate by design (they build expression trees, dicts and arrays);
# the goal is a deterministic gate that fails when a hot spot's allocation count
# grows, not zero allocations. See
# docs/superpowers/specs/2026-06-15-allocation-regression-tests-design.md.

using QuantumCumulants
using Symbolics: @variables
using Test

import QuantumCumulants.SecondQuantizedAlgebra as SQA
const QC = QuantumCumulants

# Ceilings are calibrated on Julia 1.12 (matching benchmark/Benchmarks.yaml) at
# floor * 1.15, where the floor is the deterministic minimum allocation after
# warmup. On other versions the ceilings are skipped rather than asserted: if a
# newer Julia drifts allocations, either rebaseline the constants below or tighten
# the gate to a band such as `v"1.12" <= VERSION < v"1.13"`.
const ALLOC_GATE = VERSION >= v"1.12"

const CEILING = Dict(
    :cumulant_expansion => 79_065,
    :canonicalization => 167_182,
    :evaluate => 768_715,
    :scale => 176_934,
    :System => 323_168,
    :spectrum_build => 243_128,
    :spectrum_kernel => 66_369,
)

# Warm up to absorb compilation (and first-call interning, e.g. System), then take
# the minimum over repeats: the symbolic passes are deterministic, so this is a
# stable floor with GC jitter removed.
function allocs(f)
    f()
    f()
    return minimum(@allocated(f()) for _ in 1:6)
end

function check_alloc(name, f)
    return if ALLOC_GATE
        @test allocs(f) <= CEILING[name]
    else
        @info "allocation ceiling skipped (needs Julia >= 1.12)" hot_spot = name VERSION
        @test_skip allocs(f) <= CEILING[name]
    end
end

# Fixtures =======================================================================

# Jaynes-Cummings model: cumulant expansion, canonicalization, System.
function jc_fixture()
    hf = FockSpace(:cavity); ha = NLevelSpace(:atom, (:g, :e)); h = hf ⊗ ha
    a = Destroy(h, :a); σ = Transition(h, :σ, :g, :e)
    @variables Δ g κ γ ν
    H = Δ * a' * a + g * (a' * σ + σ' * a)
    J = [a, σ, σ']; rates = [κ, γ, ν]; ops = [a' * a, σ' * σ, a * σ']
    raw = meanfield(ops, H, J; rates = rates)
    completed = complete(cumulant_expansion(raw, 2))
    return (; raw, ops, H, J, completed)
end

# Phase-invariance filter for the indexed superradiant laser.
_φ(x) = 0
_φ(::Destroy) = -1
_φ(::Create) = 1
_φ(x::Transition) = x.i - x.j
function _φ(q::SQA.QAdd)
    for (term, _) in q.arguments
        return sum(_φ(op) for op in term.ops; init = 0)
    end
    return 0
end
_φ(avg) = SQA.is_average(avg) ? _φ(SQA.undo_average(avg)) : 0
phase_invariant(x) = iszero(_φ(x))

# Indexed superradiant laser: scale and evaluate.
function superradiant_fixture()
    hc = FockSpace(:cavity); ha = NLevelSpace(:atom, 2); h = hc ⊗ ha
    a = Destroy(h, :a); σ(α, β, i) = IndexedOperator(Transition(h, :σ, α, β), i)
    @variables N Δ κ Γ R ν
    g(i) = IndexedVariable(:g, i); i = Index(h, :i, N, ha); j = Index(h, :j, N, ha)
    H = -Δ * a'a + Σ(g(i) * (a' * σ(1, 2, i) + a * σ(2, 1, i)), i)
    J = [a, σ(1, 2, i), σ(2, 1, i), σ(2, 2, i)]; rates = [κ, Γ, R, ν]
    ops = [a' * a, σ(2, 2, j)]
    completed = complete(meanfield(ops, H, J; rates = rates, order = 2); filter_func = phase_invariant)
    return (; completed, N)
end

# Phase-invariance filter for the single-atom laser spectrum (matches spectrum_kernel_test.jl).
_pi(x) = 0
_pi(::Destroy) = -1
_pi(::Create) = 1
_pi(t::Transition) = (t.i == 2 && t.j == 1) ? 1 : (t.i == 1 && t.j == 2) ? -1 : 0
function _pi(q::SQA.QAdd)
    for (term, _) in q.arguments
        return sum(_pi(op) for op in term.ops; init = 0)
    end
    return 0
end
_pi(avg) = SQA.is_average(avg) ? _pi(SQA.undo_average(avg)) : 0
pinv(x) = iszero(_pi(x))

# Single-atom laser spectrum kernel: the kernel is a pure function of the symbolic
# τ-system, so block-level checks only need a fabricated (u_end, p0).
function spectrum_fixture()
    @variables Δ g γ κ ν
    hf = FockSpace(:cavity); ha = NLevelSpace(:atom, (:g, :e)); h = hf ⊗ ha
    a = Destroy(h, :a); s = Transition(h, :σ, :g, :e)
    H = Δ * a' * a + g * (a' * s + a * s')
    eqs = complete(meanfield(a' * a, H, [a, s, s']; rates = [κ, γ, ν], order = 2); filter_func = pinv)
    c = CorrelationFunction(a', a, eqs; steady_state = true, filter_func = pinv)
    S = Spectrum(c, (Δ, g, γ, κ, ν))
    p0 = (0.0, 1.5, 0.25, 1.0, 4.0)
    u_end = Dict{Any, ComplexF64}()
    for (i, st) in enumerate(eqs.states)
        u_end[st] = (0.3 + 0.1i) + (0.2 - 0.05i) * im
    end
    return S, u_end, p0
end

# Tests ==========================================================================

@testset "allocation regressions" begin
    jc = jc_fixture()
    sr = superradiant_fixture()
    S, u_end, p0 = spectrum_fixture()

    @testset "cumulant_expansion" begin
        check_alloc(:cumulant_expansion, () -> cumulant_expansion(jc.raw, 2))
    end

    @testset "canonicalization" begin
        # Cold path: a fresh ctx (fresh cache) canonicalising every completed moment.
        # canonical_rep is memoised in ctx.cache, so a shared ctx would measure only
        # the cache lookup; rebuilding the ctx each call measures real work.
        canon_cold = function ()
            ctx = QC.build_ctx(jc.ops, jc.H, jc.J, adjoint.(jc.J))
            s = 0
            for st in jc.completed.states
                s += length(QC.canonical_rep(SQA.undo_average(st), ctx))
            end
            return s
        end
        check_alloc(:canonicalization, canon_cold)
    end

    @testset "evaluate" begin
        check_alloc(:evaluate, () -> evaluate(sr.completed; limits = (sr.N => 2)))
    end

    @testset "scale" begin
        check_alloc(:scale, () -> scale(sr.completed))
    end

    @testset "System" begin
        check_alloc(:System, () -> System(jc.completed; name = :alloc))
    end

    @testset "spectrum kernel" begin
        # Build path: each call recompiles the kernel from the symbolic τ-system.
        check_alloc(:spectrum_build, () -> QC._build_spectrum_kernel(S))
        # Reuse path: warmup populates S.cache, so the measured calls hit the cache.
        check_alloc(:spectrum_kernel, () -> QC._spectrum_kernel(S, u_end, p0))
    end
end
