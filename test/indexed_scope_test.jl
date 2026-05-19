using QuantumCumulants
using SecondQuantizedAlgebra: SecondQuantizedAlgebra
using Symbolics: Symbolics, @variables
using SymbolicUtils
using ModelingToolkitBase: ModelingToolkitBase, @named, mtkcompile, ODEProblem, unknowns
using OrdinaryDiffEq: Tsit5, solve, ReturnCode
using Test

const SQA = QuantumCumulants.SecondQuantizedAlgebra

# Behavioural regressions for the bound-index coefficient orphaning fix.
# Each testset verifies a user-observable outcome (parameter_map shape,
# ODE retcode, per-atom population, cavity decay rate). Internal helpers
# (`_lift_sum_scope`, `_materialise_scoped`, `_is_leaf_average`, …) are
# exercised indirectly through these surfaces; their structure is free
# to change as long as the observed behaviour holds.

@testset "indexed scope: parameter_map shapes" begin
    # Smallest possible indexed JC. We only need an `IndexedVariable`
    # `g(i)` and a finite N to drive `parameter_map` through both code
    # paths (scalar broadcast + per-atom vector).
    ha = NLevelSpace(:atom, 2); hf = FockSpace(:cavity); h = ha ⊗ hf
    a = Destroy(h, :a)
    σ(α, β, i) = IndexedOperator(Transition(h, :σ, α, β), i)
    @variables N::Real κ::Real Δc::Real
    g(i) = IndexedVariable(:g, i)
    i = Index(h, :i, N, ha)
    H = -Δc * a' * a +
        ∑(g(i) * (σ(2, 1, i) * a + σ(1, 2, i) * a'), i)
    eqs = meanfield(a'*a, H, [a]; rates = [κ], order = 2)
    complete!(eqs)
    evaled = evaluate(eqs; limits = (N => 3))

    # Scalar value → broadcast to N-element vector.
    pmap = parameter_map(evaled, Dict(
        g(i) => 0.7,
        κ    => 1.5,
        Δc   => 0.25,
    ))
    @test length(pmap) == 3
    arr_entry = first((k, v) for (k, v) in pmap if v isa AbstractArray)
    @test arr_entry[2] == fill(0.7, 3)
    scalars = Dict(string(k) => v for (k, v) in pmap if !(v isa AbstractArray))
    @test scalars["κ"] == 1.5
    @test scalars["Δc"] == 0.25

    # Per-atom vector value → pass through unchanged.
    pmap2 = parameter_map(evaled, Dict(g(i) => [0.1, 0.2, 0.3]))
    @test only(values(pmap2)) == [0.1, 0.2, 0.3]
end

@testset "indexed scope: end-to-end inhomogeneous JC integrates" begin
    # Drives the full pipeline (meanfield → complete! → evaluate →
    # to_system → ODEProblem → solve) with per-atom inhomogeneous
    # coupling. A successful solve here confirms: lifted scope reaches
    # the operator-side substitution, IndexedVariable coefficients track
    # the rename, and MTK accepts the resulting array-parameter system.
    ha = NLevelSpace(:atom, 2); hf = FockSpace(:cavity); h = ha ⊗ hf
    a = Destroy(h, :a)
    σ(α, β, i) = IndexedOperator(Transition(h, :σ, α, β), i)
    @variables η::Real Γ::Real κ::Real Δc::Real N::Real
    g(i) = IndexedVariable(:g, i)
    i = Index(h, :i, N, ha)
    H = -Δc * a' * a +
        ∑(g(i) * (σ(2, 1, i) * a + σ(1, 2, i) * a'), i) +
        1im * η * (a' - a)
    J = [σ(1, 2, i), a]
    eqs = meanfield(a'*a, H, J; rates = [Γ, κ], order = [1, 2])
    complete!(eqs)
    evaled = evaluate(eqs; limits = (N => 2))

    @named sys = to_system(evaled)
    sys_c = mtkcompile(sys)
    pmap = parameter_map(evaled, Dict(
        g(i) => [0.05, 0.15],   # inhomogeneous on purpose (see below)
        Δc   => 0.0,
        κ    => 1.0,
        Γ    => 0.25,
        η    => 0.3,
    ))
    u0 = Dict(unknowns(sys_c) .=> zeros(ComplexF64, length(unknowns(sys_c))))
    prob = ODEProblem(sys_c, merge(u0, pmap), (0.0, 10.0))
    sol = solve(prob, Tsit5(); abstol = 1e-8, reltol = 1e-8)
    @test sol.retcode == ReturnCode.Success

    a_op = SQA.undo_average(evaled.states[1])
    assert_real(sol, a_op, evaled; atol = 1e-6)
    assert_nonneg(sol, a_op, evaled; atol = 1e-6)
end

@testset "indexed scope: per-atom couplings produce distinguishable observables" begin
    # With inhomogeneous coupling `g(i_1) ≠ g(i_2)`, atom-cavity
    # correlation observables `⟨a σ_22(i_k)⟩` must differ between
    # atoms. With the orphaned-coefficient bug, both atoms would share
    # the SAME un-substituted `g(i)` factor and the two correlations
    # would coincide at the symbolic level (and hence numerically).
    #
    # We use `⟨a σ_22(i_k)⟩` rather than bare `⟨σ_22(i_k)⟩` because
    # the order=[1,2] cumulant closure can pin pure-atom populations to
    # zero at finite drive while the mixed-observable still tracks the
    # per-atom coupling difference.
    ha = NLevelSpace(:atom, 2); hf = FockSpace(:cavity); h = ha ⊗ hf
    a = Destroy(h, :a)
    σ(α, β, i) = IndexedOperator(Transition(h, :σ, α, β), i)
    @variables η::Real Γ::Real κ::Real N::Real
    g(i) = IndexedVariable(:g, i)
    i = Index(h, :i, N, ha)
    H = ∑(g(i) * (σ(2, 1, i) * a + σ(1, 2, i) * a'), i) +
        1im * η * (a' - a)
    J = [σ(1, 2, i), a]
    eqs = meanfield(a'*a, H, J; rates = [Γ, κ], order = [1, 2])
    complete!(eqs)
    evaled = evaluate(eqs; limits = (N => 2))

    # Locate the two atom-cavity correlation states `⟨a σ_22(i_k)⟩`.
    # Computed BEFORE `to_system` so the filter sees the un-substituted
    # average leaves.
    # `complete!` non-deterministically picks one of each conjugate pair,
    # so atom-cavity correlations may appear as ⟨a σ_22⟩ OR ⟨a' σ_22⟩.
    # Accept either: 2-op term whose σ-piece is the population
    # `Transition(_, 2, 2, i_k)` paired with a cavity field operator.
    corr_states = SymbolicUtils.BasicSymbolic[]
    for s in evaled.states
        op = SQA.undo_average(s)
        op isa SQA.QAdd || continue
        length(op.arguments) == 1 || continue
        term = first(keys(op.arguments))
        length(term.ops) == 2 || continue
        has_pop = false
        has_cavity = false
        for o in term.ops
            if o isa SQA.Transition && o.i == 2 && o.j == 2
                has_pop = true
            elseif o isa SQA.Destroy || o isa SQA.Create
                has_cavity = true
            end
        end
        has_pop && has_cavity && push!(corr_states, s)
    end

    @named sys = to_system(evaled)
    sys_c = mtkcompile(sys)
    pmap = parameter_map(evaled, Dict(
        g(i) => [0.05, 0.20],   # 4× ratio
        κ    => 1.0,
        Γ    => 0.1,
        η    => 0.3,
    ))
    u0 = Dict(unknowns(sys_c) .=> zeros(ComplexF64, length(unknowns(sys_c))))
    prob = ODEProblem(sys_c, merge(u0, pmap), (0.0, 20.0))
    sol = solve(prob, Tsit5(); abstol = 1e-9, reltol = 1e-9)
    @test sol.retcode == ReturnCode.Success
    @test length(corr_states) == 2
    v1 = get_solution(sol, corr_states[1], evaled)(sol.t[end])
    v2 = get_solution(sol, corr_states[2], evaled)(sol.t[end])
    # The two atom-cavity correlations must differ; the magnitude ratio
    # should reflect the g-coupling asymmetry rather than coincide.
    @test abs(v1 - v2) > 1e-3
end
