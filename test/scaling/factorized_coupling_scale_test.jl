using QuantumCumulants
using Symbolics: Symbolics, @variables
using SymbolicUtils: SymbolicUtils
using ModelingToolkitBase: ModelingToolkitBase, mtkcompile, ODEProblem, unknowns, parameters
using OrdinaryDiffEqTsit5: Tsit5, solve, ReturnCode
using Test

# Companion to issue #288 on the `scale` path. A two-atom coupling
# `Σ_{i≠j} J(i,j)·(a† σ⁻ᵢ σ⁻ⱼ + h.c.)` on an index-free observable `a` produces, at order 1,
# a factorised sum body `Σ_{i≠j} J(i,j)·⟨σ⁻ᵢ⟩⟨σ⁻ⱼ⟩` with BOTH summation indices bound in one
# indexed-sum node. `scale` must fold this per-leaf (each atom collapses to the representative)
# and charge the falling-factorial `N(N−1)` prefactor, NOT `undo_average` the whole body — which
# would fuse the moment product `⟨σ⟩⟨σ⟩` into the operator `σσ` and re-average to a spurious
# order-2 moment `⟨σσ⟩`, leaving the system unclosed.
@testset "scale: factorized two-index coupling stays at truncation order" begin
    ha = NLevelSpace(:atom, 2); hb = FockSpace(:cav); h = hb ⊗ ha
    a = Destroy(h, :a)
    @variables N::Real κ::Real
    i = Index(h, :i, N, ha); j = Index(h, :j, N, ha)
    J(x, y) = DoubleIndexedVariable(:J, x, y)
    sp(idx) = IndexedOperator(Transition(h, :σ, 2, 1), idx)
    sm(idx) = IndexedOperator(Transition(h, :σ, 1, 2), idx)
    H = Σ(Σ(J(i, j) * (a' * sm(i) * sm(j) + a * sp(i) * sp(j)), j, [i]), i)
    eqs_c = complete(meanfield(a, H, [a]; rates = [κ], order = 1))

    eqs_sc = scale(eqs_c; h = [2])

    # The scaled system closes and no equation exceeds the truncation order.
    @test isempty(find_missing(eqs_sc))
    @test all(eq -> get_order(eq.rhs) <= 1, eqs_sc.equations)
end

@testset "scale: factorized coupling agrees with evaluate" begin
    # The scaled (symmetry-reduced) and evaluated (unrolled at N=3) systems must give the
    # same dynamics; this pins the falling-factorial prefactor and the per-leaf folding.
    ha = NLevelSpace(:atom, 2); hb = FockSpace(:cav); h = hb ⊗ ha
    a = Destroy(h, :a)
    @variables N::Real κ::Real
    i = Index(h, :i, N, ha); j = Index(h, :j, N, ha)
    J(x, y) = DoubleIndexedVariable(:J, x, y)
    sp(idx) = IndexedOperator(Transition(h, :σ, 2, 1), idx)
    sm(idx) = IndexedOperator(Transition(h, :σ, 1, 2), idx)
    H = Σ(Σ(J(i, j) * (a' * sm(i) * sm(j) + a * sp(i) * sp(j)), j, [i]), i)
    eqs_c = complete(meanfield(a, H, [a]; rates = [κ], order = 1))

    Nv, Jv, κv = 3, 0.35, 0.7
    eqs_sc = scale(eqs_c; h = [2])
    eqs_ev = evaluate(eqs_c; limits = (N => Nv))
    @test isempty(find_missing(eqs_sc))
    @test isempty(find_missing(eqs_ev))

    sys_sc = mtkcompile(System(eqs_sc; name = :sc))
    sys_ev = mtkcompile(System(eqs_ev; name = :ev))
    init_sc = initial_values(eqs_sc, fill(ComplexF64(0.05), length(eqs_sc.equations)))
    init_ev = initial_values(eqs_ev, fill(ComplexF64(0.05), length(eqs_ev.equations)))
    psc = Dict{Any, Any}()
    for p in parameters(sys_sc)
        nm = string(p)
        nm == "J" && (psc[p] = Jv)
        nm == "κ" && (psc[p] = κv)
        nm == "N" && (psc[p] = Nv)
    end
    pev = Dict{Any, Any}()
    for p in parameters(sys_ev)
        nm = string(p)
        nm == "κ" && (pev[p] = κv)
        nm == "J" && (pev[p] = [x == y ? 0.0 : Jv for x in 1:Nv, y in 1:Nv])
    end
    sol_sc = solve(ODEProblem(sys_sc, merge(init_sc, psc), (0.0, 3.0)), Tsit5(); reltol = 1.0e-10, abstol = 1.0e-12)
    sol_ev = solve(ODEProblem(sys_ev, merge(init_ev, pev), (0.0, 3.0)), Tsit5(); reltol = 1.0e-10, abstol = 1.0e-12)
    @test sol_sc.retcode == ReturnCode.Success
    @test sol_ev.retcode == ReturnCode.Success
    for τ in (0.5, 1.0, 2.0, 3.0)
        v_sc = get_solution(sol_sc, a, eqs_sc)(τ)
        v_ev = get_solution(sol_ev, a, eqs_ev)(τ)
        @test isapprox(v_sc, v_ev; atol = 1.0e-8)
    end
end

@testset "scale: multi-average sum body spanning a scaled + a non-scaled subsystem" begin
    # A bilinear filter–atom coupling `Σ_{i,j} J(i,j)·(b_i† σ⁻_j + h.c.)` factorises at order 1
    # into `Σ_{i,j} J(i,j)·⟨b_i⟩⟨σ⁻_j⟩`: one sum node whose scope mixes a filter index `i` and
    # an atom index `j`. Scaling ONLY the atoms collapses `j` (falling-factorial prefactor) but
    # must keep `i` in a residual `Σ_i`. This exercises the `_scale_sum` residual-sum branch;
    # the result must agree with unrolling both subsystems directly.
    hf = FockSpace(:filter); ha = NLevelSpace(:atom, 2); h = hf ⊗ ha
    @variables M::Real N::Real κ::Real κf::Real
    i = Index(h, :i, M, hf); j = Index(h, :j, N, ha)
    bmode(x) = IndexedOperator(Destroy(h, :b), x)
    J(x, y) = DoubleIndexedVariable(:J, x, y)
    σm(x) = IndexedOperator(Transition(h, :σ, 1, 2), x)
    σz(x) = IndexedOperator(Transition(h, :σ, 2, 2), x)
    σp(x) = IndexedOperator(Transition(h, :σ, 2, 1), x)
    H = Σ(Σ(J(i, j) * (bmode(i)' * σm(j) + bmode(i) * σp(j)), j), i)
    eqs_c = complete(
        meanfield([bmode(i)' * bmode(i), σz(j)], H, [bmode(i), σm(j)]; rates = [κf, κ], order = 1),
    )
    atom = QuantumCumulants.SecondQuantizedAlgebra.acts_on(σz(j))  # atom subspace index

    Mv, Nv, Jv, κv, κfv = 2, 3, 0.25, 0.5, 0.3
    eqs_A = evaluate(scale(eqs_c; h = atom); limits = (M => Mv))   # scale atoms, unroll filter
    eqs_B = evaluate(eqs_c; limits = (M => Mv, N => Nv))            # unroll both
    @test isempty(find_missing(eqs_A))
    @test isempty(find_missing(eqs_B))
    @test all(eq -> get_order(eq.rhs) <= 1, eqs_A.equations)

    SQA = QuantumCumulants.SecondQuantizedAlgebra
    isarr(p) = SymbolicUtils.symtype(SymbolicUtils.unwrap(p)) <: AbstractArray
    function run(es)
        sys = mtkcompile(System(es; name = :m))
        init = initial_values(es, fill(ComplexF64(0.1), length(es.equations)))
        pm = Dict{Any, Any}()
        for p in parameters(sys)
            nm = string(p)
            nm == "κ" && (pm[p] = κv)
            nm == "κf" && (pm[p] = κfv)
            nm == "N" && (pm[p] = Nv)
            nm == "M" && (pm[p] = Mv)
            occursin("J", nm) && (pm[p] = isarr(p) ? fill(Jv, Mv, Nv) : Jv)
        end
        solve(ODEProblem(sys, merge(init, pm), (0.0, 3.0)), Tsit5(); reltol = 1.0e-10, abstol = 1.0e-12)
    end
    sol_A = run(eqs_A); sol_B = run(eqs_B)
    @test sol_A.retcode == ReturnCode.Success
    @test sol_B.retcode == ReturnCode.Success

    # Compare a filter-mode population ⟨b_i† b_i⟩ common to both representations.
    ops_A = [SQA.undo_average(s) for s in eqs_A.states]
    target = ops_A[
        findfirst(
            o -> occursin("b_", string(o)) && occursin("'", string(o)) && !occursin("σ", string(o)),
            ops_A,
        ),
    ]
    ops_B = [SQA.undo_average(s) for s in eqs_B.states]
    kB = findfirst(o -> isequal(o, target), ops_B)
    @test kB !== nothing
    for τ in (0.5, 1.0, 2.0, 3.0)
        @test isapprox(
            get_solution(sol_A, target, eqs_A)(τ),
            get_solution(sol_B, ops_B[kB], eqs_B)(τ);
            atol = 1.0e-8,
        )
    end
end
