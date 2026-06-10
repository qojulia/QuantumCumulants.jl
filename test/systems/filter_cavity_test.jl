using QuantumCumulants
using ModelingToolkitBase
using Symbolics: @variables
using OrdinaryDiffEqTsit5: ODEProblem, solve, Tsit5
using Test

@testset "filter-cavity_indexed" begin
    @variables κ g gf κf R Γ Δ ν N M
    δ_v(idx) = IndexedVariable(:δ, idx)
    hc = FockSpace(:cavity); hf = FockSpace(:filter); ha = NLevelSpace(:atom, 2)
    h = hc ⊗ hf ⊗ ha
    i = Index(h, :i, M, hf); j = Index(h, :j, N, ha)
    @qnumbers a::Destroy(h, 1)
    b(idx) = IndexedOperator(Destroy(h, :b, 2), idx)
    b(idx::Integer) = IndexedOperator(Destroy(h, :b, 2), i(2)(idx))
    σ(α, β, idx) = IndexedOperator(Transition(h, :σ, α, β, 3), idx)
    σ(α, β, idx::Integer) = IndexedOperator(Transition(h, :σ, α, β, 3), j(idx))
    H = Δ * Σ(σ(2, 2, j), j) +
        Σ(δ_v(i) * b(i)' * b(i), i) +
        gf * (Σ(a' * b(i) + a * b(i)', i)) +
        g * (Σ(a' * σ(1, 2, j) + a * σ(2, 1, j), j))
    J = [a, b(i), σ(1, 2, j), σ(2, 1, j), σ(2, 2, j)]
    rates = [κ, κf, Γ, R, ν]
    eqs = meanfield(a' * a, H, J; rates = rates, order = 2)
    eqs_c = complete!(deepcopy(eqs))
    # closure: 42 equations, 39 after scale, 552 after evaluate at M=20.
    # (552 is the conjugate-folded distinct-site count; verified equal to master
    # QuantumCumulants v0.4.3 for this exact system. The earlier 554 lock came
    # from a buggy QuantumCumulants scale that double-counted a filter conjugate pair.)
    @test length(eqs_c.equations) == 42
    eqs_sc = scale(eqs_c; h = [3])
    @test length(eqs_sc.equations) == 39
    eqs_eval = evaluate(eqs_sc; limits = Dict(M => 20))
    @test length(eqs_eval.equations) == 552
    sys_c = mtkcompile(System(eqs_eval; name = :fc_sys))
    @test length(unknowns(sys_c)) == 552
end

# Reduced M=3 variant integrated as a mean-field ODE.
@testset "filter-cavity_indexed (M=3 numerical)" begin
    @variables κ g gf κf R Γ Δ ν N M
    δ_v(idx) = IndexedVariable(:δ, idx)
    hc = FockSpace(:cavity); hf = FockSpace(:filter); ha = NLevelSpace(:atom, 2)
    h = hc ⊗ hf ⊗ ha
    i = Index(h, :i, M, hf); j = Index(h, :j, N, ha)
    @qnumbers a::Destroy(h, 1)
    b(idx) = IndexedOperator(Destroy(h, :b, 2), idx)
    b(idx::Integer) = IndexedOperator(Destroy(h, :b, 2), i(2)(idx))
    σ(α, β, idx) = IndexedOperator(Transition(h, :σ, α, β, 3), idx)
    σ(α, β, idx::Integer) = IndexedOperator(Transition(h, :σ, α, β, 3), j(idx))
    H = Δ * Σ(σ(2, 2, j), j) +
        Σ(δ_v(i) * b(i)' * b(i), i) +
        gf * (Σ(a' * b(i) + a * b(i)', i)) +
        g * (Σ(a' * σ(1, 2, j) + a * σ(2, 1, j), j))
    J = [a, b(i), σ(1, 2, j), σ(2, 1, j), σ(2, 2, j)]
    rates = [κ, κf, Γ, R, ν]
    eqs = meanfield(a' * a, H, J; rates = rates, order = 2)
    eqs_c = complete!(deepcopy(eqs))
    eqs_sc = scale(eqs_c; h = [3])
    M_ = 3
    eqs_eval = evaluate(eqs_sc; limits = Dict(M => M_))
    @test length(eqs_eval.equations) == 42
    sys_c = mtkcompile(System(eqs_eval; name = :fc3_sys))
    u0 = zeros(ComplexF64, length(eqs_eval.equations))
    init = initial_values(eqs_eval, u0)
    pmap_dict = Dict{Any, Any}(
        κ => 1.0, κf => 0.05, g => 0.5, gf => 0.05,
        R => 1.0, Γ => 0.01, Δ => 0.0, ν => 0.01, N => 10.0,
    )
    δ_vals = [-1.0, 0.0, 1.0]
    for k in 1:M_
        pmap_dict[δ_v(i(2)(k))] = δ_vals[k]
    end
    pmap = parameter_map(eqs_eval, pmap_dict)
    prob = ODEProblem(sys_c, merge(init, pmap), (0.0, 5.0))
    sol = solve(prob, Tsit5(); reltol = 1.0e-8, abstol = 1.0e-8)
    # Ground-truth ⟨a†a⟩(t=5): verified against master QuantumCumulants v0.4.3
    # (both its scaled and its fully-unrolled N=10 M=3 reference) AND QuantumCumulants's own
    # fully-unrolled system, all agreeing at 4.36032. The earlier 4.40621 lock was
    # produced by a buggy QuantumCumulants scale (filter NE-drop + δ slot mis-binding).
    @test isapprox(
        real(get_solution(sol, a' * a, eqs_eval).(sol.t[end])),
        4.360320310008423; rtol = 1.0e-4,
    )
end
