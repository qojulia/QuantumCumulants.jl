# PENDING: port of master test/test_indexed_filter_cavity.jl
#
# Status: blocked on three v1 gaps —
#   1. `evaluate(eqs; limits = (M => k), h = [hilb])` to unroll indexed sums
#      to a fixed size is not yet implemented (CHANGELOG marks as follow-up).
#   2. `scale(eqs; h = [k])` per-Hilbert-space scaling is not implemented.
#   3. The brute-force `ClusterSpace` comparison branch uses ClusterSpace,
#      which was removed in v1 (use Σ/Index directly).
#
# To enable: implement (1) and (2); rewrite the cluster-comparison branch
# to use indexed sums.

#=
using QuantumCumulants
using Symbolics: Symbolics, @variables
using ModelingToolkitBase: @named, mtkcompile, ODEProblem, unknowns
using OrdinaryDiffEq: Tsit5, solve
using Test

@testset "filter-cavity-indexing: closure + per-Hilbert-space scale" begin
    order = 2
    @variables κ::Real g::Real gf::Real κf::Real R::Real Γ::Real
    @variables Δ::Real ν::Real N::Real M::Real
    δ(i) = IndexedVariable(:δ, i)

    hc = FockSpace(:cavity)
    hf = FockSpace(:filter)
    ha = NLevelSpace(:atom, 2)
    h = hc ⊗ hf ⊗ ha

    i = Index(h, :i, M, hf)
    j = Index(h, :j, N, ha)

    @qnumbers a::Destroy(h, 1)
    b(k) = IndexedOperator(Destroy(h, :b, 2), k)
    σ(α, β, k) = IndexedOperator(Transition(h, :σ, α, β, 3), k)

    H = Δ * Σ(σ(2, 2, j), j) +
        Σ(δ(i) * b(i)' * b(i), i) +
        gf * Σ(a' * b(i) + a * b(i)', i) +
        g * Σ(a' * σ(1, 2, j) + a * σ(2, 1, j), j)

    J = [a, b(i), σ(1, 2, j), σ(2, 1, j), σ(2, 2, j)]
    rates = [κ, κf, Γ, R, ν]

    eqs = meanfield([a' * a], H, J; rates = rates, order = order)
    eqs_c = complete(eqs)
    @test length(eqs_c.equations) == 23

    M_ = 3
    eqs_eval_ = evaluate(eqs_c; limits = (M => M_), h = [2])
    @test length(eqs_eval_.equations) == 43
    eqs_sc_ = scale(eqs_eval_; h = [3])
    @test length(eqs_sc_.equations) == 42

    eqs_sc = scale(eqs_c; h = [3])
    @test length(eqs_sc.equations) == 22
    eqs_eval = evaluate(eqs_sc; limits = (M => M_))
    @test length(eqs_eval.equations) == 42
    @test length(eqs_sc_.equations) == length(eqs_eval.equations)

    # ODE comparison: scale-first vs evaluate-first must agree at steady state
    Γ_ = 1.0; Δ_ = 0Γ_; g_ = 2Γ_; κ_ = 1e3Γ_; R_ = 1e2Γ_; ν_ = 10Γ_
    gf_ = 0.1Γ_; κf_ = 0.1Γ_
    δ_ls = collect(0:(1 / M_):(1 - 1 / M_)) * 10Γ_
    ps = [Γ, κ, g, κf, gf, R, [δ(i) for i in 1:M_]..., Δ, ν, N]
    p0 = [Γ_, κ_, g_, κf_, gf_, R_, δ_ls..., Δ_, ν_, 1e4]

    @named sys = to_system(eqs_eval)
    sys_c = mtkcompile(sys)
    u0 = Dict(unknowns(sys_c) .=> zeros(ComplexF64, length(unknowns(sys_c))))
    prob = ODEProblem(sys_c, merge(u0, Dict(ps .=> p0)), (0.0, 1.0 / κf_))
    sol = solve(prob, Tsit5(); maxiters = 1e7, abstol = 1e-12, reltol = 1e-12)
    n_ind = get_solution(sol, a' * a, eqs_eval)(sol.t[end])

    @named sys2 = to_system(eqs_sc_)
    sys2_c = mtkcompile(sys2)
    u0_2 = Dict(unknowns(sys2_c) .=> zeros(ComplexF64, length(unknowns(sys2_c))))
    prob2 = ODEProblem(sys2_c, merge(u0_2, Dict(ps .=> p0)), (0.0, 1.0 / κf_))
    sol2 = solve(prob2, Tsit5(); maxiters = 1e7, abstol = 1e-12, reltol = 1e-12)
    n_ind2 = get_solution(sol2, a' * a, eqs_sc_)(sol2.t[end])

    @test n_ind ≈ n_ind2
end
=#
