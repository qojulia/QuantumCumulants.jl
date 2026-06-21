using QuantumCumulants
using SecondQuantizedAlgebra: SecondQuantizedAlgebra
using Symbolics: Symbolics, @variables
using Test
const SQA = SecondQuantizedAlgebra

@testset "closure: find_missing == 0 after complete, scale, evaluate" begin
    ha = NLevelSpace(:atom, 2); hf = FockSpace(:cavity); h = ha ⊗ hf
    a = Destroy(h, :a)
    σ(α, β, i) = IndexedOperator(Transition(h, :σ, α, β), i)
    @variables N::Real Δ::Real g::Real κ::Real
    i = Index(h, :i, N, ha)
    H = Δ * a' * a + g * ∑(σ(2, 1, i) * a + σ(1, 2, i) * a', i)
    eqs = meanfield(a' * a, H, [a]; rates = [κ], order = 2)
    complete!(eqs)
    @test isempty(find_missing(eqs))
    sc = scale(eqs)
    @test isempty(find_missing(sc))
    ev = evaluate(eqs; limits = (N => 3))
    @test isempty(find_missing(ev))
end

@testset "closure: get_adjoints=false stores conjugation-canonical reps (issue #295)" begin
    # Which member of a conjugate pair survives must not depend on objectid-seeded
    # leaf order: every stored moment is its own conjugation-canonical representative.
    ha = NLevelSpace(:atom, 2); hf = FockSpace(:cavity); h = ha ⊗ hf
    a = Destroy(h, :a)
    σ(α, β, i) = IndexedOperator(Transition(h, :σ, α, β), i)
    @variables N::Real Δ::Real g::Real κ::Real
    i = Index(h, :i, N, ha)
    H = Δ * a' * a + g * ∑(σ(2, 1, i) * a + σ(1, 2, i) * a', i)
    eqs = meanfield(a' * a, H, [a]; rates = [κ], order = 2)
    eqs_c = complete(eqs; get_adjoints = false)
    ctx = QuantumCumulants.build_ctx(eqs_c)
    # Non-vacuous only if the system has genuine conjugate pairs (non-Hermitian moments).
    @test any(
        op -> QuantumCumulants.canon_key(adjoint(op), ctx) != QuantumCumulants.canon_key(op, ctx),
        keys(eqs_c.graph.nodes),
    )
    # is_conjugate == false ⟺ op is the canonical side of its pair.
    for op in keys(eqs_c.graph.nodes)
        @test !QuantumCumulants.canonical_rep(op, ctx)[2]
    end
end

@testset "closure: chain N=10 find_missing == 0 (coordinate-consistent)" begin
    N = 10
    h_chain = ⊗([NLevelSpace(Symbol(:atom, i), (:g, :e)) for i in 1:N]...)
    σ_ch(i, j, k) = Transition(h_chain, Symbol(:σ_, k), i, j, k)
    @variables Ω γ Δ J0
    x_ = [first(@variables $(Symbol("x_$i"))) for i in 1:N]
    Jc(xi, xj) = J0 / abs(xi - xj)^3
    H = -Δ * sum(σ_ch(:e, :e, k) for k in 1:N) +
        Ω * (σ_ch(:e, :g, 1) + σ_ch(:g, :e, 1)) +
        sum(Jc(x_[k], x_[k + 1]) * (σ_ch(:e, :g, k) * σ_ch(:g, :e, k + 1) + σ_ch(:g, :e, k) * σ_ch(:e, :g, k + 1)) for k in 1:(N - 1))
    c_ops = [σ_ch(:g, :e, k) for k in 1:N]
    eqs = meanfield(σ_ch(:g, :e, 1), H, c_ops; rates = [γ for _ in 1:N], order = 2)
    complete!(eqs)
    # find_missing must agree with the resolver: the system closes at 373, so
    # find_missing must be empty (was 54 before Task 2 unified the code paths).
    @test isempty(find_missing(eqs))
end
