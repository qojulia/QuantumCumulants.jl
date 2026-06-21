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

@testset "scale: conjugation-folded count is the minimal set (issue #295)" begin
    # `scale` of a `get_adjoints=false` system must fold conjugate images: no two scaled
    # states may share a conjugation representative, regardless of objectid-seeded leaf order.
    ha = NLevelSpace(:atom, 2); hf = FockSpace(:cavity); h = ha ⊗ hf
    a = Destroy(h, :a)
    σ(α, β, i) = IndexedOperator(Transition(h, :σ, α, β), i)
    @variables N::Real Δ::Real g::Real κ::Real
    i = Index(h, :i, N, ha)
    H = Δ * a' * a + g * ∑(σ(2, 1, i) * a + σ(1, 2, i) * a', i)
    eqs = meanfield(a' * a, H, [a]; rates = [κ], order = 2)
    sc = scale(complete(eqs; get_adjoints = false))
    ctx = QuantumCumulants.build_ctx(sc)
    treatments = QuantumCumulants._treatments(sc, ctx)
    reps = [
        QuantumCumulants.canonical_rep(QuantumCumulants.undo_average(s), ctx; treatments)[1]
        for s in sc.states
    ]
    @test any(s -> QuantumCumulants.undo_average(s) isa SQA.QAdd, sc.states)  # non-vacuous
    @test allunique(reps)
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
