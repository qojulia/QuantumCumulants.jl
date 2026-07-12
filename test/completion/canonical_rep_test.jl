using QuantumCumulants
using SecondQuantizedAlgebra: SecondQuantizedAlgebra
using Symbolics: Symbolics, @variables
using Test

const SQA = QuantumCumulants.SecondQuantizedAlgebra

@testset "canonical_rep: config wrappers equal the legacy keys" begin
    ha = NLevelSpace(:atom, 2); hf = FockSpace(:cavity); h = ha ⊗ hf
    a = Destroy(h, :a)
    σ(α, β, i) = IndexedOperator(Transition(h, :σ, α, β), i)
    @variables N::Real
    i = Index(h, :i, N, ha); j = Index(h, :j, N, ha)
    ctx = QuantumCumulants.build_ctx([a], a' * a, SQA.QField[], SQA.QField[])

    op = (σ(2, 1, i) * a) * 1
    # canonical_rep ADDS the conjugation fold on top of canon_key, so under All-Free
    # treatment it equals the conjugation-min of canon_key over {op, op†}. Which side
    # wins is SQA's qadd_order_key convention, so derive it rather than hard-coding op.
    free = QuantumCumulants.all_free_treatments(ctx)
    ck = QuantumCumulants.canon_key(op, ctx)
    cka = QuantumCumulants.canon_key(adjoint(op), ctx)
    expected = SQA.qadd_order_key(cka) < SQA.qadd_order_key(ck) ? cka : ck
    @test isequal(QuantumCumulants.canonical_rep(op, ctx; treatments = free)[1], expected)
    # idempotence
    rep1 = QuantumCumulants.canonical_rep(op, ctx; treatments = free)[1]
    @test isequal(QuantumCumulants.canonical_rep(rep1 * 1, ctx; treatments = free)[1], rep1)
    # adjoint-consistency: rep is invariant under {O, O†}, sign bit flips
    r, s = QuantumCumulants.canonical_rep(op, ctx; treatments = free)
    ra, sa = QuantumCumulants.canonical_rep(adjoint(op), ctx; treatments = free)
    @test isequal(r, ra)
    @test s != sa
end

@testset "canonical_rep: Scaled treatment equals scaled_key" begin
    ha = NLevelSpace(:atom, 2); hf = FockSpace(:cavity); h = ha ⊗ hf
    a = Destroy(h, :a)
    σ(α, β, i) = IndexedOperator(Transition(h, :σ, α, β), i)
    @variables N::Real g::Real
    i = Index(h, :i, N, ha); j = Index(h, :j, N, ha)
    H = g * ∑(σ(2, 1, i) * a + σ(1, 2, i) * a', i)
    ctx = QuantumCumulants.build_ctx([a], H, SQA.QField[], SQA.QField[])
    op = (σ(2, 2, i) * σ(2, 1, j)) * 1
    atom_sp = first(ctx.symmetric)
    treatments = QuantumCumulants.with_treatment(QuantumCumulants.all_free_treatments(ctx), atom_sp, QuantumCumulants.Scaled)
    # canonical_rep ADDS the conjugation fold on top of the scaled_key node key,
    # so it equals the conjugation-min of scaled_key over {op, op†}.
    ok = QuantumCumulants.scaled_key(op, ctx)
    oka = QuantumCumulants.scaled_key(adjoint(op), ctx)
    expected = SQA.qadd_order_key(oka) < SQA.qadd_order_key(ok) ? oka : ok
    @test isequal(QuantumCumulants.canonical_rep(op, ctx; treatments)[1], expected)
    # The non-conjugate-folded config (scaled_key) is reproduced by _treatment_key.
    @test isequal(QuantumCumulants._treatment_key(op, ctx, treatments), ok)
end

@testset "treatments: complete/scale/evaluate record Free/Scaled/Concrete" begin
    ha = NLevelSpace(:atom, 2); hf = FockSpace(:cavity); h = ha ⊗ hf
    a = Destroy(h, :a)
    σ(α, β, i) = IndexedOperator(Transition(h, :σ, α, β), i)
    @variables N::Real Δ::Real g::Real κ::Real
    i = Index(h, :i, N, ha)
    H = Δ * a' * a + g * ∑(σ(2, 1, i) * a + σ(1, 2, i) * a', i)
    eqs = meanfield(a' * a, H, [a]; rates = [κ], order = 2)
    complete!(eqs)
    atom_sp = first(QuantumCumulants.build_ctx(eqs).symmetric)
    @test get(eqs.graph.treatments, atom_sp, QuantumCumulants.Free) == QuantumCumulants.Free
    @test get(scale(eqs).graph.treatments, atom_sp, QuantumCumulants.Free) == QuantumCumulants.Scaled
    @test get(evaluate(eqs; limits = (N => 3)).graph.treatments, atom_sp, QuantumCumulants.Free) == QuantumCumulants.Concrete

    eqs_scaled_bang = QuantumCumulants._copy(eqs)
    scale!(eqs_scaled_bang)
    @test eqs_scaled_bang.graph.treatments == scale(eqs).graph.treatments
end
