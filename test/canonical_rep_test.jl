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
    # All-Free coordinate reproduces canon_key exactly.
    free = QuantumCumulants.all_free_coords(ctx)
    @test isequal(QuantumCumulants.canonical_rep(op, ctx; coords = free)[1], QuantumCumulants.canon_key(op, ctx))
    # idempotence
    rep1 = QuantumCumulants.canonical_rep(op, ctx; coords = free)[1]
    @test isequal(QuantumCumulants.canonical_rep(rep1 * 1, ctx; coords = free)[1], rep1)
    # adjoint-consistency: rep is invariant under {O, O†}, sign bit flips
    r, s = QuantumCumulants.canonical_rep(op, ctx; coords = free)
    ra, sa = QuantumCumulants.canonical_rep(adjoint(op), ctx; coords = free)
    @test isequal(r, ra)
    @test s != sa
end

@testset "canonical_rep: Scaled coordinate equals orbit_key" begin
    ha = NLevelSpace(:atom, 2); hf = FockSpace(:cavity); h = ha ⊗ hf
    a = Destroy(h, :a)
    σ(α, β, i) = IndexedOperator(Transition(h, :σ, α, β), i)
    @variables N::Real g::Real
    i = Index(h, :i, N, ha); j = Index(h, :j, N, ha)
    H = g * ∑(σ(2, 1, i) * a + σ(1, 2, i) * a', i)
    ctx = QuantumCumulants.build_ctx([a], H, SQA.QField[], SQA.QField[])
    op = (σ(2, 2, i) * σ(2, 1, j)) * 1
    atom_sp = first(ctx.symmetric)
    coords = QuantumCumulants.with_coord(QuantumCumulants.all_free_coords(ctx), atom_sp, QuantumCumulants.Scaled)
    # canonical_rep ADDS the conjugation fold on top of the orbit_key node key,
    # so it equals the conjugation-min of orbit_key over {op, op†}.
    ok = QuantumCumulants.orbit_key(op, ctx)
    oka = QuantumCumulants.orbit_key(adjoint(op), ctx)
    expected = QuantumCumulants._serialize(oka) < QuantumCumulants._serialize(ok) ? oka : ok
    @test isequal(QuantumCumulants.canonical_rep(op, ctx; coords)[1], expected)
    # The non-conjugate-folded config (orbit_key) is reproduced by _coord_key.
    @test isequal(QuantumCumulants._coord_key(op, ctx, coords), ok)
end

@testset "coords: complete/scale/evaluate record Free/Scaled/Concrete" begin
    ha = NLevelSpace(:atom, 2); hf = FockSpace(:cavity); h = ha ⊗ hf
    a = Destroy(h, :a)
    σ(α, β, i) = IndexedOperator(Transition(h, :σ, α, β), i)
    @variables N::Real Δ::Real g::Real κ::Real
    i = Index(h, :i, N, ha)
    H = Δ * a' * a + g * ∑(σ(2, 1, i) * a + σ(1, 2, i) * a', i)
    eqs = meanfield(a' * a, H, [a]; rates = [κ], order = 2)
    complete!(eqs)
    atom_sp = first(QuantumCumulants.build_ctx(eqs.operators, eqs.hamiltonian, eqs.jumps, eqs.jumps_dagger).symmetric)
    @test get(eqs.coords, atom_sp, Int(QuantumCumulants.Free)) == Int(QuantumCumulants.Free)
    @test get(scale(eqs).coords, atom_sp, Int(QuantumCumulants.Free)) == Int(QuantumCumulants.Scaled)
    @test get(evaluate(eqs; limits = (N => 3)).coords, atom_sp, Int(QuantumCumulants.Free)) == Int(QuantumCumulants.Concrete)
end
