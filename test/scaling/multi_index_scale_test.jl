using QuantumCumulants
using Symbolics: @variables
using SymbolicUtils: SymbolicUtils
using Test
import QuantumCumulants as QC
import QuantumCumulants.SecondQuantizedAlgebra as SQA

# A scaled subspace with two declared atom indices (`vocab[atom] = [i, j]`) must still
# flatten every coefficient on a minted slot; otherwise an indexed `g` would survive scaling.
@testset "scale with two declared atom indices flattens minted-slot coefficients" begin
    @variables N::Real Δ::Real κ::Real Γ::Real R::Real ν::Real
    hc = FockSpace(:cavity); ha = NLevelSpace(:atom, 2); h = hc ⊗ ha
    @qnumbers a::Destroy(h)
    σ(α, β, idx) = IndexedOperator(Transition(h, :σ, α, β), idx)
    g(idx) = IndexedVariable(:g, idx)
    i = Index(h, :i, N, ha); j = Index(h, :j, N, ha); k = Index(h, :k, N, ha)

    # Coupling sits on the bound index k; the tracked two-body moment carries free i and j.
    H = -Δ * a' * a + Σ(g(k) * (a' * σ(1, 2, k) + a * σ(2, 1, k)), k)
    J = [a, σ(1, 2, k), σ(2, 1, k), σ(2, 2, k)]
    eqs = meanfield([a' * a, σ(2, 1, i) * σ(1, 2, j)], H, J; rates = [κ, Γ, R, ν], order = 2)

    # Precondition: the atom subspace really has two declared reps (else the test is vacuous).
    @test length(eqs.graph.ctx.vocab[2]) == 2

    eqs_c = complete(eqs)
    eqs_sc = scale(eqs_c)
    @test isempty(find_missing(eqs_sc))

    # The Scaled minting path was exercised: minted slots appear in the keys.
    minted = Symbol[SQA.index_name(idx) for k in keys(eqs_sc.graph.nodes) for idx in SQA.get_indices(k)]
    @test :i_1 in minted

    # Core property: no indexed coefficient survives scaling, on any equation.
    any_indexed_var(x) =
        (
        x = SymbolicUtils.unwrap(x); QC._is_indexed_var(x) ||
            (SymbolicUtils.iscall(x) && any(any_indexed_var, SymbolicUtils.arguments(x)))
    )
    @test !any(any_indexed_var(eq.rhs) for eq in eqs_sc.equations)
end

# Regression: two mutually-`≠` collapsed indices must give N(N−1), not (N−1)^2, and a
# symbolic range must not throw.
@testset "_sum_scope_prefactor: falling-factorial count, symbolic range safe" begin
    ha = NLevelSpace(:atom, 2)
    σ(α, β, k) = IndexedOperator(Transition(ha, :σ, α, β), k)
    moment(ri, rj) = SQA.undo_average(
        average(
            Σ(
                Σ(
                    σ(2, 1, Index(ha, :i, ri, ha)) * σ(1, 2, Index(ha, :j, rj, ha)),
                    Index(ha, :j, rj, ha), [Index(ha, :i, ri, ha)]
                ), Index(ha, :i, ri, ha)
            )
        ),
    )
    @test isequal(QC._sum_scope_prefactor(moment(5, 5), Set([1])), 20)   # 5 * 4

    @variables N::Int
    pf = QC._sum_scope_prefactor(moment(N, N), Set([1]))                 # no throw on symbolic N
    @test isequal(Symbolics.simplify(pf - N * (N - 1); expand = true), 0)

    # Single collapsed index with an external `≠` partner keeps the N−1 count.
    ext = Index(ha, :k, N, ha)
    ib = Index(ha, :i, N, ha)
    op1 = SQA.undo_average(average(Σ(σ(2, 1, ext) * σ(1, 2, ib), ib, [ext])))
    @test isequal(Symbolics.simplify(QC._sum_scope_prefactor(op1, Set([1])) - (N - 1); expand = true), 0)
end
