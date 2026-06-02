using QuantumCumulants
using SecondQuantizedAlgebra
const SQA = SecondQuantizedAlgebra
using Symbolics: @variables
using Test

# The cross-atom NE policy is selected by a STRUCTURAL property of the channel
# set, `_has_dephasing_channel` (does any indexed diagonal atomic jump `σ^{αα}`
# appear?), replacing the old `user_concretes` naming heuristic. These tests
# lock the two properties that motivated that change (see TODO §3):
#   1. the classification is correct and depends only on the jump operators;
#   2. the resulting closure is invariant to how the user writes a readout atom
#      index (free `j` vs slot-minted `j(1)`), the exact case the old heuristic
#      mis-classified.

const QC = QuantumCumulants

@testset "dephasing-channel discriminator: classification" begin
    @variables N
    hc = FockSpace(:cavity)
    ha = NLevelSpace(:atom, 2)
    h = hc ⊗ ha
    i = Index(h, :i, N, ha)
    a = Destroy(h, :a, 1)
    σ(α, β, idx) = IndexedOperator(Transition(h, :σ, α, β, 2), idx)

    # `σ^{22}_i` is an indexed diagonal (dephasing) jump → population system.
    @test QC._has_dephasing_channel([a, σ(1, 2, i), σ(2, 1, i), σ(2, 2, i)])
    @test QC._has_dephasing_channel([σ(2, 2, i)])
    # Only ladder/decay atomic jumps and a bosonic loss → concrete-site system.
    @test !QC._has_dephasing_channel([a, σ(1, 2, i)])
    @test !QC._has_dephasing_channel([a])
    # A NON-indexed diagonal transition is not a population-defining channel
    # (a scalar system has no cross-atom moments to fold).
    @test !QC._has_dephasing_channel([Transition(h, :σ, 2, 2, 2)])

    # Per-operator predicate.
    @test QC._is_diag_atom_jump(σ(2, 2, i))
    @test !QC._is_diag_atom_jump(σ(1, 2, i))      # off-diagonal (decay)
    @test !QC._is_diag_atom_jump(a)               # bosonic
    @test !QC._is_diag_atom_jump(Transition(h, :σ, 2, 2, 2))  # not indexed
end

@testset "dephasing-channel discriminator: readout-shape invariance" begin
    # Population laser (σ^{22} dephasing present). The derived closure must be
    # identical whether the readout atom index is written free (`j`) or
    # slot-minted (`j(1)`); the discriminator keys on the channel set, not on
    # the readout index's syntactic shape.
    @variables N Δ κ Γ R ν
    hc = FockSpace(:cavity)
    ha = NLevelSpace(:atom, 2)
    h = hc ⊗ ha
    g(idx) = IndexedVariable(:g, idx)
    i = Index(h, :i, N, ha)
    j = Index(h, :j, N, ha)
    a = Destroy(h, :a, 1)
    σ(α, β, idx) = IndexedOperator(Transition(h, :σ, α, β, 2), idx)

    H = -Δ * a' * a + Σ(g(i) * (a' * σ(1, 2, i) + a * σ(2, 1, i)), i)
    J = [a, σ(1, 2, i), σ(2, 1, i), σ(2, 2, i)]
    rates = [κ, Γ, R, ν]

    readout_free = σ(2, 2, j)
    readout_slot = IndexedOperator(Transition(h, :σ, 2, 2, 2), j(1))

    c_free = complete(meanfield([a' * a, readout_free], H, J; rates = rates, order = 2))
    c_slot = complete(meanfield([a' * a, readout_slot], H, J; rates = rates, order = 2))

    @test length(c_free.equations) == length(c_slot.equations)
    @test length(scale(c_free).equations) == length(scale(c_slot).equations)
end
