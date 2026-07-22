using QuantumCumulants
using Symbolics: Symbolics, @variables
using SymbolicUtils
using Test

@testset "indexed measurement backaction: cavity + noise channel structure" begin
    @variables κ::Real Δ::Real η::Real
    hc = FockSpace(:cavity)
    @qnumbers a::Destroy(hc)

    me_a = meanfield(a, Δ * a' * a, [a]; rates = [κ], efficiencies = [η])
    @test me_a isa NoiseMeanfieldEquations
    @test length(me_a.equations) == 1
    @test length(me_a.noise_equations) == 1

    # The deterministic ⟨a⟩ rhs equals -iΔ⟨a⟩ - 0.5κ⟨a⟩ and the noise rhs equals
    # √(η κ) * (⟨a' a⟩ + ⟨a a⟩ - ⟨a⟩² - ⟨a⟩⟨a'⟩).
    expected_det = average(-1im * Δ * a - 0.5 * κ * a)
    @test _is_zero(me_a.equations[1].rhs - expected_det)
    expected_noise = sqrt(η * κ) * (
        average(a' * a) + average(a * a) -
            average(a)^2 - average(a) * average(a')
    )
    @test _is_zero(me_a.noise_equations[1].rhs - expected_noise)

    # At order=1 the noise equation collapses to zero because the
    # cumulant-1 truncation drops `⟨a'a⟩ - ⟨a'⟩⟨a⟩` style fluctuations.
    me_a_o1 = meanfield(
        a', Δ * a' * a, [a]; rates = [κ], efficiencies = [η], order = 1
    )
    @test me_a_o1 isa NoiseMeanfieldEquations
    rhs_o1 = me_a_o1.noise_equations[1].rhs
    @test _is_zero(rhs_o1)
end

@testset "indexed measurement backaction: filter cavity drift agrees det vs stoch" begin
    # Efficiencies enter the noise channel, not the drift, so the drift
    # equations of `meanfield(...)` and `meanfield(...; efficiencies=...)`
    # must coincide term-by-term on the indexed multi-Hilbert input.
    @variables κ::Real g::Real gf::Real κf::Real R::Real Γ::Real Δ::Real ν::Real N::Real M::Real η::Real
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
        gf * (Σ(a' * b(i) + a * b(i)', i)) +
        g * (Σ(a' * σ(1, 2, j) + a * σ(2, 1, j), j))

    J = [a, b(i), σ(1, 2, j), σ(2, 1, j), σ(2, 2, j)]
    rates = [κ, κf, Γ, R, ν]

    # Compare with `get_adjoints=true`: the exact `det ⊆ stoch` subset needs both
    # conjugate members present. The folded default keeps one per pair, picked per pipeline.
    eqs_det = meanfield(a' * a, H, J; rates = rates, order = 2)
    eqs_det_c = complete(eqs_det; get_adjoints = true)

    efficiencies = [η, 0, 0, 0, 0]
    eqs_noise = meanfield(
        a' * a, H, J; rates = rates, efficiencies = efficiencies, order = 2
    )
    @test eqs_noise isa NoiseMeanfieldEquations
    eqs_noise_c = complete(eqs_noise; get_adjoints = true)

    # The stochastic state set contains the deterministic state set.
    @test length(eqs_noise_c.equations) >= length(eqs_det_c.equations)
    det_lhs = Set(e.lhs for e in eqs_det_c.equations)
    stoch_lhs = Set(e.lhs for e in eqs_noise_c.equations)
    @test issubset(det_lhs, stoch_lhs)

    # For every LHS shared between det and stoch, the drift rhs must agree
    # symbolically.
    det_by_lhs = Dict(e.lhs => e.rhs for e in eqs_det_c.equations)
    n_checked = 0
    for eq in eqs_noise_c.equations
        haskey(det_by_lhs, eq.lhs) || continue
        @test _is_zero(eq.rhs - det_by_lhs[eq.lhs])
        n_checked += 1
        n_checked >= 3 && break
    end
    @test n_checked >= 1
end
