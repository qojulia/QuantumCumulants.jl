using QuantumCumulants
using Symbolics: @variables
using ModelingToolkitBase: mtkcompile
using OrdinaryDiffEq: ODEProblem, Tsit5, solve
using Test

@testset "dephasing vs concrete-site: closure responds to the channel set" begin
    # The σ^{22} population drift differs when the indexed diagonal dephasing jump
    # σ^{22}_i is in the channel set versus the ladder-only concrete-site set.
    @variables N Δ κ Γ R ν
    hc = FockSpace(:cavity); ha = NLevelSpace(:atom, 2); h = hc ⊗ ha
    i = Index(h, :i, N, ha)
    a = Destroy(h, :a, 1)
    σ(α, β, idx) = IndexedOperator(Transition(h, :σ, α, β, 2), idx)
    H = -Δ * a' * a + Σ(IndexedVariable(:g, i) * (a' * σ(1, 2, i) + a * σ(2, 1, i)), i)

    c_deph = complete(
        meanfield(
            [a' * a, σ(2, 2, i)], H,
            [a, σ(1, 2, i), σ(2, 1, i), σ(2, 2, i)]; rates = [κ, Γ, R, ν], order = 2
        )
    )
    c_conc = complete(
        meanfield(
            [a' * a, σ(2, 2, i)], H,
            [a, σ(1, 2, i)]; rates = [κ, Γ], order = 2
        )
    )

    @test isempty(find_missing(c_deph))
    @test isempty(find_missing(c_conc))

    # The σ^{22} population equation drift differs between the two channel sets.
    σ22_drift(eqs) = first(
        eq.rhs for eq in eqs.equations
            if occursin("σ", string(eq.lhs)) && occursin("₂₂", string(eq.lhs))
    )
    @test !isequal(σ22_drift(c_deph), σ22_drift(c_conc))
end

@testset "dephasing-channel discriminator: readout-shape invariance" begin
    # The derived closure must be identical whether the readout atom index is
    # written free (`j`) or slot-minted (`j(1)`).
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

@testset "dephasing vs concrete-site: different steady-state physics" begin
    # The dephasing channel set has an incoherent pump (the σ^{21} raising jump,
    # rate R), so ⟨σ^{22}⟩ settles at R/(R+Γ) > 0; the concrete-site set has only
    # decay, so ⟨σ^{22}⟩ → 0. Single atom via evaluate(N=>1).
    @variables N Δ g κ Γ R ν
    hc = FockSpace(:cavity); ha = NLevelSpace(:atom, 2); h = hc ⊗ ha
    i = Index(h, :i, N, ha)
    a = Destroy(h, :a, 1)
    σ(α, β, idx) = IndexedOperator(Transition(h, :σ, α, β, 2), idx)
    H = -Δ * a' * a + g * Σ(a' * σ(1, 2, i) + a * σ(2, 1, i), i)

    function σ22_steady(J, rate_syms, rate_vals)
        eqs_c = complete(meanfield([a' * a, σ(2, 2, i)], H, J; rates = rate_syms, order = 1))
        eqs_ev = evaluate(eqs_c; limits = (N => 1))
        sys = mtkcompile(System(eqs_ev; name = :s))
        pmap = parameter_map(eqs_ev, [Δ => 0.0, g => 1.0, rate_vals...])
        u0 = initial_values(eqs_ev, ComplexF64[0.0 for _ in eqs_ev.states])
        prob = ODEProblem(sys, merge(u0, pmap), (0.0, 80.0))
        sol = solve(prob, Tsit5(); reltol = 1.0e-9, abstol = 1.0e-11, saveat = 80.0)
        σ22_lhs = first(eq.lhs for eq in eqs_ev.equations if startswith(string(eq.lhs), "⟨σ"))
        return real(get_solution(sol, undo_average(σ22_lhs), eqs_ev)(80.0))
    end

    Γv, Rv = 0.25, 4.0
    deph_ss = σ22_steady(
        [a, σ(1, 2, i), σ(2, 1, i), σ(2, 2, i)],
        [κ, Γ, R, ν], [κ => 1.0, Γ => Γv, R => Rv, ν => 1.0]
    )
    conc_ss = σ22_steady([a, σ(1, 2, i)], [κ, Γ], [κ => 1.0, Γ => Γv])

    # Dephasing/pumped set: ⟨σ22⟩ ≈ R/(R+Γ); concrete set: ⟨σ22⟩ ≈ 0.
    @test isapprox(deph_ss, Rv / (Rv + Γv); atol = 1.0e-2)
    @test isapprox(conc_ss, 0.0; atol = 1.0e-6)
    @test !isapprox(deph_ss, conc_ss; atol = 1.0e-2)
end
