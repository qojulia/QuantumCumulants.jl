using QuantumCumulants
using SymbolicUtils: SymbolicUtils
using ModelingToolkitBase: @named, mtkcompile, ODEProblem, unknowns
using OrdinaryDiffEq: Tsit5, solve, ReturnCode
using Test

@testset "find_missing on JC" begin
    hc = FockSpace(:cavity); ha = NLevelSpace(:atom, 2)
    h = hc ⊗ ha
    @qnumbers a::Destroy(h, 1)
    σ(i, j) = Transition(h, :σ, i, j, 2)
    @variables Δ g κ γ
    H = Δ * a' * a + g * (a * σ(2, 1) + a' * σ(1, 2))
    eqs = meanfield([a, σ(2, 2)], H, [a, σ(1, 2)]; rates = [κ, γ], order = 2)
    missing_states = find_missing(eqs)
    @test !isempty(missing_states)
    @test all(
        SymbolicUtils.iscall(m) && QuantumCumulants.SecondQuantizedAlgebra.is_average(m)
            for m in missing_states
    )
end

@testset "complete! closes JC at order 2" begin
    hc = FockSpace(:cavity); ha = NLevelSpace(:atom, 2)
    h = hc ⊗ ha
    @qnumbers a::Destroy(h, 1)
    σ(i, j) = Transition(h, :σ, i, j, 2)
    @variables Δ g κ γ
    H = Δ * a' * a + g * (a * σ(2, 1) + a' * σ(1, 2))
    eqs = meanfield([a, σ(2, 2)], H, [a, σ(1, 2)]; rates = [κ, γ], order = 2)
    complete!(eqs)
    @test isempty(find_missing(eqs))
    @test length(eqs.equations) >= 2
end

@testset "complete is non-mutating" begin
    hc = FockSpace(:cavity)
    @qnumbers a::Destroy(hc)
    @variables ω κ U
    H = ω * a' * a + U * a' * a' * a * a
    eqs = meanfield([a], H, [a]; rates = [κ], order = 2)
    n_before = length(eqs.equations)
    eqs2 = complete(eqs)
    @test length(eqs.equations) == n_before
    @test isempty(find_missing(eqs2))
end

@testset "find_missing: closure size is invariant of user index naming" begin
    # The closure size depends only on the algebra and the cumulant order, not on
    # the symbol the user picked for the atom index, so `k` and `l` must agree.
    hc = FockSpace(:cavity); ha = NLevelSpace(:atom, 2); h = hc ⊗ ha
    @qnumbers a::Destroy(h, 1)
    σ(α, β, idx) = IndexedOperator(Transition(h, :σ, α, β, 2), idx)
    @variables N::Real Δ::Real g::Real κ::Real Γ::Real R::Real ν::Real

    function closure_size(idx_name::Symbol)
        idx = Index(h, idx_name, N, ha)
        H = -Δ * a' * a + g * (Σ(a' * σ(1, 2, idx), idx) + Σ(a * σ(2, 1, idx), idx))
        J = [a, σ(1, 2, idx), σ(2, 1, idx), σ(2, 2, idx)]
        eqs = complete(
            meanfield(
                [a' * a, σ(2, 2, idx)], H, J;
                rates = [κ, Γ, R, ν], order = 1
            )
        )
        return length(eqs.equations)
    end

    @test closure_size(:k) == closure_size(:l)
end

@testset "complete: conjugate observable is the adjoint state" begin
    # Driven cavity: ⟨a⟩ closes to a single state; ⟨a'⟩ closes to its adjoint.
    hc = FockSpace(:cavity)
    @qnumbers a::Destroy(hc)
    @variables ω κ η
    H = ω * a' * a + 1im * η * (a' - a)

    eqs_a = complete(meanfield(a, H, [a]; rates = [κ]))
    eqs_ad = complete(meanfield(a', H, [a]; rates = [κ]))

    @test length(eqs_a.states) == 1
    @test length(eqs_ad.states) == 1
    a_op = undo_average(only(eqs_a.states))
    ad_op = undo_average(only(eqs_ad.states))
    @test isequal(ad_op, adjoint(a_op))
end

@testset "find_missing get_adjoints=false: one rep per conjugate pair" begin
    # With `get_adjoints=false`, only one representative per conjugate pair is
    # emitted; the conjugate partner is rewritten at codegen time.
    hc = FockSpace(:resonator); ha = NLevelSpace(:atom, 2); h = hc ⊗ ha
    @variables N::Real ωc::Real ωa::Real g::Real κ::Real γ::Real
    j = Index(h, :j, N, ha); k = Index(h, :k, N, ha)
    @qnumbers a::Destroy(h, 1)
    σ(α, β, idx) = IndexedOperator(Transition(h, :σ, α, β, 2), idx)
    H = ωc * a' * a + ωa * Σ(σ(2, 2, j), j) +
        g * a' * Σ(σ(1, 2, j), j) + g * a * Σ(σ(2, 1, j), j)
    J = [a, σ(1, 2, j), σ(2, 1, j)]
    eqs = meanfield(
        [a', a' * a, σ(2, 2, k), σ(1, 2, k), a * a], H, J;
        rates = [κ, γ, γ], order = 2,
    )
    scaled_full = scale(complete(eqs; get_adjoints = true))
    scaled = scale(complete(eqs; get_adjoints = false))
    # `get_adjoints=false` closure is strictly smaller (conjugate partners absorbed).
    @test length(scaled.states) < length(scaled_full.states)
    @named sys = System(scaled)
    @test mtkcompile(sys) isa Any
end

@testset "complete get_adjoints true/false: identical physics" begin
    # Folding ⟨X†⟩ onto conj(⟨X⟩) at codegen must leave every observable
    # unchanged. Driven Jaynes-Cummings, order 2.
    hc = FockSpace(:cavity); ha = NLevelSpace(:atom, 2); h = hc ⊗ ha
    @qnumbers a::Destroy(h, 1)
    σ(α, β) = Transition(h, :σ, α, β, 2)
    @variables Δ g κ γ η
    H = -Δ * a' * a + η * (a' + a) + g * (a' * σ(1, 2) + a * σ(2, 1))
    base = meanfield([a' * a, σ(2, 2)], H, [a, σ(1, 2)]; rates = [κ, γ], order = 2)
    eqs_t = complete(base; get_adjoints = true)
    eqs_f = complete(base; get_adjoints = false)
    @test length(eqs_f.states) < length(eqs_t.states)

    function _solve(eqs)
        sys = mtkcompile(System(eqs; name = :s))
        u0 = Dict(unknowns(sys) .=> zeros(ComplexF64, length(unknowns(sys))))
        p = Dict(Δ => 0.0, g => 1.5, κ => 0.7, γ => 0.3, η => 1.2)
        return solve(
            ODEProblem(sys, merge(u0, p), (0.0, 3.0)), Tsit5();
            reltol = 1.0e-10, abstol = 1.0e-12
        )
    end
    sol_t = _solve(eqs_t); sol_f = _solve(eqs_f)
    @test sol_t.retcode == ReturnCode.Success
    @test sol_f.retcode == ReturnCode.Success
    for τ in (0.5, 1.0, 2.0, 3.0)
        @test isapprox(
            get_solution(sol_t, a' * a, eqs_t)(τ),
            get_solution(sol_f, a' * a, eqs_f)(τ); atol = 1.0e-8
        )
        @test isapprox(
            get_solution(sol_t, a, eqs_t)(τ),
            get_solution(sol_f, a, eqs_f)(τ); atol = 1.0e-8
        )
    end
end

@testset "complete: closes when every user index is bound" begin
    # Every declared index is consumed by an H sum / jump; the system must
    # still close.
    @variables κ::Real ω::Real N::Real M::Real
    hc = FockSpace(:cavity); hf = FockSpace(:filter)
    ha = NLevelSpace(:atom, 2)
    h = hc ⊗ hf ⊗ ha
    @qnumbers a::Destroy(h, 1)
    b(k) = IndexedOperator(Destroy(h, :b, 2), k)
    i_idx = Index(h, :i, M, hf)
    H = ω * Σ(b(i_idx)' * b(i_idx), i_idx)
    J = [a, b(i_idx)]
    eqs = complete(meanfield([a' * a], H, J; rates = [κ, κ], order = 1))
    @test isempty(find_missing(eqs))
    @test length(eqs.equations) >= 1
end

@testset "complete: irrelevant-NE leaves do not over-count the closure" begin
    # `complete` must be idempotent on its own output: a second pass adds no
    # new states.
    @variables N::Real Δ::Real g::Real κ::Real Γ::Real
    hc = FockSpace(:cavity); ha = NLevelSpace(:atom, 2); h = hc ⊗ ha
    @qnumbers a::Destroy(h, 1)
    σ(α, β, k) = IndexedOperator(Transition(h, :σ, α, β, 2), k)
    i = Index(h, :i, N, ha)
    H = Δ * Σ(σ(2, 2, i), i) + g * Σ(a' * σ(1, 2, i) + a * σ(2, 1, i), i)
    J = [a, σ(1, 2, i), σ(2, 1, i), σ(2, 2, i)]
    eqs_c = complete(meanfield([a' * a, σ(2, 2, i)], H, J; rates = [κ, Γ, Γ, Γ], order = 2))
    n1 = length(eqs_c.equations)
    eqs_cc = complete(eqs_c)
    @test length(eqs_cc.equations) == n1
    @test isempty(find_missing(eqs_c))
end

@testset "complete!: free LHS atom index does not collide with H-bound i" begin
    # When H is summed over `i` and the user-declared free LHS index reuses the
    # name `i`, the cross-atom commutator must be preserved, so the closure grows
    # past the single-atom case.
    @variables Δ::Real g::Real κ::Real Γ::Real N::Real
    hc = FockSpace(:cavity); ha = NLevelSpace(:atom, 2); h = hc ⊗ ha
    @qnumbers a::Destroy(h, 1)
    σ(α, β, k) = IndexedOperator(Transition(h, :σ, α, β, 2), k)
    i_idx = Index(h, :i, N, ha)

    H = Δ * Σ(σ(2, 2, i_idx), i_idx) +
        g * Σ(a' * σ(1, 2, i_idx) + a * σ(2, 1, i_idx), i_idx)
    J = [a, σ(1, 2, i_idx)]
    eqs = meanfield([a' * a, σ(2, 2, i_idx)], H, J; rates = [κ, Γ], order = 2)
    eqs_c = complete(eqs)

    # The order-2 indexed JC laser closes at >= 7 equations (atom-atom cross
    # states present).
    @test length(eqs_c.equations) >= 7
end
