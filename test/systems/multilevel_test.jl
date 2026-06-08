using QuantumCumulants
using Symbolics: Symbolics, @variables
using SymbolicUtils
using Test

@testset "V-level atom: cumulant_expansion + complete to 16-eq closure" begin
    hf = FockSpace(:cavity)
    ha = NLevelSpace(:atom, 3)
    h = hf ⊗ ha
    @variables g::Real Δ_2::Real Δ_3::Real Ω_2::Real Ω_3::Real
    @variables Δ_c::Real Γ_2::Real Γ_3::Real κ::Real
    a = Destroy(h, :a)
    σ(i, j) = Transition(h, :σ, i, j)

    H_atom = -Δ_2 * σ(2, 2) - Δ_3 * σ(3, 3) +
        Ω_2 * (σ(2, 1) + σ(1, 2)) + Ω_3 * (σ(3, 1) + σ(1, 3))
    H_cav = -Δ_c * a' * a + g * (a' * σ(1, 2) + a * σ(2, 1))
    H = H_atom + H_cav

    J = [a, σ(1, 2), σ(1, 3)]
    rates = [κ, Γ_2, Γ_3]

    he_n = meanfield(a' * a, H, J; rates = rates)
    he_avg = cumulant_expansion(he_n, 2)

    he_c_full = complete(he_avg)
    @test isempty(find_missing(he_c_full))

    he_c_canon = complete(he_avg; get_adjoints = false)
    @test isempty(find_missing(he_c_canon; get_adjoints = false))
    @test length(he_c_canon.equations) == 16
    @test length(he_c_full.equations) >= length(he_c_canon.equations)
end

@testset "N-atom two-level laser: complete + filter_func" begin
    N = 3  # smaller N keeps the closed system tractable
    @variables κ
    Δ = [Symbolics.variable(Symbol(:Δ_, i)) for i in 1:N]
    g = [Symbolics.variable(Symbol(:g_, i)) for i in 1:N]
    γ = [Symbolics.variable(Symbol(:γ_, i)) for i in 1:N]
    ν = [Symbolics.variable(Symbol(:ν_, i)) for i in 1:N]

    h_cavity = FockSpace(:cavity)
    h_atoms = [NLevelSpace(Symbol(:atom, i), (:g, :e)) for i in 1:N]
    h = ⊗(h_cavity, h_atoms...)

    a = Destroy(h, :a)
    σ(i, j, k) = Transition(h, Symbol(:σ_, k), i, j, k + 1)

    H = sum(Δ[i] * σ(:e, :e, i) for i in 1:N) +
        sum(g[i] * (a' * σ(:g, :e, i) + a * σ(:e, :g, i)) for i in 1:N)
    J = [a; [σ(:g, :e, k) for k in 1:N]; [σ(:e, :g, k) for k in 1:N]]
    Jdagger = adjoint.(J)
    rates = [κ, γ..., ν...]

    SQA = QuantumCumulants.SecondQuantizedAlgebra
    ϕ(::Destroy) = -1
    ϕ(::Create) = 1
    function ϕ(t::Transition)
        t.i == t.j && return 0
        return t.i == :e ? 1 : -1
    end
    ϕ(q::SQA.QAdd) = sum(ϕ(arg) for (arg, _) in q.arguments)
    ϕ(q::SQA.QTerm) = sum(ϕ(op) for op in q.ops)
    function ϕ(avg)
        avg isa SymbolicUtils.BasicSymbolic && SQA.is_average(avg) || return 0
        return ϕ(SQA.undo_average(avg))
    end
    phase_invariant(x) = iszero(ϕ(x))

    he = meanfield(a' * a, H, J; Jdagger = Jdagger, rates = rates, order = 2)
    # The 2N+1 laser closure is the conjugate-FOLDED basis: ⟨a'σ_k⟩ is counted
    # once per atom (its conjugate ⟨aσ_k⟩ is the same unknown, recovered at
    # codegen), so complete with get_adjoints=false. The ground populations
    # ⟨σ_gg_k⟩ fold to 1 - ⟨σ_ee_k⟩ via the moment-level completeness reduction,
    # so only ⟨σ_ee_k⟩ is tracked.
    complete!(he; filter_func = phase_invariant, get_adjoints = false)
    @test isempty(find_missing(he; filter_func = phase_invariant, get_adjoints = false))
    # Closure basis: ⟨a'a⟩ + ⟨σee_k⟩ + ⟨a'σ_k⟩ gives 2N+1 equations.
    n_eqs = 2N + 1
    @test length(he.equations) == n_eqs

    he_ops = meanfield(
        [
            a' * a; [σ(:e, :e, k) for k in 1:N];
            [a' * σ(:g, :e, k) for k in 1:N]
        ], H, J;
        Jdagger = Jdagger, rates = rates, order = 2
    )
    @test length(he_ops.equations) == n_eqs
end
