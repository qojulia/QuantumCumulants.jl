using QCNew
using SecondQuantizedAlgebra: SecondQuantizedAlgebra
using Symbolics: Symbolics, @variables
using SymbolicUtils: SymbolicUtils
using Test

@testset "indexed meanfield: collective atomic emission" begin
    hc = FockSpace(:cavity); ha = NLevelSpace(:atom, 2); h = hc ⊗ ha
    @qnumbers a::Destroy(h, 1)
    σ(x, y, k) = IndexedOperator(Transition(h, :σ, x, y, 2), k)
    @variables N Δ g κ γ
    i = Index(h, :i, N, ha)
    H = Δ * a' * a + Σ(g * (a * σ(2, 1, i) + a' * σ(1, 2, i)), i)
    eqs = meanfield([a], H, [a, σ(1, 2, i)]; rates = [κ, γ])
    @test eqs isa MeanFieldEquations
    # `meanfield` derives one equation per seed operator; one seed (`a`) → 1.
    @test length(eqs.equations) == 1
end

@testset "indexed average appears in RHS" begin
    hc = FockSpace(:cavity); ha = NLevelSpace(:atom, 2); h = hc ⊗ ha
    @qnumbers a::Destroy(h, 1)
    σ(x, y, k) = IndexedOperator(Transition(h, :σ, x, y, 2), k)
    @variables N Δ g κ γ
    i = Index(h, :i, N, ha)
    H = Δ * a' * a + Σ(g * (a * σ(2, 1, i) + a' * σ(1, 2, i)), i)
    eqs = meanfield([a], H, [a, σ(1, 2, i)]; rates = [κ, γ])
    @test occursin("σ", repr(eqs.equations[1].rhs))
end

@testset "indexed meanfield: two-atom collective JC drift" begin
    @variables Δc::Real η::Real Δa::Real κ::Real Γ::Real
    N = 2
    hc = FockSpace(:cavity)
    ha = NLevelSpace(Symbol(:atom), 2)
    h = hc ⊗ ha

    i_ind = Index(h, :i, N, ha)
    k_ind = Index(h, :k, N, ha)

    g(k) = IndexedVariable(:g, k)

    @qnumbers a::Destroy(h)
    σ(i, j, k) = IndexedOperator(Transition(h, :σ, i, j), k)

    Hc = Δc * a' * a + η * (a' + a)
    Ha = Δa * Σ(σ(2, 2, i_ind), i_ind)
    Hi = Σ(g(i_ind) * (a' * σ(1, 2, i_ind) + a * σ(2, 1, i_ind)), i_ind)
    H = Hc + Ha + Hi

    J = [a, σ(1, 2, i_ind)]
    rates = [κ, Γ]

    ops = [a, σ(2, 2, k_ind), σ(1, 2, k_ind)]
    eqs = meanfield(ops, H, J; rates = rates, order = 2)
    # `meanfield` derives one equation per seed; three seeds (`a`, `σ22`,
    # `σ12`) → 3.
    @test length(eqs.equations) == 3
    @test eqs isa MeanFieldEquations
    ceqs = complete(eqs)
    @test ceqs isa MeanFieldEquations
    @test length(ceqs.equations) >= 3

    # `complete` is idempotent on a closed system.
    ceqs_again = complete(ceqs)
    @test length(ceqs_again.equations) == length(ceqs.equations)

    # Every index in the completed system traces back to the user's declared
    # vocabulary, including slot reps `Symbol(base, "_", k)` derived from a
    # declared index.
    declared = Set([i_ind.name, k_ind.name])
    function _is_user_derived(name::Symbol)
        name in declared && return true
        s = String(name)
        for d in declared
            ds = String(d) * "_"
            startswith(s, ds) || continue
            tail = s[(length(ds) + 1):end]
            isempty(tail) && continue
            all(isdigit, tail) && return true
        end
        return false
    end
    for op in ceqs.operators
        for idx in SecondQuantizedAlgebra.get_indices(op)
            @test _is_user_derived(idx.name)
        end
    end
end

@testset "indexed jumps: per-jump dissipator is summed over jump index" begin
    # Indexed jumps induce independent-decay semantics (Σ_i D[J_i] applied
    # separately for each atom), so the i ≠ j contribution vanishes by
    # commutation and no cross-site σ_i*σ_j products survive on the RHS.
    hc = FockSpace(:cavity); ha = NLevelSpace(:atom, 2); h = hc ⊗ ha
    @qnumbers a::Destroy(h)
    σ(α, β, k) = IndexedOperator(Transition(h, :σ, α, β), k)
    @variables N Δ R Γ ν
    g(k) = IndexedVariable(:g, k)
    i = Index(h, :i, N, ha); j = Index(h, :j, N, ha)

    H = -Δ * a' * a + Σ(g(i) * (a' * σ(1, 2, i) + a * σ(2, 1, i)), i)
    J = [σ(1, 2, i), σ(2, 1, i), σ(2, 2, i)]
    rates = [Γ, R, ν]

    eqs = meanfield([σ(2, 2, j)], H, J; rates = rates, order = 2)
    rhs = eqs.equations[1].rhs

    # The dissipator collapses to
    #   R + (-R - Γ) ⟨σ_j^{22}⟩ + i g_j ⟨a' σ_j^{12}⟩ - i g_j ⟨a σ_j^{21}⟩.
    expected = R + (-R - Γ) * average(σ(2, 2, j)) +
        Symbolics.IM * g(j) * average(a' * σ(1, 2, j)) -
        Symbolics.IM * g(j) * average(a * σ(2, 1, j))
    diff = rhs - expected
    @test SymbolicUtils._iszero(SymbolicUtils.simplify(diff; expand = true))
end

@testset "indexed jumps: σ_jj equation has no σ²/σ³ leak after scale" begin
    SQA = SecondQuantizedAlgebra
    hc = FockSpace(:cavity); ha = NLevelSpace(:atom, 2); h = hc ⊗ ha
    @qnumbers a::Destroy(h)
    σ(α, β, k) = IndexedOperator(Transition(h, :σ, α, β), k)
    @variables N Δ κ Γ R ν g
    i = Index(h, :i, N, ha); j = Index(h, :j, N, ha)
    g_idx(k) = IndexedVariable(:g, k)

    H = -Δ * a' * a + Σ(g_idx(i) * (a' * σ(1, 2, i) + a * σ(2, 1, i)), i)
    J = [a, σ(1, 2, i), σ(2, 1, i), σ(2, 2, i)]
    rates = [κ, Γ, R, ν]

    φ(x) = 0
    φ(::Destroy) = -1
    φ(::Create) = 1
    φ(x::SQA.QSym) = x isa SQA.Transition ? x.i - x.j : 0
    function φ(q::SQA.QAdd)
        for (term, _) in q.arguments
            p = 0
            for op in term.ops
                p += φ(op)
            end
            return p
        end
        return 0
    end
    φ(avg) = SQA.is_average(avg) ? φ(SQA.undo_average(avg)) : 0
    phase_invariant(x) = iszero(φ(x))

    eqs = meanfield([a' * a, σ(2, 2, j)], H, J; rates = rates, order = 2)
    eqs_c = complete(eqs; filter_func = phase_invariant)
    eqs_sc = scale(eqs_c)

    # Pick the σ_{j_1}^{22} equation. `scale` renames free indices to suffixed
    # slots; `j` is the only free LHS atom index, so canonical slot 1 is `j(1)`.
    j1 = j(1)
    σ22_op = average(σ(2, 2, j1))
    σ22_idx = findfirst(eq -> isequal(eq.lhs, σ22_op), eqs_sc.equations)
    @test σ22_idx !== nothing
    σ22_rhs = eqs_sc.equations[σ22_idx].rhs

    # `IndexedVariable(:g, i)` flattens to the scalar `g` under scale.
    expected = R + (-R - Γ) * average(σ(2, 2, j1)) +
        Symbolics.IM * g * average(a' * σ(1, 2, j1)) -
        Symbolics.IM * g * average(a * σ(2, 1, j1))
    diff = σ22_rhs - expected
    @test SymbolicUtils._iszero(SymbolicUtils.simplify(diff; expand = true))
end


@testset "indexed meanfield: alpha-equivalent leaves are one state" begin
    # With `H` summed over `i_ind` and `ops` using `k_ind`, the leaf
    # ⟨σ_{i_ind}^{12}⟩ and the user-declared state ⟨σ_{k_ind}^{12}⟩ denote the
    # same observable; `find_missing` must not introduce a duplicate state.
    @variables Δc::Real η::Real Δa::Real κ::Real Γ::Real
    N = 2
    hc = FockSpace(:cavity); ha = NLevelSpace(:atom, 2); h = hc ⊗ ha
    i_ind = Index(h, :i, N, ha); k_ind = Index(h, :k, N, ha)
    @qnumbers a::Destroy(h)
    σ(i, j, k) = IndexedOperator(Transition(h, :σ, i, j), k)
    g(k) = IndexedVariable(:g, k)

    H = Δc * a' * a + η * (a' + a) + Δa * Σ(σ(2, 2, i_ind), i_ind) +
        Σ(g(i_ind) * (a' * σ(1, 2, i_ind) + a * σ(2, 1, i_ind)), i_ind)
    ops = [a, σ(2, 2, k_ind), σ(1, 2, k_ind)]
    J = [a, σ(1, 2, i_ind)]
    eqs = meanfield(ops, H, J; rates = [κ, Γ], order = 2)

    # `find_missing` must not report the i-indexed single-σ states as missing;
    # they are alpha-equivalent to the declared k-indexed ones.
    missing_after_first = QCNew.find_missing(eqs)
    for m in missing_after_first
        s = repr(m)
        @test !occursin(r"⟨σ_i₁₂⟩", s)
        @test !occursin(r"⟨σ_i₂₂⟩", s)
        @test !occursin(r"⟨σ_i₂₁⟩", s)
    end
end
