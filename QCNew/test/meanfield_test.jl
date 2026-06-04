using QCNew
using Symbolics: Symbolics, @variables, expand
using SymbolicUtils
using Test

function _iz(x)
    x isa Number && return iszero(x)
    x isa SymbolicUtils.BasicSymbolic || return iszero(x)
    return SymbolicUtils._iszero(SymbolicUtils.simplify(x; expand = true))
end

@testset "meanfield: damped cavity" begin
    hc = FockSpace(:cavity)
    @qnumbers a::Destroy(hc)
    @variables ω κ

    H = ω * a' * a
    eqs = meanfield([a], H, [a]; rates = [κ])

    @test eqs isa MeanFieldEquations
    @test length(eqs.equations) == 1
    @test length(eqs.states) == 1
    @test isequal(eqs.states[1], eqs.equations[1].lhs)
    @test isequal(eqs.operators[1], a * 1)
    @test eqs.hamiltonian === H
    @test isequal(eqs.rates, [κ])
    @test eqs.order === nothing
    @test occursin("⟨a⟩", repr(eqs.equations[1].rhs))
end

@testset "meanfield with order=2" begin
    hc = FockSpace(:cavity)
    @qnumbers a::Destroy(hc)
    @variables ω κ U
    H = ω * a' * a + U * a' * a' * a * a
    eqs = meanfield([a], H, [a]; rates = [κ], order = 2)
    @test eqs.order == [2]
    @test get_order(eqs.equations[1].rhs) <= 2
end

@testset "meanfield: closed two-level" begin
    ha = NLevelSpace(:atom, 2)
    σ(i, j) = Transition(ha, :σ, i, j)
    @variables Δ γ

    H = Δ * σ(2, 2)
    eqs = meanfield([σ(2, 1), σ(2, 2)], H, [σ(1, 2)]; rates = [γ])

    @test length(eqs.equations) == 2
    @test isequal(eqs.states[1], eqs.equations[1].lhs)
    @test isequal(eqs.states[2], eqs.equations[2].lhs)
end

@testset "meanfield: Jaynes-Cummings commutators (parameters)" begin
    hf = FockSpace(:cavity)
    ha = NLevelSpace(:atom, (:g, :e))
    h = hf ⊗ ha
    a = Destroy(h, :a)
    σ = Transition(h, :σ, :g, :e)
    σee = Transition(h, :σ, :e, :e, 2)
    σgg = Transition(h, :σ, :g, :g, 2)

    @variables ωc::Real ωa::Real g::Real
    H = ωc * a' * a + ωa * σ' * σ + g * (a' * σ + σ' * a)

    da = commutator(1im * H, a)
    @test iszero(simplify(da - ((-1im * g) * σ + (-1im * ωc) * a)))

    # `commutator` does not apply the level-completeness rewrite σ_gg + σ_ee = 1,
    # so [iH, σ] is the unreduced form `(-iωa)σ + (-ig)a·σ_gg + (ig)a·σ_ee`.
    ds = commutator(1im * H, σ)
    @test iszero(
        simplify(
            ds - (
                (-1im * g) * a * σgg + (-1im * ωa) * σ +
                    (1im * g) * a * σee
            )
        )
    )

    dn = commutator(1im * H, a' * a)
    @test iszero(simplify(dn - ((-1im * g) * a' * σ + (1im * g) * a * σ')))
end

@testset "meanfield: Jaynes-Cummings commutators (numeric)" begin
    hf = FockSpace(:cavity)
    ha = NLevelSpace(:atom, (:g, :e))
    h = hf ⊗ ha
    a = Destroy(h, :a)
    σ = Transition(h, :σ, :g, :e)
    σee = Transition(h, :σ, :e, :e, 2)
    σgg = Transition(h, :σ, :g, :g, 2)

    (ωc, ωa, g) = (1.1341, 0.4321, 2.15013)
    H = ωc * a' * a + ωa * σ' * σ + g * (a' * σ + σ' * a)

    da = commutator(im * H, a)
    @test iszero(simplify(da - (-im * ωc * a + (-im * g) * σ)))

    ds = commutator(im * H, σ)
    @test iszero(
        simplify(
            ds - (
                (-im * g) * a * σgg + (-im * ωa) * σ +
                    (im * g) * a * σee
            )
        )
    )

    dn = commutator(im * H, a' * a)
    @test iszero(simplify(dn - ((-im * g) * a' * σ + (im * g) * a * σ')))

    # `meanfield` applies SQA's `expand_completeness`, so analytic targets that
    # reference σ_gg explicitly must do the same to land in the same basis before
    # subtracting.
    he = meanfield([a, σ, a' * a], H)
    @test _iz(he.equations[1].rhs - average(expand_completeness(da)))
    @test _iz(he.equations[2].rhs - average(expand_completeness(ds)))
    @test _iz(he.equations[3].rhs - average(expand_completeness(dn)))
end

@testset "meanfield: Lossy JC drift matches analytic" begin
    hf = FockSpace(:cavity)
    ha = NLevelSpace(:atom, (:g, :e))
    h = hf ⊗ ha
    a = Destroy(h, :a)
    σ = Transition(h, :σ, :g, :e)
    σee = Transition(h, :σ, :e, :e, 2)
    σgg = Transition(h, :σ, :g, :g, 2)

    (ωc, ωa, g) = (1.1341, 0.4321, 2.15013)
    κ, γ = 3.333, 0.1313131313
    H = ωc * a' * a + ωa * σ' * σ + g * (a' * σ + σ' * a)
    J = [a, σ]
    he_diss = meanfield([a, σ, σ' * σ], H, J; rates = [κ, γ])

    @test _iz(
        he_diss.equations[1].rhs -
            average((-im * ωc - 0.5κ) * a + (-im * g) * σ)
    )
    # `meanfield` applies SQA's `expand_completeness`, so σ_gg gets
    # rewritten as `1 - σ_ee` on every RHS. Match that basis here.
    @test _iz(
        he_diss.equations[2].rhs -
            average(
            expand_completeness(
                (-im * g) * a * σgg +
                    (-im * ωa - 0.5γ) * σ +
                    (im * g) * a * σee
            )
        )
    )
    @test _iz(
        he_diss.equations[3].rhs -
            average((-γ) * σee + (im * g) * a' * σ + (-im * g) * a * σ')
    )
end

@testset "meanfield: single-atom laser drift matches analytic" begin
    hf = FockSpace(:cavity)
    ha = NLevelSpace(:atom, (:g, :e))
    h = hf ⊗ ha
    a = Destroy(h, :a)
    σ = Transition(h, :σ, :g, :e)
    σee = Transition(h, :σ, :e, :e, 2)
    σgg = Transition(h, :σ, :g, :g, 2)

    (ωc, ωa, g) = (1.1341, 0.4321, 2.15013)
    κ, γ, ν = 3.333, 0.1313131313, 3.44444444
    H = ωc * a' * a + ωa * σ' * σ + g * (a' * σ + σ' * a)
    J = [a, σ, σ']
    he_laser = meanfield([a' * a, σ' * σ, a' * σ], H, J; rates = [κ, γ, ν])

    @test _iz(
        he_laser.equations[1].rhs -
            average((-κ) * a' * a + (-im * g) * a' * σ + (im * g) * a * σ')
    )
    # `meanfield` applies SQA's `expand_completeness`, so σ_gg gets
    # rewritten as `1 - σ_ee` and the analytic form ν*σ_gg - γ*σ_ee
    # collapses to ν + (-ν - γ)*σ_ee. Use the same expansion here.
    @test _iz(
        he_laser.equations[2].rhs -
            average(
            expand_completeness(
                ν * σgg - γ * σee +
                    (im * g) * a' * σ + (-im * g) * a * σ'
            )
        )
    )
end

@testset "meanfield: Jaynes-Cummings system shape" begin
    hf = FockSpace(:cavity)
    ha = NLevelSpace(:atom, (:g, :e))
    h = hf ⊗ ha

    a = Destroy(h, :a)
    σ = Transition(h, :σ, :g, :e)

    @variables ωc ωa g κ γ ν
    H = ωc * a' * a + ωa * σ' * σ + g * (a' * σ + σ' * a)

    he = meanfield([a, σ, a' * a], H)
    @test length(he.equations) == 3
    @test eltype(he.equations) == Symbolics.Equation

    he_diss = meanfield([a, σ, σ' * σ], H, [a, σ]; rates = [κ, γ])
    @test length(he_diss.equations) == 3

    he_laser = meanfield([a' * a, σ' * σ, a' * σ], H, [a, σ, σ']; rates = [κ, γ, ν])
    @test length(he_laser.equations) == 3
end

@testset "meanfield: numbers in commutator and Hamiltonian" begin
    @variables Δ
    h2 = NLevelSpace(:a, 2)
    s(i, j) = Transition(h2, :s, i, j)
    H = Δ * s(2, 2) - Δ * s(1, 1)
    H_s = simplify(H)
    @test meanfield(s(1, 2), H_s) isa MeanFieldEquations
end

@testset "meanfield: Pauli closed two-spin" begin
    hs1 = PauliSpace(:Spin1)
    hs2 = PauliSpace(:Spin2)
    h = hs1 ⊗ hs2
    σ(i, axis) = Pauli(h, Symbol(:σ_, i), axis, i)

    @variables J Δ_1 Δ_2
    H = Δ_1 * σ(1, 3) + Δ_2 * σ(2, 3) + J * σ(1, 1) * σ(2, 1)

    eqs = meanfield([σ(1, 3)], H)
    eqs_c = complete(eqs; order = 2)
    @test length(eqs_c.equations) == 6

    ops2 = [σ(1, 1), σ(1, 2), σ(1, 3), σ(2, 1), σ(2, 2), σ(2, 3)]
    eqs2 = meanfield(ops2, H)
    eqs2_c = complete(eqs2; order = 2)
    @test length(eqs2_c.equations) == 14
end

@testset "meanfield: collective decay rate matrix" begin
    N = 2
    @variables G11::Real G12::Real G21::Real G22::Real δ::Real
    h = ⊗([NLevelSpace(Symbol(:atom, i), 2) for i in 1:N]...)
    σ_(i, j, k) = Transition(h, Symbol("σ_{$k}"), i, j, k)
    H = δ * (σ_(2, 2, 1) + σ_(2, 2, 2))
    ops_0 = [σ_(1, 2, 1)]

    # Native matrix-rates form: `J` is a `Vector{Vector{Op}}`, and `rates`
    # is a `Vector{Matrix}` carrying the cross-rate matrix per bath.
    J_nested = [[σ_(1, 2, 1), σ_(1, 2, 2)]]
    rates_matrix = [[G11 G12; G21 G22]]
    eqs_native = meanfield(ops_0, H, J_nested; rates = rates_matrix)
    eqs_native_c = complete(eqs_native)
    @test eqs_native_c isa MeanFieldEquations
    @test isempty(find_missing(eqs_native_c))

    # Flattened-jumps form should produce the same algebraic content.
    JumpOp = Transition[]
    JumpOpConj = Transition[]
    for i in 1:N, j in 1:N
        push!(JumpOp, σ_(1, 2, i))
        push!(JumpOpConj, σ_(2, 1, j))
    end
    rates_flat = [G11, G21, G12, G22]
    eqs_flat = meanfield(ops_0, H, JumpOp; Jdagger = JumpOpConj, rates = rates_flat)
    eqs_flat_c = complete(eqs_flat)
    @test length(eqs_native_c.equations) == length(eqs_flat_c.equations)

    # Diagonal-only rate case via simple `rates = [γ, γ]` form must also
    # close.
    @variables γ::Real
    eqs2 = meanfield(ops_0, H, [σ_(1, 2, 1), σ_(1, 2, 2)]; rates = [γ, γ])
    eqs_c2 = complete(eqs2)
    @test eqs_c2 isa MeanFieldEquations
end
