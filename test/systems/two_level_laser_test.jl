using QuantumCumulants
using Symbolics: Symbolics, @variables, substitute
using SymbolicUtils
using ModelingToolkitBase: ModelingToolkitBase, @named, mtkcompile, ODEProblem, unknowns
using OrdinaryDiffEq: Tsit5, solve, ReturnCode
using Test

@testset "JC model: cumulant_expansion ≡ meanfield(order=2)" begin
    # Verify that cumulant-expansion at order 2 commutes with the order kwarg
    hf = FockSpace(:cavity)
    ha = NLevelSpace(:atom, (:g, :e))
    h = hf ⊗ ha
    a = Destroy(h, :a)
    σ = Transition(h, :σ, :g, :e)
    @variables Δ::Real g::Real κ::Real γ::Real ν::Real

    H = Δ * a' * a + g * (a' * σ + σ' * a)
    J = [a, σ, σ']
    he = meanfield([a' * a, σ' * σ, a * σ'], H, J; rates = [κ, γ, ν])
    he_exp = cumulant_expansion(he, 2)
    he_o2 = meanfield([a' * a, σ' * σ, a * σ'], H, J; rates = [κ, γ, ν], order = 2)
    @test length(he_exp.equations) == length(he_o2.equations)
    for (e1, e2) in zip(he_exp.equations, he_o2.equations)
        @test isequal(e1.lhs, e2.lhs)
        diff = e1.rhs - e2.rhs
        @test SymbolicUtils._iszero(SymbolicUtils.simplify(diff; expand = true))
    end
end

@testset "Single-atom laser: phase-invariant filter closes system" begin
    hf = FockSpace(:cavity)
    ha = NLevelSpace(:atom, (:g, :e))
    h = hf ⊗ ha
    a = Destroy(h, :a)
    σ(i, j) = Transition(h, :σ, i, j)
    @variables Δ g κ γ ν

    H = Δ * a' * a + g * (a' * σ(:g, :e) + σ(:e, :g) * a)
    J = [a, σ(:g, :e), σ(:e, :g)]

    ϕ(::Destroy) = -1
    ϕ(::Create) = 1
    function ϕ(t::Transition)
        t.i == t.j && return 0
        return t.i == :e ? 1 : -1
    end
    SQA = QuantumCumulants.SecondQuantizedAlgebra
    ϕ(q::SQA.QAdd) = sum(ϕ(arg) for (arg, _) in q.arguments)
    ϕ(q::SQA.QTerm) = sum(ϕ(op) for op in q.ops)
    function ϕ(avg)
        avg isa SymbolicUtils.BasicSymbolic && SQA.is_average(avg) || return 0
        return ϕ(SQA.undo_average(avg))
    end
    phase_invariant(x) = iszero(ϕ(x))

    he = meanfield(a' * a, H, J; rates = [κ, γ, ν], order = 2)
    complete!(he; filter_func = phase_invariant)
    @test isempty(find_missing(he; filter_func = phase_invariant))
    @test length(he.equations) >= 3

    @named sys = System(he)
    @test sys isa ModelingToolkitBase.System
end

@testset "Single-atom laser: phase-invariant ODE + sanity bounds" begin
    hf = FockSpace(:cavity)
    ha = NLevelSpace(:atom, (:g, :e))
    h = hf ⊗ ha
    a = Destroy(h, :a)
    σ = Transition(h, :σ, :g, :e)

    @variables Δ::Real g::Real κ::Real γ::Real ν::Real
    H = Δ * a' * a + g * (a' * σ + σ' * a)
    J = [a, σ, σ']
    he_laser = meanfield([a' * a, σ' * σ, a * σ'], H, J; rates = [κ, γ, ν])
    he_avg = cumulant_expansion(he_laser, 2)
    he_comp = complete(he_avg)

    ps = (Δ, g, γ, κ, ν)
    p0 = ps .=> ComplexF64[1, 1.5, 0.25, 1, 4]
    @named sys = System(he_comp)
    sys_c = mtkcompile(sys)
    u0 = unknowns(sys_c) .=> zeros(ComplexF64, length(unknowns(sys_c)))
    tmax = 10.0
    prob = ODEProblem(sys_c, merge(Dict(u0), Dict(p0)), (0.0, tmax))
    sol = solve(prob, Tsit5())
    @test sol.retcode == ReturnCode.Success

    # Photon population stays real and non-negative; atomic excitation real and in [0, 1].
    assert_real(sol, a' * a, he_comp)
    assert_nonneg(sol, a' * a, he_comp)
    assert_real(sol, σ' * σ, he_comp)
    assert_population(sol, σ' * σ, he_comp)
end
