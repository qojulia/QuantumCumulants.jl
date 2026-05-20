using QuantumCumulants
using Symbolics: Symbolics, @variables
using SymbolicUtils
using ModelingToolkitBase: ModelingToolkitBase, @named, mtkcompile, ODEProblem, unknowns
using OrdinaryDiffEq: Tsit5, solve, ReturnCode
using Test

# v1 surface: mixed-order cumulant on an indexed JC laser. Master's expected
# closure size (8 equations after complete! with get_adjoints=false; 18 after
# evaluate at N=3) is reproduced once cumulant_expansion under a per-Hilbert
# order vector recurses into sub-blocks whose subspace has lower per-space
# order (see src/cumulant.jl).

@testset "mixed-order: indexed 1-atom-cavity collective closes under order=[1,2]" begin
    ha = NLevelSpace(:atom, 2)
    hf = FockSpace(:cavity)
    h = ha ⊗ hf

    a = Destroy(h, :a)
    σ(α, β, i) = IndexedOperator(Transition(h, :σ, α, β), i)

    @variables η::Real Γ::Real κ::Real Δc::Real N::Real ξ::Real ν::Real
    g(i) = IndexedVariable(:g, i)
    Δa(i) = IndexedVariable(:Δa, i)

    i = Index(h, :i, N, ha)
    j = Index(h, :j, N, ha)

    H = -Δc * a' * a - ∑(Δa(i) * σ(2, 2, i), i) +
        ∑(g(i) * (σ(2, 1, i) * a + σ(1, 2, i) * a'), i) +
        1im * η * (a' - a)

    J = [σ(1, 2, i), a, a' * a, σ(2, 2, i)]
    rates = [Γ, κ, ξ, ν]

    eqs = meanfield(a * σ(2, 2, j), H, J; rates = rates, order = [1, 2])
    complete!(eqs; get_adjoints = false)
    @test eqs.order == [1, 2]
    @test length(eqs.equations) == 8
    @test isempty(find_missing(eqs; get_adjoints = false))

    SQA = QuantumCumulants.SecondQuantizedAlgebra
    for state in eqs.states
        @test SQA.is_average(state)
    end
end

@testset "mixed-order: evaluate on order=[1,2] with N=3 closes" begin
    ha = NLevelSpace(:atom, 2)
    hf = FockSpace(:cavity)
    h = ha ⊗ hf

    a = Destroy(h, :a)
    σ(α, β, i) = IndexedOperator(Transition(h, :σ, α, β), i)

    @variables η::Real Γ::Real κ::Real Δc::Real N::Real ξ::Real ν::Real
    g(i) = IndexedVariable(:g, i)
    Δa(i) = IndexedVariable(:Δa, i)

    i = Index(h, :i, N, ha)
    j = Index(h, :j, N, ha)

    H = -Δc * a' * a - ∑(Δa(i) * σ(2, 2, i), i) +
        ∑(g(i) * (σ(2, 1, i) * a + σ(1, 2, i) * a'), i) +
        1im * η * (a' - a)

    J = [σ(1, 2, i), a, a' * a, σ(2, 2, i)]
    rates = [Γ, κ, ξ, ν]

    eqs = meanfield(a * σ(2, 2, j), H, J; rates = rates, order = [1, 2])
    complete!(eqs; get_adjoints = false)
    @test length(eqs.equations) == 8
    evaled = evaluate(eqs; limits = (N => 3))
    @test length(evaled.equations) == 18
    @test isempty(find_missing(evaled; get_adjoints = false))
    @test evaled.order == [1, 2]

    SQA = QuantumCumulants.SecondQuantizedAlgebra
    for state in evaled.states
        op = SQA.undo_average(state)
        len = if op isa SQA.QAdd && length(op.arguments) == 1
            length(first(op.arguments)[1].ops)
        else
            1
        end
        @test len <= 2
    end
end

@testset "mixed-order: ODE integration, photon-number physicality (N=2)" begin
    # Strengthening: master's count assertion (`length(eqs) == 8`,
    # `length(evaled) == 18`) is sensitive to v1's cumulant expansion
    # tagging (a valid alternative closure derives more equations). What
    # we CAN verify across both branches: the resulting ODE integrates,
    # the steady-state cavity number is finite and non-negative, and the
    # trajectory is physical (real photon number, non-negative). Uses the
    # per-atom inhomogeneous g(i), Δa(i) IndexedVariables from master.
    ha = NLevelSpace(:atom, 2)
    hf = FockSpace(:cavity)
    h = ha ⊗ hf

    a = Destroy(h, :a)
    σ(α, β, i) = IndexedOperator(Transition(h, :σ, α, β), i)

    @variables η::Real Γ::Real κ::Real Δc::Real N::Real ν::Real
    g(i) = IndexedVariable(:g, i)
    Δa(i) = IndexedVariable(:Δa, i)

    i = Index(h, :i, N, ha)

    H = -Δc * a' * a - ∑(Δa(i) * σ(2, 2, i), i) +
        ∑(g(i) * (σ(2, 1, i) * a + σ(1, 2, i) * a'), i) +
        1im * η * (a' - a)
    J = [σ(1, 2, i), a, σ(2, 2, i)]
    rates = [Γ, κ, ν]

    eqs = meanfield(a' * a, H, J; rates = rates, order = [1, 2])
    complete!(eqs)
    evaled = evaluate(eqs; limits = (N => 2))
    @test isempty(find_missing(evaled; get_adjoints = false))

    SQA = QuantumCumulants.SecondQuantizedAlgebra
    @named sys = System(evaled)
    sys_c = mtkcompile(sys)
    # Below-threshold parameter point. Above threshold (strong η, large ν)
    # the cumulant closure is not bounded; this is a physics constraint,
    # not a numerical one.
    pmap = parameter_map(
        evaled, Dict(
            g(i) => 0.1,
            Δa(i) => 0.0,
            Δc => 0.0,
            κ => 1.0,
            Γ => 0.25,
            ν => 0.0,
            η => 0.5,
        )
    )
    u0 = Dict(unknowns(sys_c) .=> zeros(ComplexF64, length(unknowns(sys_c))))
    prob = ODEProblem(sys_c, merge(u0, pmap), (0.0, 20.0))
    sol = solve(prob, Tsit5(); abstol = 1.0e-8, reltol = 1.0e-8)
    @test sol.retcode == ReturnCode.Success

    a_op = SQA.undo_average(evaled.states[1])
    assert_real(sol, a_op, evaled; atol = 1.0e-6)
    assert_nonneg(sol, a_op, evaled; atol = 1.0e-6)
end
