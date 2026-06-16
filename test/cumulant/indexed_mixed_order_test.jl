using QuantumCumulants
using Symbolics: Symbolics, @variables
using SymbolicUtils
using ModelingToolkitBase: ModelingToolkitBase, @named, mtkcompile, ODEProblem, unknowns
using OrdinaryDiffEqTsit5: Tsit5, solve, ReturnCode
using Test

# Mixed-order cumulant on an indexed JC laser: a per-Hilbert order vector
# closes at 8 equations after complete! (get_adjoints=false), 18 after
# evaluate at N=3.

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
    @test isempty(find_missing(eqs))

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
    @test isempty(find_missing(evaled))
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
    # The resulting ODE integrates and the cavity photon number stays real and
    # non-negative along the trajectory, with per-atom inhomogeneous g(i),
    # Δa(i).
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
    @test isempty(find_missing(evaled))

    SQA = QuantumCumulants.SecondQuantizedAlgebra
    @named sys = System(evaled)
    sys_c = mtkcompile(sys)
    # Below-threshold parameter point; above threshold the cumulant closure is
    # not bounded.
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
