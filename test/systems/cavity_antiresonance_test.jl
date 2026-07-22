using QuantumCumulants
using ModelingToolkitBase
using Symbolics: @variables
using OrdinaryDiffEqTsit5: ODEProblem, solve, Tsit5
using Test

@testset "cavity_antiresonance_indexed" begin
    hc = FockSpace(:cavity); ha = NLevelSpace(:atom, 2); h = hc ⊗ ha
    @variables N Δc η Δa κ
    g_v(idx) = IndexedVariable(:g, idx)
    Γ_v(i, j) = DoubleIndexedVariable(:Γ, i, j)
    Ω_v(i, j) = DoubleIndexedVariable(:Ω, i, j; identical = false)
    i = Index(h, :i, N, ha); j = Index(h, :j, N, ha)
    @qnumbers a::Destroy(h)
    σ(x, y, idx) = IndexedOperator(Transition(h, :σ, x, y), idx)
    Hc = Δc * a' * a + η * (a' + a)
    Ha = Δa * Σ(σ(2, 2, i), i) +
        Σ(Σ(Ω_v(i, j) * σ(2, 1, i) * σ(1, 2, j), j, [i]), i)
    Hi = Σ(g_v(i) * (a' * σ(1, 2, i) + a * σ(2, 1, i)), i)
    H = Hc + Ha + Hi
    J = [a, σ(1, 2, i)]
    rates = [κ, Γ_v(i, j)]
    eqs = meanfield(a, H, J; rates = rates, order = 1)
    complete!(eqs)
    @test length(eqs.equations) == 3
    eqs_ = evaluate(eqs; limits = (N => 2))
    @test length(eqs_.equations) == 5
    sys_c = mtkcompile(System(eqs_; name = :ca_sys))
    @test length(unknowns(sys_c)) == 5
    u0 = zeros(ComplexF64, length(eqs_.equations))
    init = initial_values(eqs_, u0)
    pmap_dict = Dict{Any, Any}()
    for p in ModelingToolkitBase.parameters(sys_c)
        nm = string(p)
        if nm == "Δc"
            pmap_dict[p] = 0.0
        elseif nm == "Δa"
            pmap_dict[p] = 0.0
        elseif nm == "η"
            pmap_dict[p] = 0.1
        elseif nm == "κ"
            pmap_dict[p] = 1.0
        elseif nm == "g"
            pmap_dict[p] = [1.0, 1.0]
        elseif nm == "Γ"
            pmap_dict[p] = [1.0 1.0; 1.0 1.0]
        elseif nm == "Ω"
            pmap_dict[p] = [0.0 0.5; 0.5 0.0]
        end
    end
    prob = ODEProblem(sys_c, merge(init, pmap_dict), (0.0, 30.0))
    sol = solve(prob, Tsit5(); reltol = 1.0e-9, abstol = 1.0e-9)
    @test isapprox(
        abs2(get_solution(sol, a, eqs_).(sol.t[end])),
        0.001984; rtol = 1.0e-4,
    )
end
