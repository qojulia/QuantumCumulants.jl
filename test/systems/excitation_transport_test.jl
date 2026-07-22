using QuantumCumulants
using ModelingToolkitBase
using Symbolics: @variables
using OrdinaryDiffEqTsit5: ODEProblem, solve, Tsit5
using Test

@testset "excitation-transport-chain (N=10 structural)" begin
    N = 10
    h_chain = ⊗([NLevelSpace(Symbol(:atom, i), (:g, :e)) for i in 1:N]...)
    σ_ch(i, j, k) = Transition(h_chain, Symbol(:σ_, k), i, j, k)
    @variables Ω γ Δ J0
    x_ = [first(@variables $(Symbol("x_$i"))) for i in 1:N]
    Jc(xi, xj) = J0 / abs(xi - xj)^3
    H = -Δ * sum(σ_ch(:e, :e, k) for k in 1:N) +
        Ω * (σ_ch(:e, :g, 1) + σ_ch(:g, :e, 1)) +
        sum(
        Jc(x_[k], x_[k + 1]) * (
                σ_ch(:e, :g, k) * σ_ch(:g, :e, k + 1) +
                σ_ch(:g, :e, k) * σ_ch(:e, :g, k + 1)
            )
            for k in 1:(N - 1)
    )
    c_ops = [σ_ch(:g, :e, k) for k in 1:N]
    eqs = meanfield(σ_ch(:g, :e, 1), H, c_ops; rates = [γ for _ in 1:N], order = 2)
    complete!(eqs)
    # Order-2 closure: 245 conjugate-folded moments (`get_adjoints=false`). complete!'s
    # iteration cap must reach the fixpoint (the N=10 chain needs > 200 BFS
    # node-expansions), else closure truncates instead of erroring.
    @test isempty(find_missing(eqs))
    @test length(eqs.equations) == 245
    sys_c = mtkcompile(System(eqs; name = :chain_sys))
    @test length(unknowns(sys_c)) == 245
end

# N=4 chain integrated as a mean-field ODE, pinning the end-of-chain population.
@testset "excitation-transport-chain (N=4 numerical)" begin
    N = 4
    h_chain = ⊗([NLevelSpace(Symbol(:atom, i), (:g, :e)) for i in 1:N]...)
    σ_ch(i, j, k) = Transition(h_chain, Symbol(:σ_, k), i, j, k)
    @variables Ω γ Δ J0
    x_ = [first(@variables $(Symbol("x4_$i"))) for i in 1:N]
    Jc(xi, xj) = J0 / abs(xi - xj)^3
    H = -Δ * sum(σ_ch(:e, :e, k) for k in 1:N) +
        Ω * (σ_ch(:e, :g, 1) + σ_ch(:g, :e, 1)) +
        sum(
        Jc(x_[k], x_[k + 1]) * (
                σ_ch(:e, :g, k) * σ_ch(:g, :e, k + 1) +
                σ_ch(:g, :e, k) * σ_ch(:e, :g, k + 1)
            )
            for k in 1:(N - 1)
    )
    c_ops = [σ_ch(:g, :e, k) for k in 1:N]
    eqs = meanfield(σ_ch(:g, :e, 1), H, c_ops; rates = [γ for _ in 1:N], order = 2)
    complete!(eqs)
    @test length(eqs.equations) == 38
    sys_c = mtkcompile(System(eqs; name = :chain4_sys))
    @test length(unknowns(sys_c)) == 38
    u0 = zeros(ComplexF64, length(eqs.equations))
    init = initial_values(eqs, u0)
    ps = Dict(
        Ω => 2.0, γ => 1.0, Δ => 0.0, J0 => 1.0,
        x_[1] => 0.0, x_[2] => 1.0, x_[3] => 2.0, x_[4] => 3.0,
    )
    prob = ODEProblem(sys_c, merge(init, ps), (0.0, 15.0))
    sol = solve(prob, Tsit5(); reltol = 1.0e-9, abstol = 1.0e-9)
    @test isapprox(
        real(get_solution(sol, σ_ch(:e, :e, N), eqs).(sol.t[end])),
        0.09034045248945875; rtol = 1.0e-4,
    )
end
