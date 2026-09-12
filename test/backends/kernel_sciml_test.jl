using QuantumCumulants
using ModelingToolkitBase: mtkcompile
using OrdinaryDiffEqLowOrderRK: RK4, solve
using SciMLBase: ODEProblem, remake
using Symbolics: @variables
using Test

const QC = QuantumCumulants

function sciml_pauli_fixture()
    n = 3
    h = ⊗([PauliSpace(Symbol(:sciml_spin, i)) for i in 1:n]...)
    σx(i) = Pauli(h, :σ, 1, i)
    σy(i) = Pauli(h, :σ, 2, i)
    σz(i) = Pauli(h, :σ, 3, i)
    σm(i) = (σx(i) - 1im * σy(i)) / 2
    @variables J hx γ
    H = -J * sum(σz(i) * σz(i + 1) for i in 1:(n - 1)) - hx * sum(σx(i) for i in 1:n)
    eqs = meanfield(
        [σz(i) for i in 1:n],
        H,
        [σm(i) for i in 1:n];
        rates = fill(γ, n),
        order = 2,
    )
    complete!(eqs)
    return eqs, J, hx, γ
end

@testset "KernelBackend ODEProblem has a concrete numeric parameter payload" begin
    eqs, J, hx, γ = sciml_pauli_fixture()
    ps = Dict(J => 1.0, hx => 1.0, γ => 0.2)
    u0 = zeros(ComplexF64, length(eqs.states))
    prob = ODEProblem(eqs, u0, (0.0, 1.0), ps; backend = KernelBackend())

    @test prob.p isa QC.KernelParameters{ComplexF64}
    @test prob.p.values isa Vector{ComplexF64}
    @test prob.p.coeffs isa Vector{ComplexF64}
    @test prob.f.f isa QC.KernelRHS
    @test prob.f.f.kernel isa QC.MomentKernel{ComplexF64}
    @test !hasfield(typeof(prob.f.f.kernel), :coeffs)
end

@testset "parameter updates match fresh construction without relowering" begin
    eqs, J, hx, γ = sciml_pauli_fixture()
    ps = Dict(J => 1.0, hx => 1.0, γ => 0.2)
    ps2 = Dict(J => 1.7, hx => 0.4, γ => 0.31)
    u0 = zeros(ComplexF64, length(eqs.states))
    u = ComplexF64[0.1cos(3.7i) + 0.05im * sin(1.3i) for i in eachindex(eqs.states)]

    prob = ODEProblem(eqs, u0, (0.0, 1.0), ps; backend = KernelBackend())
    kernel = prob.f.f.kernel
    update_parameters!(prob, ps2)
    fresh = ODEProblem(eqs, u0, (0.0, 1.0), ps2; backend = KernelBackend())

    @test prob.f.f.kernel === kernel
    du_updated = similar(u)
    du_fresh = similar(u)
    prob.f(du_updated, u, prob.p, 0.0)
    fresh.f(du_fresh, u, fresh.p, 0.0)
    @test du_updated == du_fresh
end

@testset "remake parameter state is isolated" begin
    eqs, J, hx, γ = sciml_pauli_fixture()
    ps = Dict(J => 1.0, hx => 1.0, γ => 0.2)
    u0 = zeros(ComplexF64, length(eqs.states))
    u = ComplexF64[0.1cos(2.1i) + 0.03im * sin(0.8i) for i in eachindex(eqs.states)]

    prob = ODEProblem(eqs, u0, (0.0, 1.0), ps; backend = KernelBackend())
    prob2 = remake(prob; p = copy(prob.p))
    update_parameters!(prob2, Dict(J => 2.0))

    du1 = similar(u)
    du2 = similar(u)
    prob.f(du1, u, prob.p, 0.0)
    prob2.f(du2, u, prob2.p, 0.0)
    @test du1 != du2

    fresh = ODEProblem(
        eqs,
        u0,
        (0.0, 1.0),
        Dict(J => 2.0, hx => 1.0, γ => 0.2);
        backend = KernelBackend(),
    )
    duf = similar(u)
    fresh.f(duf, u, fresh.p, 0.0)
    @test du2 == duf
end

@testset "direct trajectory agrees with ModelingToolkit" begin
    eqs, J, hx, γ = sciml_pauli_fixture()
    ps = Dict(J => 1.0, hx => 0.7, γ => 0.2)
    u0 = ComplexF64[0.2cos(i) + 0.05im * sin(0.3i) for i in eachindex(eqs.states)]
    tspan = (0.0, 1.0)
    saveat = 0.1

    direct = solve(
        ODEProblem(eqs, u0, tspan, ps; backend = KernelBackend()),
        RK4();
        saveat,
        abstol = 1.0e-10,
        reltol = 1.0e-10,
    )

    sys = mtkcompile(System(eqs; name = :kernel_reference))
    values = merge(initial_values(eqs, u0), Dict(parameter_map(eqs, ps)))
    mtk_prob = ODEProblem(sys, parameter_map(sys, values), tspan)
    reference = solve(mtk_prob, RK4(); saveat, abstol = 1.0e-10, reltol = 1.0e-10)

    for state in eqs.states
        got = get_solution(direct, state, eqs).(direct.t)
        ref = get_solution(reference, state, eqs).(reference.t)
        @test maximum(abs, got .- ref) < 1.0e-6
    end
end

@testset "direct solution lookup recovers folded conjugates" begin
    h = FockSpace(:solution_cavity)
    a = Destroy(h, :a)
    @variables Δ Ω κ U
    H = Δ * a' * a + U * a' * a' * a * a + Ω * (a + a')
    eqs = meanfield([a], H, [a]; rates = [κ], order = 2)
    complete!(eqs; get_adjoints = false)
    ps = Dict(Δ => -1.0, Ω => 1.3, κ => 1.0, U => 0.1)

    sol = solve(
        ODEProblem(
            eqs,
            zeros(ComplexF64, length(eqs.states)),
            (0.0, 0.5),
            ps;
            backend = KernelBackend(),
        ),
        RK4();
        saveat = 0.1,
    )
    a_traj = get_solution(sol, a, eqs).(sol.t)
    adag_traj = get_solution(sol, a', eqs).(sol.t)
    @test adag_traj == conj.(a_traj)
end

@testset "u0 vector and state-keyed dictionary agree" begin
    eqs, J, hx, γ = sciml_pauli_fixture()
    ps = Dict(J => 1.0, hx => 1.0, γ => 0.2)
    u0 = ComplexF64[0.1i for i in eachindex(eqs.states)]
    dict0 = Dict(eqs.states[i] => u0[i] for i in eachindex(u0))

    pvec = ODEProblem(eqs, u0, (0.0, 1.0), ps; backend = KernelBackend())
    pdict = ODEProblem(eqs, dict0, (0.0, 1.0), ps; backend = KernelBackend())
    @test pvec.u0 == pdict.u0
end
