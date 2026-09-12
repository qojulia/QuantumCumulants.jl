using QuantumCumulants
using OrdinaryDiffEqLowOrderRK: RK4, solve
using OrdinaryDiffEqRosenbrock: OrdinaryDiffEqRosenbrock, Rodas5P
using Random: MersenneTwister
using SciMLBase: ODEFunction, ODEProblem, ReturnCode, remake
using SymbolicUtils: SymbolicUtils
using Symbolics: Symbolics, @variables
using Test

const QC = QuantumCumulants

struct CountingJacobian{J}
    jac::J
    calls::Base.RefValue{Int}
end

function (counter::CountingJacobian)(J, u, p, t)
    counter.calls[] += 1
    return counter.jac(J, u, p, t)
end

function jacobian_fixture()
    n = 3
    h = ⊗([PauliSpace(Symbol(:jac_spin, i)) for i in 1:n]...)
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
    complete!(eqs; get_adjoints = true)
    return eqs, Dict(J => 1.0, hx => 1.0, γ => 0.2), J, hx, γ
end

function jacobian_substitutions(eqs, ps, u)
    subs = Dict{Any, Any}(Symbolics.unwrap(k) => v for (k, v) in ps)
    for (i, state) in enumerate(eqs.states)
        subs[Symbolics.unwrap(state)] = u[i]
    end
    for eq in eqs.equations, var in Symbolics.get_variables(Symbolics.unwrap(eq.rhs))
        v = Symbolics.unwrap(var)
        if SymbolicUtils.issym(v) &&
                Base.nameof(v) === :im &&
                SymbolicUtils.symtype(v) === Number
            subs[v] = im
        end
    end
    return subs
end

@testset "analytic Jacobian agrees with Symbolics derivatives" begin
    eqs, ps, _, _, _ = jacobian_fixture()
    n = length(eqs.states)
    u0 = zeros(ComplexF64, n)
    u = ComplexF64[0.1cos(3.7i) + 0.05im * sin(1.3i) for i in 1:n]
    prob = ODEProblem(
        eqs,
        u0,
        (0.0, 1.0),
        ps;
        backend = KernelBackend(),
        jac = true,
    )

    Jmat = copy(prob.f.jac_prototype)
    prob.f.jac(Jmat, u, prob.p, 0.0)
    subs = jacobian_substitutions(eqs, ps, u)
    rng = MersenneTwister(1)
    maxrel = 0.0
    for _ in 1:30
        i = rand(rng, 1:n)
        j = rand(rng, 1:n)
        derivative = Symbolics.derivative(
            eqs.equations[i].rhs,
            Symbolics.wrap(Symbolics.unwrap(eqs.states[j])),
        )
        ref = ComplexF64(
            SymbolicUtils.unwrap_const(
                Symbolics.substitute(Symbolics.unwrap(derivative), subs),
            ),
        )
        maxrel = max(maxrel, abs(Jmat[i, j] - ref) / max(abs(ref), 1.0e-12))
    end
    @test maxrel < 1.0e-10
end

@testset "Jacobian plan does not alter the RHS execution plan" begin
    eqs, ps, _, _, _ = jacobian_fixture()
    u0 = zeros(ComplexF64, length(eqs.states))
    plain = ODEProblem(eqs, u0, (0.0, 1.0), ps; backend = KernelBackend())
    with_jac = ODEProblem(
        eqs,
        u0,
        (0.0, 1.0),
        ps;
        backend = KernelBackend(),
        jac = true,
    )

    a = plain.f.f.kernel
    b = with_jac.f.f.kernel
    @test a.parent == b.parent
    @test a.leaf == b.leaf
    @test a.rowptr == b.rowptr
    @test a.monomial == b.monomial
    @test a.coeff_id == b.coeff_id
    @test with_jac.f.f.jacobian isa QC.MomentJacobianKernel{ComplexF64}
end

@testset "Jacobian follows parameter updates" begin
    eqs, ps, J, hx, γ = jacobian_fixture()
    u0 = zeros(ComplexF64, length(eqs.states))
    u = ComplexF64[0.08cos(2.7i) + 0.04im * sin(0.9i) for i in eachindex(eqs.states)]
    changed = Dict(J => 1.7, hx => 0.4, γ => 0.31)

    prob = ODEProblem(
        eqs,
        u0,
        (0.0, 1.0),
        ps;
        backend = KernelBackend(),
        jac = true,
    )
    update_parameters!(prob, changed)
    fresh = ODEProblem(
        eqs,
        u0,
        (0.0, 1.0),
        changed;
        backend = KernelBackend(),
        jac = true,
    )

    Ja = copy(prob.f.jac_prototype)
    Jb = copy(fresh.f.jac_prototype)
    prob.f.jac(Ja, u, prob.p, 0.0)
    fresh.f.jac(Jb, u, fresh.p, 0.0)
    @test Ja == Jb
end

@testset "folded closures reject an invalid complex Jacobian" begin
    h = FockSpace(:jac_folded)
    a = Destroy(h, :a)
    @variables Δ Ω κ U
    H = Δ * a' * a + U * a' * a' * a * a + Ω * (a + a')
    eqs = meanfield([a], H, [a]; rates = [κ], order = 2)
    complete!(eqs; get_adjoints = false)
    ps = Dict(Δ => -1.0, Ω => 1.3, κ => 1.0, U => 0.1)

    @test_throws QC.HolomorphicJacobianError ODEProblem(
        eqs,
        zeros(ComplexF64, length(eqs.states)),
        (0.0, 1.0),
        ps;
        backend = KernelBackend(),
        jac = true,
    )
end

@testset "Jacobian mode is explicit and works with an implicit solver" begin
    eqs, ps, _, _, _ = jacobian_fixture()
    u0 = zeros(ComplexF64, length(eqs.states))
    @test_throws ArgumentError ODEProblem(
        eqs,
        u0,
        (0.0, 1.0),
        ps;
        backend = KernelBackend(),
        jac = :fd,
    )

    prob = ODEProblem(
        eqs,
        u0,
        (0.0, 0.2),
        ps;
        backend = KernelBackend(),
        jac = true,
    )
    jac_calls = Ref(0)
    counted_jac = CountingJacobian(prob.f.jac, jac_calls)
    counted_f = ODEFunction{true}(
        prob.f.f;
        jac = counted_jac,
        jac_prototype = copy(prob.f.jac_prototype),
    )
    counted_prob = remake(prob; f = counted_f)
    sol = solve(
        counted_prob,
        Rodas5P(autodiff = OrdinaryDiffEqRosenbrock.AutoFiniteDiff());
        saveat = 0.1,
    )
    @test sol.retcode == ReturnCode.Success
    @test jac_calls[] > 0

    reference = solve(
        ODEProblem(eqs, u0, (0.0, 0.2), ps; backend = KernelBackend()),
        RK4();
        saveat = 0.1,
        abstol = 1.0e-10,
        reltol = 1.0e-10,
    )
    @test maximum(abs, sol.u[end] .- reference.u[end]) < 1.0e-4
end
