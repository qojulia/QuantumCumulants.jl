using QuantumCumulants
using Symbolics: @variables
using SymbolicUtils: SymbolicUtils
using ModelingToolkitBase: ModelingToolkitBase, brownians, mtkcompile, unknowns
using OrdinaryDiffEqTsit5: ODEProblem, solve, Tsit5
using Test

# Local zero-test for symbolic noise-drift comparisons (mirrors noise_test.jl).
function _iz(x)
    x isa Number && return iszero(x)
    x isa SymbolicUtils.BasicSymbolic || return iszero(x)
    return SymbolicUtils._iszero(SymbolicUtils.simplify(x; expand = true))
end

# Two independent driven-damped cavities, each loss channel monitored. Independent
# subsystems make the multi-channel noise structure unambiguous: monitoring one mode
# may not inject noise into the other, and a linear system keeps the order-2 cumulants
# exact, so the conditional ensemble mean equals the deterministic mean-field exactly.

@testset "multi-channel noise: monitored channels superpose independently" begin
    hc2 = FockSpace(:c1) ⊗ FockSpace(:c2)
    a = Destroy(hc2, :a, 1); b = Destroy(hc2, :b, 2)
    @variables Δa::Real Δb::Real κa::Real κb::Real Ea::Real Eb::Real ηa::Real ηb::Real
    H = Δa * a' * a + Δb * b' * b + Ea * (a + a') + Eb * (b + b')
    ops = [a, b, a' * a, b' * b, a * a, b * b]

    # Same `ops` ⇒ identical state order across the three calls, so index i aligns.
    m12 = meanfield(ops, H, [a, b]; rates = [κa, κb], efficiencies = [ηa, ηb], order = 2)
    ma = meanfield(ops, H, [a, b]; rates = [κa, κb], efficiencies = [ηa, 0], order = 2)
    mb = meanfield(ops, H, [a, b]; rates = [κa, κb], efficiencies = [0, ηb], order = 2)
    @test all(
        isequal(m12.states[i], ma.states[i]) && isequal(m12.states[i], mb.states[i])
            for i in eachindex(m12.states)
    )

    # Independence: the two-channel noise drift is exactly the sum of the single-channel
    # drifts (no cross terms), per equation.
    @test all(
        _iz(
                m12.noise_equations[i].rhs -
                (ma.noise_equations[i].rhs + mb.noise_equations[i].rhs)
            )
            for i in eachindex(m12.noise_equations)
    )

    # Each channel drives noise (neither single-channel drift is trivially zero).
    @test count(!_iz(ma.noise_equations[i].rhs) for i in eachindex(ma.noise_equations)) > 0
    @test count(!_iz(mb.noise_equations[i].rhs) for i in eachindex(mb.noise_equations)) > 0

    # Brownian count tracks the number of active (nonzero-efficiency) channels.
    # `complete` closes the cross-mode moments the noise drift references.
    @test length(brownians(System(complete(m12); name = :two_active))) == 2
    @test length(brownians(System(complete(ma); name = :one_active))) == 1
end

@testset "multi-channel noise: two-channel SDE ensemble matches deterministic" begin
    using StochasticDiffEq: SDEProblem, EM, EnsembleProblem
    using StochasticDiffEq.SciMLBase: ReturnCode
    import Random

    hc2 = FockSpace(:c1) ⊗ FockSpace(:c2)
    a = Destroy(hc2, :a, 1); b = Destroy(hc2, :b, 2)
    @variables Δa::Real Δb::Real κa::Real κb::Real Ea::Real Eb::Real ηa::Real ηb::Real
    H = Δa * a' * a + Δb * b' * b + Ea * (a + a') + Eb * (b + b')
    ops = [a, b, a' * a, b' * b, a * a, b * b]
    eqs = complete(
        meanfield(ops, H, [a, b]; rates = [κa, κb], efficiencies = [ηa, ηb], order = 2);
        get_adjoints = false,
    )

    sys_st = mtkcompile(System(eqs; name = :cav_st))
    sys_det = mtkcompile(System(MeanfieldEquations(eqs); name = :cav_det))

    p0 = Dict(
        Δa => 0.0, Δb => 0.0, κa => 2.0, κb => 3.0,
        Ea => 1.0, Eb => 0.7, ηa => 0.8, ηb => 0.5,
    )
    T_end = 2.0
    tspan = range(0.0, T_end, length = 41)
    u0_det = Dict{Any, Any}(unknowns(sys_det) .=> zeros(ComplexF64, length(unknowns(sys_det))))
    u0_st = Dict{Any, Any}(unknowns(sys_st) .=> zeros(ComplexF64, length(unknowns(sys_st))))
    pd = parameter_map(sys_det, merge(u0_det, Dict{Any, Any}(p0)))
    ps = parameter_map(sys_st, merge(u0_st, Dict{Any, Any}(p0)))

    sol_det = solve(ODEProblem(sys_det, pd, (0.0, T_end)), Tsit5(); reltol = 1.0e-9, abstol = 1.0e-9, saveat = tspan)
    na_det = real.(get_solution(sol_det, a' * a, eqs).(tspan))
    nb_det = real.(get_solution(sol_det, b' * b, eqs).(tspan))

    Random.seed!(7)
    traj = 300
    sol = solve(
        EnsembleProblem(SDEProblem(sys_st, ps, (0.0, T_end))),
        EM(); dt = T_end / 2000, trajectories = traj, saveat = tspan,
    )
    na = zeros(length(tspan)); nb = zeros(length(tspan)); nsucc = 0
    for i in 1:traj
        s = sol.u[i]
        s.retcode == ReturnCode.Success || continue
        na .+= real.(get_solution(s, a' * a, eqs).(tspan))
        nb .+= real.(get_solution(s, b' * b, eqs).(tspan))
        nsucc += 1
    end
    @test nsucc > 0.95 * traj
    na ./= nsucc; nb ./= nsucc

    # Ensemble mean of the conditional trajectories tracks the deterministic
    # mean-field; Monte-Carlo error at 300 trajectories is well below the tolerance.
    @test maximum(abs, na .- na_det) < 0.02
    @test maximum(abs, nb .- nb_det) < 0.02
end
