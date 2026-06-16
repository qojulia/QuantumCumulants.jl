using QuantumCumulants
using ModelingToolkitBase
using Symbolics: @variables
using OrdinaryDiffEqTsit5: ODEProblem, solve, Tsit5
using Test

import QuantumCumulants.SecondQuantizedAlgebra as SQA

_phase_inv(x) = 0
_phase_inv(::Destroy) = -1
_phase_inv(::Create) = 1
function _phase_inv(τ::Transition)
    return (τ.i == 2 && τ.j == 1) ? 1 :
        (τ.i == 1 && τ.j == 2) ? -1 : 0
end
function _phase_inv(q::SQA.QAdd)
    for (term, _) in q.arguments
        p = 0
        for op in term.ops
            p += _phase_inv(op)
        end
        return p
    end
    return 0
end
function _phase_inv(avg)
    SQA.is_average(avg) || return 0
    return _phase_inv(SQA.undo_average(avg))
end
phase_invariant(x) = iszero(_phase_inv(x))

@testset "superradiant_laser_indexed" begin
    hc = FockSpace(:cavity); ha = NLevelSpace(:atom, 2); h = hc ⊗ ha
    @qnumbers a::Destroy(h)
    σ(α, β, idx) = IndexedOperator(Transition(h, :σ, α, β), idx)
    @variables N Δ κ Γ R ν
    g_v(idx) = IndexedVariable(:g, idx)
    i = Index(h, :i, N, ha); j = Index(h, :j, N, ha)
    H = -Δ * a' * a + Σ(g_v(i) * (a' * σ(1, 2, i) + a * σ(2, 1, i)), i)
    J = [a, σ(1, 2, i), σ(2, 1, i), σ(2, 2, i)]
    rates = [κ, Γ, R, ν]
    ops = [a' * a, σ(2, 2, j)]
    eqs = meanfield(ops, H, J; rates = rates, order = 2)
    eqs_c = complete(eqs; filter_func = phase_invariant)
    @test length(eqs_c.equations) == 6
    eqs_sc = scale(eqs_c)
    @test length(eqs_sc.equations) == 5
    sys_c = mtkcompile(System(eqs_sc; name = :sr_sys))
    @test length(unknowns(sys_c)) == 5
    u0 = zeros(ComplexF64, length(eqs_sc.equations))
    init = initial_values(eqs_sc, u0)
    pmap = parameter_map(
        eqs_sc, Dict(
            N => 50.0, Δ => 0.0, g_v(i) => 1.0, κ => 1.0, Γ => 0.25,
            R => 4.0, ν => 0.0,
        )
    )
    prob = ODEProblem(sys_c, merge(init, pmap), (0.0, 30.0))
    sol = solve(prob, Tsit5(); reltol = 1.0e-8, abstol = 1.0e-8)
    @test isapprox(
        real(get_solution(sol, a' * a, eqs_sc).(sol.t[end])),
        91.96090096151146; rtol = 1.0e-4,
    )
    corr = CorrelationFunction(
        a', a, eqs_c; steady_state = true, filter_func = phase_invariant,
    )
    corr_sc = scale(corr)
    @test length(corr_sc.eqs.equations) >= 1
    ps_tup = (Δ, κ, Γ, R, ν, N, g_v(i))
    p0 = (0.0, 1.0, 0.25, 4.0, 0.0, 50.0, 1.0)
    S = Spectrum(corr_sc, ps_tup)
    s_vals = S([-0.5, 0.0, 0.5], sol.u[end], p0)
    @test all(isfinite, s_vals)
    # spectrum symmetry S(ω) = S(-ω)
    @test isapprox(s_vals[1], s_vals[3]; rtol = 1.0e-6)
end
