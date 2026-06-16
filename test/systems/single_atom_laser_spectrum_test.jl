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

@testset "single-atom-laser-spectrum" begin
    @variables Δ g γ κ ν
    hf = FockSpace(:cavity); ha = NLevelSpace(:atom, (:g, :e))
    h = hf ⊗ ha
    a = Destroy(h, :a); s = Transition(h, :σ, :g, :e)
    H = Δ * a' * a + g * (a' * s + a * s')
    eq_n = meanfield(a' * a, H, [a, s, s']; rates = [κ, γ, ν], order = 2)
    eqs = complete(eq_n; filter_func = phase_invariant)
    @test length(eqs.equations) == 4
    sys_c = mtkcompile(System(eqs; name = :sal_sys))
    @test length(unknowns(sys_c)) == 4
    u0 = zeros(ComplexF64, length(eqs.equations))
    init = initial_values(eqs, u0)
    ps = Dict(Δ => 0.0, g => 1.5, γ => 0.25, κ => 1.0, ν => 4.0)
    prob = ODEProblem(sys_c, merge(init, ps), (0.0, 20.0))
    sol = solve(prob, Tsit5(); reltol = 1.0e-10, abstol = 1.0e-10)
    @test isapprox(
        real(get_solution(sol, a' * a, eqs).(sol.t[end])),
        1.4835417514686384; rtol = 1.0e-5,
    )
    c = CorrelationFunction(
        a', a, eqs; steady_state = true, filter_func = phase_invariant,
    )
    @test length(c.eqs.equations) == 2
    ps_tup = (Δ, g, γ, κ, ν)
    p0 = (0.0, 1.5, 0.25, 1.0, 4.0)
    sys_steady = mtkcompile(System(eqs; name = :sal_steady))
    u0_ss = zeros(ComplexF64, length(unknowns(sys_steady)))
    dict_ss = merge(Dict(unknowns(sys_steady) .=> u0_ss), Dict(ps_tup .=> p0))
    prob_ss = ODEProblem(sys_steady, dict_ss, (0.0, 20.0))
    sol_ss = solve(prob_ss, Tsit5(); reltol = 1.0e-10, abstol = 1.0e-10)
    S = Spectrum(c, ps_tup)
    s_vals = S([-0.5, 0.0, 0.5], sol_ss.u[end], p0)
    @test isapprox(s_vals[2], 12.01916023864571; rtol = 1.0e-4)
    @test isapprox(s_vals[1], 2.683110742979243; rtol = 1.0e-4)
    # spectrum symmetry S(ω) = S(-ω)
    @test isapprox(s_vals[1], s_vals[3]; rtol = 1.0e-10)
end
