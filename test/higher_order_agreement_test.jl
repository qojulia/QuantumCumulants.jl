using QuantumCumulants
using Symbolics: Symbolics, @variables
using SymbolicUtils
using ModelingToolkitBase: @named, mtkcompile, ODEProblem, unknowns
using OrdinaryDiffEq: Tsit5, solve, ReturnCode
using Test

# Steady-state agreement between 4th- and 6th-order cumulant closures of
# the damped JC laser. Split out of higher_order_test.jl because the two
# MTK compiles plus two long ODE solves dominate the file's wall-clock;
# isolating it lets ParallelTestRunner schedule it independently of the
# closure-only testsets.

@testset "higher-order: order=4 vs order=6 steady-state agreement" begin
    @variables Δ::Real g::Real γ::Real κ::Real ν::Real
    hf = FockSpace(:cavity)
    ha = NLevelSpace(:atom, 2)
    h = hf ⊗ ha
    a = Destroy(h, :a)
    s = Transition(h, :σ, 1, 2)

    H = Δ * a' * a + g * (a' * s + a * s')
    J = [a, s, s']
    rates = [κ, γ, ν]

    SQA = QuantumCumulants.SecondQuantizedAlgebra
    ϕ(::Number) = 0
    ϕ(::SQA.Destroy) = -1
    ϕ(::SQA.Create) = 1
    function ϕ(t::SQA.Transition)
        t.i == t.j && return 0
        return t.i == 2 ? 1 : -1
    end
    ϕ(q::SQA.QAdd) = isempty(q.arguments) ? 0 : sum(ϕ(term) for (term, _) in q.arguments)
    ϕ(t::SQA.QTerm) = isempty(t.ops) ? 0 : sum(ϕ(op) for op in t.ops)
    function ϕ(avg)
        avg isa SymbolicUtils.BasicSymbolic && SQA.is_average(avg) || return 0
        return ϕ(SQA.undo_average(avg))
    end
    phase_invariant(x) = iszero(ϕ(x))

    he4 = complete(meanfield(a' * a, H, J; rates = rates);
                    order = 4, filter_func = phase_invariant)
    he6 = complete(meanfield(a' * a, H, J; rates = rates);
                    order = 6, filter_func = phase_invariant)

    @named sys4 = to_system(he4); sys4_c = mtkcompile(sys4)
    @named sys6 = to_system(he6); sys6_c = mtkcompile(sys6)
    u04 = Dict(unknowns(sys4_c) .=> zeros(ComplexF64, length(unknowns(sys4_c))))
    u06 = Dict(unknowns(sys6_c) .=> zeros(ComplexF64, length(unknowns(sys6_c))))
    ps = [Δ, g, γ, κ, ν]
    pn = [1.0, 1.5, 0.25, 1.0, 4.0]
    prob4 = ODEProblem(sys4_c, merge(u04, Dict(ps .=> pn)), (0.0, 50.0))
    prob6 = ODEProblem(sys6_c, merge(u06, Dict(ps .=> pn)), (0.0, 50.0))
    sol4 = solve(prob4, Tsit5(); abstol = 1e-9, reltol = 1e-9)
    sol6 = solve(prob6, Tsit5(); abstol = 1e-9, reltol = 1e-9)
    @test sol4.retcode == ReturnCode.Success
    @test sol6.retcode == ReturnCode.Success

    n_ss_4 = real(get_solution(sol4, a' * a, he4)(sol4.t[end]))
    n_ss_6 = real(get_solution(sol6, a' * a, he6)(sol6.t[end]))
    @test isfinite(n_ss_4)
    @test isfinite(n_ss_6)
    @test n_ss_4 >= 0
    @test n_ss_6 >= 0
    @test abs(n_ss_4 - n_ss_6) / max(abs(n_ss_6), 1e-3) < 0.05

    assert_real(sol4, a' * a, he4; atol = 1e-6)
    assert_nonneg(sol4, a' * a, he4; atol = 1e-6)
    assert_real(sol6, a' * a, he6; atol = 1e-6)
    assert_nonneg(sol6, a' * a, he6; atol = 1e-6)
end
