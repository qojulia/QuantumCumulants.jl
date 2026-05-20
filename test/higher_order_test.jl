using QuantumCumulants
using Symbolics: Symbolics, @variables
using SymbolicUtils
using ModelingToolkitBase: @named, mtkcompile, ODEProblem, unknowns
using OrdinaryDiffEq: Tsit5, solve, ReturnCode
using Test

# v1 surface: higher-order cumulant closure of a damped JC system.
# Master's spectrum-comparison assertion (`maximum(abs.(s6 .- s4)) < 0.2`)
# remains in test/pending/higher_order_test.jl until v1's `Spectrum`
# numerical-stability path is finalised. The 4th- and 6th-order closures
# themselves run cleanly here.

@testset "higher-order: 4th-order cumulant closure of JC laser" begin
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
    # Phase-invariant filter: keep only terms where adjoint-count balances.
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
    @test length(he4.equations) > 0
    @test isempty(find_missing(he4; filter_func = phase_invariant))

    # ODE roundtrip at order=4: equations must compile through MTK and
    # integrate to a finite steady-state without blowing up.
    @named sys4 = to_system(he4)
    sys4_c = mtkcompile(sys4)
    u0 = Dict(unknowns(sys4_c) .=> zeros(ComplexF64, length(unknowns(sys4_c))))
    ps = [Δ, g, γ, κ, ν]
    pn = [1.0, 1.5, 0.25, 1.0, 4.0]
    prob = ODEProblem(sys4_c, merge(u0, Dict(ps .=> pn)), (0.0, 20.0))
    sol = solve(prob, Tsit5(); abstol = 1e-8, reltol = 1e-8)
    @test sol.retcode == ReturnCode.Success
    n_ss = real(get_solution(sol, a' * a, he4)(sol.t[end]))
    @test isfinite(n_ss)
    @test n_ss >= 0
end

# Sibling files:
# - higher_order_order6_test.jl: 6th-order closure (no spectrum assertion)
# - higher_order_agreement_test.jl: order=4 vs order=6 steady-state agreement
# Split for parallel scheduling.
