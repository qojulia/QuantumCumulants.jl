using QuantumCumulants
using Symbolics: Symbolics, @variables
using ModelingToolkitBase: ModelingToolkitBase, @named, mtkcompile, ODEProblem, unknowns
using OrdinaryDiffEq: Tsit5, solve, ReturnCode
using Test

@testset "to_system: damped cavity end-to-end" begin
    hc = FockSpace(:cavity)
    @qnumbers a::Destroy(hc)
    @variables ω κ
    H = ω * a' * a
    eqs = meanfield([a], H, [a]; rates = [κ])
    complete!(eqs)
    sys = to_system(eqs; name = :damped_cav)

    sys_c = mtkcompile(sys)
    @test length(unknowns(sys_c)) == 1

    u0 = initial_values(eqs; defaults = Dict(average(a) => 1.0 + 0.0im))
    ω_val, κ_val = 2.0, 0.5
    p = Dict(ω => ω_val, κ => κ_val)
    prob = ODEProblem(sys_c, merge(u0, p), (0.0, 5.0))
    sol = solve(prob, Tsit5(); reltol = 1.0e-10, abstol = 1.0e-12)
    @test sol.retcode == ReturnCode.Success

    # Verify against analytical solution: a(t) = a(0) * exp((-iω - κ/2)*t).
    # ⟨a'a⟩ is NOT in this closure (the linear Heisenberg equation for ⟨a⟩
    # closes on itself), so the analytical match here is the physicality
    # check; |⟨a⟩|² decreases monotonically as exp(-κ t).
    for τ in (0.0, 1.0, 2.0, 5.0)
        num_val = get_solution(sol, a, eqs)(τ)
        ana_val = exp((-im * ω_val - κ_val / 2) * τ)
        @test abs(num_val - ana_val) < 1.0e-9
    end
    # |⟨a⟩(t)|² decays monotonically (Cauchy-Schwarz-stable lower bound).
    mags = [abs2(get_solution(sol, a, eqs)(t)) for t in sol.t]
    @test all(diff(mags) .<= 1.0e-9)
end
