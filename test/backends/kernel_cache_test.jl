using QuantumCumulants
using Symbolics: Symbolics, @variables
using SciMLBase: SciMLBase, ODEProblem
using JLD2: jldopen
using Test

# Two-spin order-1 model (small fixture; the cache machinery is what is under test).
Ns = 2
hs = ⊗([PauliSpace(Symbol(:s, i)) for i in 1:Ns]...)
sz(i) = Pauli(hs, :σ, 3, i)
sx(i) = Pauli(hs, :σ, 1, i)
sy(i) = Pauli(hs, :σ, 2, i)
sm(i) = (sx(i) - 1im * sy(i)) / 2
@variables J hx γ
Hs = -J * sz(1) * sz(2) - hx * (sx(1) + sx(2))
eqs_s = meanfield([sz(i) for i in 1:Ns], Hs, [sm(i) for i in 1:Ns]; rates = [γ, γ], order = 1)
complete!(eqs_s)
ps = Dict(J => 1.0, hx => 1.0, γ => 0.2)
pd2 = Dict(J => 1.7, hx => 0.4, γ => 0.31)
nst = length(eqs_s.states)
u0s = zeros(ComplexF64, nst)
u = ComplexF64[0.1cos(3.7i) + 0.05im * sin(1.3i) for i in 1:nst]

@testset "cache round trip, sweep after load, corruption" begin
    path = joinpath(mktempdir(), "qc_kernels.jld2")
    kb = KernelBackend(cache = path)
    prob1 = ODEProblem(eqs_s, u0s, (0.0, 1.0), ps; backend = kb)      # miss: lowers, stores
    @test isfile(path)
    prob2 = ODEProblem(eqs_s, u0s, (0.0, 1.0), ps; backend = kb)      # hit: loads, no lowering
    du1, du2 = similar(u0s), similar(u0s)
    prob1.f(du1, u, prob1.p, 0.0)
    prob2.f(du2, u, prob2.p, 0.0)
    @test du1 == du2                                                   # bit-exact round trip

    update_parameters!(prob2, pd2)                                     # a hit is fully sweepable
    fresh = ODEProblem(eqs_s, u0s, (0.0, 1.0), pd2; backend = KernelBackend())
    prob2.f(du2, u, prob2.p, 0.0)
    fresh.f(du1, u, fresh.p, 0.0)
    @test du1 == du2

    # corrupt file degrades to warn + fresh lowering, never a hard error
    write(path, rand(UInt8, 512))
    prob3 = @test_logs (:warn, r"cache") match_mode = :any ODEProblem(
        eqs_s, u0s, (0.0, 1.0), ps; backend = kb
    )
    prob3.f(du2, u, prob3.p, 0.0)
    fresh_ps = ODEProblem(eqs_s, u0s, (0.0, 1.0), ps; backend = KernelBackend())
    fresh_ps.f(du1, u, fresh_ps.p, 0.0)
    @test du1 == du2
end

@testset "tampered entry is a miss" begin
    # write an entry under the CORRECT digest but with altered stored text: the
    # byte-exact verify must treat it as a miss and lower fresh
    path = joinpath(mktempdir(), "qc_tampered.jld2")
    kb = KernelBackend(cache = path)
    probe = ODEProblem(eqs_s, u0s, (0.0, 1.0), ps; backend = kb)       # stores the real entry
    digest = only(jldopen(keys, path, "r"))
    rm(path)
    jldopen(path, "w") do f
        f["$digest/text"] = "TAMPERED"
    end
    prob = ODEProblem(eqs_s, u0s, (0.0, 1.0), ps; backend = kb)        # miss: lowers fresh
    du1, du2 = similar(u0s), similar(u0s)
    prob.f(du1, u, prob.p, 0.0)
    probe.f(du2, u, probe.p, 0.0)
    @test du1 == du2
end
