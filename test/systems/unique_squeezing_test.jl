using QuantumCumulants
using ModelingToolkitBase
using ModelingToolkitBase: t_nounits as t
using Symbolics: @variables
using OrdinaryDiffEqTsit5: ODEProblem, solve, Tsit5
using Test

@testset "unique_squeezing" begin
    hf = FockSpace(:harmonic); ha = NLevelSpace(:spin, 2); h = hf ⊗ ha
    @variables ω Ω ωd η κ g γ N ξ
    @qnumbers a::Destroy(h)
    σ(x, y, idx) = IndexedOperator(Transition(h, :σ, x, y), idx)
    b = a * cosh(ξ) + a' * sinh(ξ)
    i = Index(h, :i, N, ha); j = Index(h, :j, N, ha)
    Hf = ω * a' * a + η * (b' * exp(-1im * ωd * t) + b * exp(1im * ωd * t))
    Ha = Ω * Σ(σ(2, 2, i) - σ(1, 1, i), i) / 2
    Hi = g * Σ((σ(1, 2, i) + σ(2, 1, i)) * (a + a'), i) / 2
    H = Hf + Ha + Hi
    J = [b, σ(1, 2, i)]
    rates = [κ, γ]
    eqs = meanfield([a, a' * a, σ(2, 2, j)], H, J; rates = rates, order = 2)
    eqs_c = complete(eqs)
    @test length(eqs_c.equations) == 13
    eqs_sc = scale(eqs_c)
    @test length(eqs_sc.equations) == 12
    sys_c = mtkcompile(System(eqs_sc; name = :us_sys))
    @test length(unknowns(sys_c)) == 12
    ω_ = 1.0; Ω_ = 2.0e3; η_ = 4.0; κ_ = 1.0; γ_ = 1.0
    N_ = 1
    gc_ = sqrt(Ω_ * ω_ / N_)
    g_ = 0.9 * gc_
    ωd_ = sqrt(1 - g_^2 / gc_^2) * ω_
    ξ_ = (1 / 4) * log(1 - N_ * g_^2 / (ω_ * Ω_))
    u0 = zeros(ComplexF64, length(eqs_sc.equations))
    init = initial_values(eqs_sc, u0)
    ps = Dict(
        ω => ω_, Ω => Ω_, ωd => ωd_, g => g_, η => η_,
        κ => κ_, γ => γ_, N => Float64(N_), ξ => ξ_,
    )
    tend = Float64(2π / ωd_)
    prob = ODEProblem(sys_c, merge(init, ps), (0.0, tend))
    sol = solve(prob, Tsit5(); reltol = 1.0e-8, abstol = 1.0e-8)
    @test isapprox(
        real(get_solution(sol, a' * a, eqs_sc).(sol.t[end])),
        30.149750050640534; rtol = 1.0e-3,
    )
end

@testset "unique_squeezing (free-j shape, N=100 plateau)" begin
    hf = FockSpace(:harmonic); ha = NLevelSpace(:spin, 2); h = hf ⊗ ha
    @variables ω Ω ωd η κ g γ N ξ
    @qnumbers a::Destroy(h)
    σ(x, y, idx) = IndexedOperator(Transition(h, :σ, x, y), idx)
    b = a * cosh(ξ) + a' * sinh(ξ)
    i = Index(h, :i, N, ha); j = Index(h, :j, N, ha)
    Hf = ω * a' * a + η * (b' * exp(-1im * ωd * t) + b * exp(1im * ωd * t))
    Ha = Ω * Σ(σ(2, 2, i) - σ(1, 1, i), i) / 2
    Hi = g * Σ((σ(1, 2, i) + σ(2, 1, i)) * (a + a'), i) / 2
    H = Hf + Ha + Hi
    J = [b, σ(1, 2, i)]
    rates = [κ, γ]
    eqs = meanfield([a, a' * a, σ(2, 2, j)], H, J; rates = rates, order = 2)
    eqs_c = complete(eqs)
    @test length(eqs_c.equations) == 13
    eqs_sc = scale(eqs_c)
    @test length(eqs_sc.equations) == 12
    sys_c = mtkcompile(System(eqs_sc; name = :us_freej_sys))
    @test length(unknowns(sys_c)) == 12
    ω_ = 1.0; Ω_ = 2.0e3; η_ = 4.0; κ_ = 1.0; γ_ = 1.0
    N_ = 100
    gc_ = sqrt(Ω_ * ω_ / N_)
    g_ = 0.9 * gc_
    ωd_ = sqrt(1 - g_^2 / gc_^2) * ω_
    ξ_ = (1 / 4) * log(1 - N_ * g_^2 / (ω_ * Ω_))
    u0 = zeros(ComplexF64, length(eqs_sc.equations))
    u0d = Dict{Any, Any}(unknowns(sys_c) .=> u0)
    ps = Dict{Any, Any}(
        ω => ω_, Ω => Ω_, ωd => ωd_, g => g_, η => η_,
        κ => κ_, γ => γ_, N => Float64(N_), ξ => ξ_,
    )
    tend = Float64(4π / ωd_)
    prob = ODEProblem(sys_c, merge(u0d, ps), (0.0, tend))
    sol = solve(
        prob, Tsit5(); saveat = π / 30ωd_,
        reltol = 1.0e-10, abstol = 1.0e-10,
    )
    t_ = sol.t
    adag_a = get_solution(sol, a' * a, eqs_sc).(t_)
    aa = get_solution(sol, a * a, eqs_sc).(t_)
    adag_adag = get_solution(sol, a' * a', eqs_sc).(t_)
    a_ = get_solution(sol, a, eqs_sc).(t_)
    adag = get_solution(sol, a', eqs_sc).(t_)
    sqx = real.(adag_adag + aa + 2 * adag_a .+ 1 - (adag + a_) .^ 2)
    sqy = real.(adag_adag + aa - 2 * adag_a .- 1 - (adag - a_) .^ 2)
    # N=100 squeezing plateau, Gietka et al., PRL 131, 223604: X ≈ 2.29, P ≈ 0.44.
    @test isapprox(sqx[end], 2.292622635966603; rtol = 1.0e-3)
    @test isapprox(-sqy[end], 0.43703195626758884; rtol = 1.0e-3)
end
