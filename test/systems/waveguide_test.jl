using QuantumCumulants
using ModelingToolkitBase
using OrdinaryDiffEqTsit5: ODEProblem, solve, Tsit5
using Test

@testset "waveguide" begin
    M_p = 2; M_np = 2; M = M_p + M_np
    h_spin(i) = SpinSpace(Symbol("spin_$(i)"))
    h = tensor([h_spin(i) for i in 1:M]...)
    Sx(i) = Spin(h, Symbol("S$(i)"), 1, i)
    Sy(i) = Spin(h, Symbol("S$(i)"), 2, i)
    Sz(i) = Spin(h, Symbol("S$(i)"), 3, i)
    Sm(i) = (Sx(i) - 1im * Sy(i))
    Sp(i) = (Sx(i) + 1im * Sy(i))
    _Ω_cache = Dict()
    _Γ_cache = Dict()
    Ωp(i, j) = (
        k = i > j ? (j, i) : (i, j);
        get!(_Ω_cache, k, (ModelingToolkitBase.@variables $(Symbol("Ω_$(k[1])_$(k[2])")))[1])
    )
    Γp(i, j) = (
        k = i > j ? (j, i) : (i, j);
        get!(_Γ_cache, k, (ModelingToolkitBase.@variables $(Symbol("Γ_$(k[1])_$(k[2])")))[1])
    )
    H = sum((i ≠ j) * Ωp(i, j) * Sp(i) * Sm(j) for i in 1:M for j in 1:M)
    J = [Sm(c1) for c1 in 1:M]
    rates = [Γp(c1, c2) for c1 in 1:M, c2 in 1:M]
    S(i) = [Sx(i), Sy(i), Sz(i)]
    SiSi(i) = [Sx(i)Sx(i), Sx(i)Sy(i), Sx(i)Sz(i), Sy(i)Sy(i), Sy(i)Sz(i), Sz(i)Sz(i)]
    ops = []
    for i in 1:M
        push!(ops, S(i)...)
    end
    for i in 1:M
        push!(ops, SiSi(i)...)
    end
    for i in 1:M, j in i:M
        if i ≠ j
            for α in 1:3, β in 1:3
                push!(ops, S(i)[α] * S(j)[β])
            end
        end
    end
    @test length(ops) == 90
    eqs = meanfield(ops, H, J; rates = rates, order = 2)
    @test length(eqs.equations) == 90
    sys_c = mtkcompile(System(eqs; name = :wg_sys))
    @test length(unknowns(sys_c)) == 90
    u0 = zeros(ComplexF64, length(eqs.equations))
    for k in 1:M
        for (idx, s) in enumerate(eqs.states)
            if string(s) == string(average(Sz(k)))
                u0[idx] = 0.5 + 0im
            elseif string(s) == string(average(Sz(k) * Sz(k)))
                u0[idx] = 0.25 + 0im
            end
        end
    end
    init = initial_values(eqs, u0)
    pmap_dict = Dict{Any, Any}()
    for i in 1:M, j in (i + 1):M
        pmap_dict[Ωp(i, j)] = 0.1 * (i + j)
    end
    for i in 1:M, j in 1:M
        pmap_dict[Γp(i, j)] = (i == j ? 1.0 : 0.2)
    end
    prob = ODEProblem(sys_c, merge(init, pmap_dict), (0.0, 2.0))
    sol = solve(prob, Tsit5(); reltol = 1.0e-9, abstol = 1.0e-9)
    @test isapprox(
        real(get_solution(sol, Sz(1), eqs).(sol.t[end])),
        -0.029879167913446643; rtol = 1.0e-4,
    )
    @test isapprox(
        real(get_solution(sol, Sz(M), eqs).(sol.t[end])),
        -0.030352188453010635; rtol = 1.0e-4,
    )
end
