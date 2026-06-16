using QuantumCumulants
using Symbolics: Symbolics, @variables
using SymbolicUtils
using ModelingToolkitBase: ModelingToolkitBase, @named, mtkcompile, ODEProblem, unknowns
using OrdinaryDiffEqTsit5: Tsit5, solve, ReturnCode
using Test

# Filter-cavity system combining two indexed Hilbert subspaces; exercises the
# full-system scale + evaluate path for closure and ODE integrability.

@testset "indexed_filter_cavity: meanfield + complete closes the system" begin
    order = 2
    @variables κ::Real g::Real gf::Real κf::Real R::Real Γ::Real
    @variables Δ::Real ν::Real N::Real M::Real
    δ(i) = IndexedVariable(:δ, i)

    hc = FockSpace(:cavity)
    hf = FockSpace(:filter)
    ha = NLevelSpace(:atom, 2)
    h = hc ⊗ hf ⊗ ha

    i = Index(h, :i, M, hf)
    j = Index(h, :j, N, ha)

    @qnumbers a::Destroy(h, 1)
    b(k) = IndexedOperator(Destroy(h, :b, 2), k)
    σ(α, β, k) = IndexedOperator(Transition(h, :σ, α, β, 3), k)

    H = Δ * Σ(σ(2, 2, j), j) +
        Σ(δ(i) * b(i)' * b(i), i) +
        gf * Σ(a' * b(i) + a * b(i)', i) +
        g * Σ(a' * σ(1, 2, j) + a * σ(2, 1, j), j)

    J = [a, b(i), σ(1, 2, j), σ(2, 1, j), σ(2, 2, j)]
    rates = [κ, κf, Γ, R, ν]

    eqs = meanfield([a' * a], H, J; rates = rates, order = order)
    eqs_c = complete(eqs)
    @test length(eqs_c.equations) >= 1
    @test isempty(find_missing(eqs_c))
end

@testset "indexed_filter_cavity: evaluate over the filter index" begin
    order = 2
    @variables κ::Real g::Real gf::Real κf::Real R::Real Γ::Real
    @variables Δ::Real ν::Real N::Real M::Real
    δ(i) = IndexedVariable(:δ, i)

    hc = FockSpace(:cavity)
    hf = FockSpace(:filter)
    ha = NLevelSpace(:atom, 2)
    h = hc ⊗ hf ⊗ ha

    i = Index(h, :i, M, hf)
    j = Index(h, :j, N, ha)

    @qnumbers a::Destroy(h, 1)
    b(k) = IndexedOperator(Destroy(h, :b, 2), k)
    σ(α, β, k) = IndexedOperator(Transition(h, :σ, α, β, 3), k)

    H = Δ * Σ(σ(2, 2, j), j) +
        Σ(δ(i) * b(i)' * b(i), i) +
        gf * Σ(a' * b(i) + a * b(i)', i) +
        g * Σ(a' * σ(1, 2, j) + a * σ(2, 1, j), j)

    J = [a, b(i), σ(1, 2, j), σ(2, 1, j), σ(2, 2, j)]
    rates = [κ, κf, Γ, R, ν]

    eqs_c = complete(meanfield([a' * a], H, J; rates = rates, order = order))
    n_before = length(eqs_c.equations)
    evaled = evaluate(eqs_c; limits = (M => 3))
    @test length(evaled.equations) >= n_before
end

@testset "indexed_filter_cavity: per-Hilbert-space evaluate/scale commute" begin
    # Hybrid system with two indexed subspaces, distinct couplings on the
    # filter side (`δ(i)`) and permutation-symmetric atoms. The two
    # materialisation orders must give equation-count agreement:
    #  - evaluate(h=[2]) unrolls the filter, then scale(h=[3]) collapses atoms
    #  - scale(h=[3]) collapses atoms, then evaluate(M=>3) unrolls the filter
    order = 2
    @variables κ::Real g::Real gf::Real κf::Real R::Real Γ::Real
    @variables Δ::Real ν::Real N::Real M::Real
    δ(i) = IndexedVariable(:δ, i)

    hc = FockSpace(:cavity)
    hf = FockSpace(:filter)
    ha = NLevelSpace(:atom, 2)
    h = hc ⊗ hf ⊗ ha

    i = Index(h, :i, M, hf)
    j = Index(h, :j, N, ha)

    @qnumbers a::Destroy(h, 1)
    b(k) = IndexedOperator(Destroy(h, :b, 2), k)
    σ(α, β, k) = IndexedOperator(Transition(h, :σ, α, β, 3), k)

    H = Δ * Σ(σ(2, 2, j), j) +
        Σ(δ(i) * b(i)' * b(i), i) +
        gf * Σ(a' * b(i) + a * b(i)', i) +
        g * Σ(a' * σ(1, 2, j) + a * σ(2, 1, j), j)

    J = [a, b(i), σ(1, 2, j), σ(2, 1, j), σ(2, 2, j)]
    rates = [κ, κf, Γ, R, ν]

    eqs_c = complete(meanfield([a' * a], H, J; rates = rates, order = order))

    M_ = 3
    eqs_eval_ = evaluate(eqs_c; limits = (M => M_), h = [2])
    eqs_sc_ = scale(eqs_eval_; h = [3])

    eqs_sc = scale(eqs_c; h = [3])
    eqs_eval = evaluate(eqs_sc; limits = (M => M_))

    @test length(eqs_sc_.equations) == length(eqs_eval.equations)
end

@testset "indexed_filter_cavity: ODE solve, photon-number physicality (cavity only)" begin
    # Closed atom+cavity system (no filter modes); ⟨a'a⟩ is real and
    # non-negative along the trajectory.
    @variables κ::Real g::Real R::Real Γ::Real Δ::Real ν::Real N::Real

    hc = FockSpace(:cavity)
    ha = NLevelSpace(:atom, 2)
    h = hc ⊗ ha

    j = Index(h, :j, N, ha)

    @qnumbers a::Destroy(h, 1)
    σ(α, β, k) = IndexedOperator(Transition(h, :σ, α, β, 2), k)

    H = Δ * Σ(σ(2, 2, j), j) +
        g * Σ(a' * σ(1, 2, j) + a * σ(2, 1, j), j)
    J = [a, σ(1, 2, j), σ(2, 1, j), σ(2, 2, j)]
    rates = [κ, Γ, R, ν]

    eqs_c = complete(meanfield([a' * a], H, J; rates = rates, order = 2))
    eqs_ev = evaluate(eqs_c; limits = (N => 2))
    @test isempty(find_missing(eqs_ev))

    SQA = QuantumCumulants.SecondQuantizedAlgebra
    @named sys_ev = System(eqs_ev)
    sys_c = mtkcompile(sys_ev)
    pars = ModelingToolkitBase.parameters(sys_c)
    pmap = Dict{Any, Any}()
    for p in pars
        ps = string(p)
        if occursin("κ", ps)
            pmap[p] = 0.1
        elseif occursin("g", ps)
            pmap[p] = 1.0
        elseif occursin("R", ps)
            pmap[p] = 0.3
        elseif occursin("Γ", ps)
            pmap[p] = 0.05
        elseif occursin("ν", ps)
            pmap[p] = 0.0
        elseif occursin("Δ", ps)
            pmap[p] = 0.0
        else
            pmap[p] = 0.0
        end
    end
    u0 = Dict(unknowns(sys_c) .=> zeros(ComplexF64, length(unknowns(sys_c))))
    prob = ODEProblem(sys_c, merge(u0, pmap), (0.0, 20.0))
    sol = solve(prob, Tsit5(); abstol = 1.0e-9, reltol = 1.0e-9)
    @test sol.retcode == ReturnCode.Success

    a_op = SQA.undo_average(eqs_ev.states[1])
    assert_real(sol, a_op, eqs_ev; atol = 1.0e-7)
    assert_nonneg(sol, a_op, eqs_ev; atol = 1.0e-7)
end
