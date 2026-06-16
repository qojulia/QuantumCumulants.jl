using QuantumCumulants
using Symbolics: Symbolics, @variables, expand
using SymbolicUtils
using ModelingToolkitBase: ModelingToolkitBase, brownians, equations
using Test

function _iz(x)
    x isa Number && return iszero(x)
    x isa SymbolicUtils.BasicSymbolic || return iszero(x)
    return SymbolicUtils._iszero(SymbolicUtils.simplify(x; expand = true))
end

@testset "forward noise meanfield: shape" begin
    hc = FockSpace(:cavity)
    @qnumbers a::Destroy(hc)
    @variables ω κ η
    H = ω * a' * a
    eqs = meanfield([a, a' * a], H, [a]; rates = [κ], efficiencies = [η], order = 2)
    @test eqs isa NoiseMeanfieldEquations
    @test length(eqs.equations) == 2
    @test length(eqs.noise_equations) == 2
    @test eqs.direction isa Forward
end

@testset "backward noise meanfield: shape" begin
    hc = FockSpace(:cavity)
    @qnumbers a::Destroy(hc)
    @variables ω κ η
    H = ω * a' * a
    eqs = meanfield(
        [a, a' * a], H, [a];
        rates = [κ], efficiencies = [η],
        direction = Backward(), order = 2
    )
    @test eqs isa NoiseMeanfieldEquations
    @test eqs.direction isa Backward
end

@testset "measurement backaction: damped cavity matches analytic noise term" begin
    @variables ω::Real κ::Real η::Real
    hc = FockSpace(:cavity)
    a = Destroy(hc, :a)
    eqs = meanfield(a, ω * a' * a, [a]; rates = [κ], efficiencies = [η], order = 2)

    test_noise = sqrt(κ * η) *
        (average(a' * a + a * a) - average(a') * average(a) - average(a)^2)
    @test isequal(eqs[1].lhs, average(a))
    @test _iz(test_noise - eqs.noise_equations[1].rhs)
end

@testset "measurement backaction: zero efficiencies zero out noise" begin
    @variables ω::Real κ::Real γ::Real ωa::Real
    ha = NLevelSpace(:atom, 2)
    hc = FockSpace(:cavity)
    h = ha ⊗ hc
    a = Destroy(h, :a)
    σ = Transition(h, :σ, 2, 1)

    eqs = meanfield(
        [a, σ], ω * a' * a + ωa * σ' * σ, [a, σ];
        rates = [κ, γ], efficiencies = [0, 0]
    )
    @test _iz(eqs.noise_equations[1].rhs)
    @test _iz(eqs.noise_equations[2].rhs)
end

@testset "measurement backaction: full drift + noise match analytic" begin
    @variables Δ::Real κ::Real η::Real
    hc = FockSpace(:cavity)
    a = Destroy(hc, :a)

    me_a = meanfield(a, Δ * a' * a, [a]; rates = [κ], efficiencies = [η], order = 2)
    # Build comparison targets via `average` of an operator-level expression so
    # the imaginary unit is composed with `Symbolics.IM` and stays comparable
    # under `simplify`/`expand`.
    test_eq_a = average(-1im * Δ * a - (1 // 2) * κ * a)
    test_noise_eq_a = sqrt(η * κ) *
        (average(a' * a) + average(a * a) - average(a)^2 - average(a) * average(a'))
    @test _iz(me_a.equations[1].rhs - test_eq_a)
    @test _iz(me_a.noise_equations[1].rhs - test_noise_eq_a)

    # At order 1 the noise contribution to ⟨a'⟩ under jump a vanishes.
    me_a_o1 = meanfield(a', Δ * a' * a, [a]; rates = [κ], efficiencies = [η], order = 1)
    @test _iz(me_a_o1.noise_equations[1].rhs)
end

@testset "measurement backaction: drift matches deterministic when det/stoch are split" begin
    # The drift equations are identical whether or not noise is requested.
    @variables Δ::Real κ::Real η::Real
    hc = FockSpace(:cavity)
    a = Destroy(hc, :a)
    H = Δ * a' * a

    det = meanfield(a, H, [a]; rates = [κ], order = 2)
    stoch = meanfield(a, H, [a]; rates = [κ], efficiencies = [η], order = 2)
    @test _iz(stoch.equations[1].rhs - det.equations[1].rhs)
end

@testset "MeanfieldEquations(::NoiseMeanfieldEquations) + parameter_map(sys, …)" begin
    # `MeanfieldEquations` strips the noise drift (and parameters living only
    # inside it); `parameter_map(sys, …)` filters a shared user dict down to the
    # parameters each system actually uses.
    @variables Δ::Real κ::Real η::Real
    hc = FockSpace(:cavity)
    a = Destroy(hc, :a)
    H = Δ * a' * a
    # Track the full second-order set so the noise drift's `⟨a'·a⟩` / `⟨a·a⟩`
    # references reduce to states and the SDE compile sees a linear Brownian.
    eqs_noise = meanfield(
        [a, a' * a, a * a], H, [a];
        rates = [κ], efficiencies = [η], order = 2
    )
    eqs_det = MeanfieldEquations(eqs_noise)
    @test eqs_det isa MeanfieldEquations
    @test eqs_det.direction === eqs_noise.direction
    @test length(eqs_det.equations) == length(eqs_noise.equations)

    sys_det = ModelingToolkitBase.mtkcompile(System(eqs_det; name = :det))
    sys_st = ModelingToolkitBase.mtkcompile(System(eqs_noise; name = :stoch))

    # η lives only in the noise drift, so it must be stripped from the
    # deterministic map and kept for the stochastic one.
    full_pmap = Dict(Δ => 1.0, κ => 0.5, η => 0.3)
    det_pmap = parameter_map(sys_det, full_pmap)
    st_pmap = parameter_map(sys_st, full_pmap)
    @test !haskey(det_pmap, η)
    @test haskey(det_pmap, Δ) && haskey(det_pmap, κ)
    @test haskey(st_pmap, η)
end

@testset "System on noise eqs: substitutes conjugate of state on RHS" begin
    # Driven damped cavity: the drift of ⟨a'a⟩ references ⟨a'⟩, which is the
    # conjugate of state ⟨a⟩ rather than a separate state, so the SDE codegen
    # must rewrite it as `conj(u_for_a(t))`.
    @variables Δ::Real κ::Real η::Real η_d::Real
    hc = FockSpace(:cavity)
    a = Destroy(hc, :a)
    H = Δ * a' * a + η_d * (a + a')
    eqs = meanfield(
        [a, a' * a], H, [a];
        rates = [κ], efficiencies = [η], order = 2
    )

    @test length(eqs.states) == 2
    @test eqs isa NoiseMeanfieldEquations

    sys = System(complete(eqs); name = :driven_cav_noise)
    @test sys isa ModelingToolkitBase.System
end

@testset "System on noise eqs: independent monitored channels get independent brownians" begin
    hc = FockSpace(:cavity)
    ha = NLevelSpace(:atom, 2)
    h = hc ⊗ ha
    a = Destroy(h, :a, 1)
    σ(i, j) = Transition(h, :σ, i, j, 2)

    @variables Δ::Real κ::Real γ::Real η1::Real η2::Real
    H = Δ * a' * a
    eqs = meanfield(
        [a, σ(2, 2)], H, [a, σ(1, 2)];
        rates = [κ, γ], efficiencies = [η1, η2], order = 2
    )

    sys = System(complete(eqs); name = :two_channel_noise)
    ws = brownians(sys)

    @test length(ws) == 2
    eq_text = join(string.(equations(sys)), "\n")
    @test occursin(string(ws[1]), eq_text)
    @test occursin(string(ws[2]), eq_text)
end
