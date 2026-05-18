using QuantumCumulants
using Symbolics: Symbolics, @variables, expand
using SymbolicUtils
using ModelingToolkitBase: ModelingToolkitBase
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
    @test eqs isa NoiseMeanFieldEquations
    @test length(eqs.equations) == 2
    @test length(eqs.noise_equations) == 2
    @test eqs.direction isa Forward
end

@testset "backward noise meanfield: shape" begin
    hc = FockSpace(:cavity)
    @qnumbers a::Destroy(hc)
    @variables ω κ η
    H = ω * a' * a
    eqs = meanfield([a, a' * a], H, [a];
                    rates = [κ], efficiencies = [η],
                    direction = Backward(), order = 2)
    @test eqs isa NoiseMeanFieldEquations
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

    eqs = meanfield([a, σ], ω * a' * a + ωa * σ' * σ, [a, σ];
                    rates = [κ, γ], efficiencies = [0, 0])
    @test _iz(eqs.noise_equations[1].rhs)
    @test _iz(eqs.noise_equations[2].rhs)
end

@testset "measurement backaction: full drift + noise match analytic" begin
    @variables Δ::Real κ::Real η::Real
    hc = FockSpace(:cavity)
    a = Destroy(hc, :a)

    me_a = meanfield(a, Δ * a' * a, [a]; rates = [κ], efficiencies = [η], order = 2)
    # Build comparison targets by averaging an operator-level expression rather
    # than multiplying averages by complex scalars. The latter promotes through
    # `Complex{Num}` and Symbolics materialises the literal `complex(re, im)`
    # symbolic term, which is opaque to `simplify`/`expand`. Going via
    # `average(op-with-complex-coeff)` routes through SQA's `average(::QAdd)`,
    # which composes the imaginary unit using `Symbolics.IM`.
    test_eq_a = average(-1im * Δ * a - (1//2) * κ * a)
    test_noise_eq_a = sqrt(η * κ) *
        (average(a' * a) + average(a * a) - average(a)^2 - average(a) * average(a'))
    @test _iz(me_a.equations[1].rhs - test_eq_a)
    @test _iz(me_a.noise_equations[1].rhs - test_noise_eq_a)

    # Order-1: noise contribution vanishes for ⟨a'⟩ under jump a (cumulant-trunc'd)
    me_a_o1 = meanfield(a', Δ * a' * a, [a]; rates = [κ], efficiencies = [η], order = 1)
    @test _iz(me_a_o1.noise_equations[1].rhs)
end

@testset "measurement backaction: drift matches deterministic when det/stoch are split" begin
    # The drift equations should be identical whether we ask for noise or not.
    @variables Δ::Real κ::Real η::Real
    hc = FockSpace(:cavity)
    a = Destroy(hc, :a)
    H = Δ * a' * a

    det = meanfield(a, H, [a]; rates = [κ], order = 2)
    stoch = meanfield(a, H, [a]; rates = [κ], efficiencies = [η], order = 2)
    @test _iz(stoch.equations[1].rhs - det.equations[1].rhs)
end

@testset "to_system on noise eqs: substitutes conjugate of state on RHS" begin
    # Driven damped cavity. The drift of ⟨a'a⟩ references ⟨a'⟩, but ⟨a'⟩ is
    # the conjugate of state ⟨a⟩ and is therefore not itself added as a
    # separate state. The SDE codegen path must rewrite ⟨a'⟩ as
    # `conj(u_for_a(t))` before substituting averages into MTK variables;
    # without that rewrite it would emit a literal `avg(a')` call and the
    # solver would later fail with "AvgFunc not callable".
    @variables Δ::Real κ::Real η::Real η_d::Real
    hc = FockSpace(:cavity)
    a = Destroy(hc, :a)
    H = Δ * a' * a + η_d * (a + a')
    eqs = meanfield([a, a' * a], H, [a];
                    rates = [κ], efficiencies = [η], order = 2)

    @test length(eqs.states) == 2
    @test eqs isa NoiseMeanFieldEquations

    sys = to_system(eqs; name = :driven_cav_noise)
    @test sys isa ModelingToolkitBase.System
end
