using QuantumCumulants
using Symbolics: @variables
using LinearAlgebra: Diagonal
using Test

const SU = QuantumCumulants.SymbolicUtils
const Sym = QuantumCumulants.Symbolics

# `phase_invariant` is the shared U(1) filter from runtests.jl `init_code`.

# Build a spectrum without solving: the kernel is a pure function of the symbolic
# τ-system, so block-level checks only need a fabricated (u_end, p0).
function _laser_spectrum()
    @variables Δ g γ κ ν
    hf = FockSpace(:cavity); ha = NLevelSpace(:atom, (:g, :e)); h = hf ⊗ ha
    a = Destroy(h, :a); s = Transition(h, :σ, :g, :e)
    H = Δ * a' * a + g * (a' * s + a * s')
    eqs = complete(
        meanfield(a' * a, H, [a, s, s']; rates = [κ, γ, ν], order = 2);
        filter_func = phase_invariant,
    )
    c = CorrelationFunction(a', a, eqs; steady_state = true, filter_func = phase_invariant)
    S = Spectrum(c, (Δ, g, γ, κ, ν))
    p0 = (0.0, 1.5, 0.25, 1.0, 4.0)
    u_end = Dict{Any, ComplexF64}()
    for (i, st) in enumerate(eqs.states)
        u_end[st] = (0.3 + 0.1i) + (0.2 - 0.05i) * im
    end
    return S, u_end, p0
end

@testset "spectrum kernel: CSE equivalence" begin
    S, u_end, p0 = _laser_spectrum()
    resolve = QuantumCumulants._ss_resolver(S.c, u_end)
    kern0 = QuantumCumulants._build_spectrum_kernel(S; cse = false)
    kern1 = QuantumCumulants._build_spectrum_kernel(S; cse = true)
    A0, _, b0 = QuantumCumulants._eval_spectrum_blocks(kern0, S, u_end, p0, resolve)
    A1, _, b1 = QuantumCumulants._eval_spectrum_blocks(kern1, S, u_end, p0, resolve)
    @test isapprox(A0, A1; atol = 1.0e-12)
    @test isapprox(b0, b1; atol = 1.0e-12)
    # the off-diagonal coupling (the IM-bearing terms) must be nonzero: a regression guard
    # against the imaginary unit being captured as a build_function input and zeroed.
    @test any(!iszero, A0 - Diagonal(A0))
end

@testset "spectrum kernel: cache reuse" begin
    S, u_end, p0 = _laser_spectrum()
    @test S.cache[] === nothing
    A, rhs_b, n = QuantumCumulants._spectrum_kernel(S, u_end, p0)
    @test S.cache[] isa QuantumCumulants.SpectrumKernel
    kern = S.cache[]
    # a second call reuses the same compiled kernel object and reproduces the system.
    A2, rhs_b2, n2 = QuantumCumulants._spectrum_kernel(S, u_end, p0)
    @test S.cache[] === kern
    @test n == n2
    @test isapprox(A, A2; atol = 1.0e-12)
    @test isapprox(rhs_b, rhs_b2; atol = 1.0e-12)
end

@testset "spectrum kernel: linearity assertion" begin
    # A coefficient block that still depends on a state placeholder means the τ-system is
    # nonlinear in its states; the build must reject it rather than emit a wrong A.
    x1 = SU.unwrap(Sym.variable(:__spec_x1; T = Number))
    x2 = SU.unwrap(Sym.variable(:__spec_x2; T = Number))
    @variables p::Real
    pu = SU.unwrap(p)
    state_ph = Set(Any[x1, x2])
    good_block = Any[pu + 2.0, pu]               # linear: free of state placeholders
    bad_block = Any[pu * x1, pu]                 # x1 leftover -> nonlinear
    @test QuantumCumulants._assert_linear_blocks(Any[good_block], state_ph) === nothing
    @test_throws ArgumentError QuantumCumulants._assert_linear_blocks(
        Any[bad_block], state_ph,
    )
end
