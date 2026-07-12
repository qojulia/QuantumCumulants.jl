using Test
using CumulantHomotopy
using QuantumCumulants
using Symbolics
using SymbolicUtils
using SecondQuantizedAlgebra
const SQA = SecondQuantizedAlgebra

# Numerically evaluate the original (complex) meanfield RHS at a reconstructed moment
# solution, to check the realified/HC solution against the untouched hierarchy.
function complex_residual(eqs, moment_values, params)
    psub = Dict(SymbolicUtils.unwrap(k) => v for (k, v) in params)
    mvals = Dict(SymbolicUtils.unwrap(k) => v for (k, v) in moment_values)
    isavg(v) = SymbolicUtils.iscall(v) && SymbolicUtils.operation(v) isa SQA.AvgFunc
    conjleaf(v) = SymbolicUtils.unwrap(SQA.average(SQA.undo_average(v)'))
    function ev(ex)
        ex = SymbolicUtils.unwrap(ex)
        ex isa Number && return ComplexF64(ex)
        ex === Symbolics.IM && return im
        if isavg(ex)
            haskey(mvals, ex) && return mvals[ex]
            return conj(mvals[conjleaf(ex)])          # conjugate moment
        end
        SymbolicUtils.isconst(ex) && return ComplexF64(ex.val)
        if !SymbolicUtils.iscall(ex)
            haskey(psub, ex) && return ComplexF64(psub[ex])
            error("unknown leaf $ex")
        end
        op = SymbolicUtils.operation(ex)
        return op([ev(a) for a in SymbolicUtils.arguments(ex)]...)
    end
    return [abs(ev(eq.rhs)) for eq in eqs.equations]
end

@testset "CumulantHomotopy" begin
    @test isdefined(CumulantHomotopy, :stationary_state)
    @test isdefined(CumulantHomotopy, :stationary_sequence)

    @testset "driven Kerr resonator (order 2)" begin
        @variables Δ::Real U::Real κ::Real F::Real
        h = FockSpace(:cavity)
        @qnumbers a::Destroy(h)
        H = Δ * a' * a + U * a' * a' * a * a + F * (a + a')
        eqs = complete(meanfield(a, H, [a]; rates = [κ], order = 2))

        # Realification produces a square real polynomial system (SPEC §3.2).
        sys = realify(eqs)
        @test sys isa StationaryPolynomialSystem
        @test length(sys.equations) == length(sys.variables) == 5
        @test Set(string.(sys.parameters)) == Set(["Δ", "U", "κ", "F"])

        params = Dict(Δ => 1.0, U => 0.2, κ => 1.0, F => 2.0)
        res = stationary_states(sys, params; show_progress = false)

        # A bistable working point: three real (physical) roots.
        @test length(res) == 3

        for sol in res
            # Every physical root solves the *original* complex hierarchy.
            @test maximum(complex_residual(eqs, sol, params)) < 1e-8
            # Photon number ⟨a'a⟩ is real.
            n = sol[SQA.average(a' * a)]
            @test abs(imag(n)) < 1e-8
        end

        # The largest photon-number branch is a genuine (positive) occupation.
        ns = [real(sol[SQA.average(a' * a)]) for sol in res]
        @test maximum(ns) > 0

        # Enumerating all complex roots returns at least as many solutions.
        res_all = stationary_states(sys, params; only_physical = false, show_progress = false)
        @test length(res_all) >= length(res)
    end
end
