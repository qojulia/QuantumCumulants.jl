using QuantumCumulants
using Symbolics: @variables
using Test

const QC = QuantumCumulants

function assert_same_closure_graph(reference, candidate)
    rkeys = collect(keys(reference.nodes))
    ckeys = collect(keys(candidate.nodes))
    @test length(rkeys) == length(ckeys)
    for i in eachindex(rkeys)
        @test isequal(rkeys[i], ckeys[i])
        r = reference.nodes[rkeys[i]]
        c = candidate.nodes[ckeys[i]]
        @test isequal(r.drift, c.drift)
        @test isequal(r.op_drift, c.op_drift)
        @test isequal(r.noise, c.noise)
        @test isequal(r.op_noise, c.op_noise)
        @test r.order == c.order
        @test r.aon == c.aon
    end
    return nothing
end

function jc_closure_graph()
    hc = FockSpace(:cavity)
    ha = NLevelSpace(:atom, 2)
    h = hc ⊗ ha
    a = Destroy(h, :a, 1)
    σ(i, j) = Transition(h, :σ, i, j, 2)
    @variables Δ::Real g::Real Ω::Real κ::Real γ::Real
    H = Δ * a' * a + g * (a' * σ(1, 2) + a * σ(2, 1)) + Ω * (a + a')
    eqs = meanfield(
        [a' * a, σ(2, 2)], H, [a, σ(1, 2)]; rates = [κ, γ], order = 2,
    )
    return eqs.graph
end

@testset "closure: Polyester frontier matches serial semantics" begin
    g = jc_closure_graph()
    cases = (
        ("full adjoints", (; get_adjoints = true)),
        ("folded adjoints", (; get_adjoints = false)),
        ("non-foldable adjoints", (; get_adjoints = false, foldable = _ -> false)),
        ("filtered", (; get_adjoints = false, filter = phase_invariant)),
    )
    for (name, kwargs) in cases
        @testset "$name" begin
            serial = QC._closure(g, QC._derive_frontier_serial; kwargs...)
            for _ in 1:3
                parallel = QC._closure(g, QC._derive_frontier_parallel; kwargs...)
                assert_same_closure_graph(serial, parallel)
            end
        end
    end
end

@testset "closure: indexed collective-rate frontier is deterministic" begin
    h = NLevelSpace(:atom, 2)
    @variables N::Real
    i = Index(h, :i, N, h)
    j = Index(h, :j, N, h)
    σ(x, y, k) = IndexedOperator(Transition(h, :σ, x, y), k)
    eqs = meanfield(
        [σ(1, 2, i), σ(2, 2, i)],
        0 * Σ(σ(2, 2, i), i),
        [σ(1, 2, i)];
        rates = [DoubleIndexedVariable(:Γ, i, j)],
        order = 2,
    )
    serial = QC._closure(eqs.graph, QC._derive_frontier_serial; get_adjoints = false)
    for _ in 1:5
        parallel = QC._closure(eqs.graph, QC._derive_frontier_parallel; get_adjoints = false)
        assert_same_closure_graph(serial, parallel)
    end
end

@testset "closure: stochastic frontier includes noise leaves exactly" begin
    h = FockSpace(:cavity)
    a = Destroy(h, :a)
    @variables Δ::Real κ::Real η::Real
    eqs = meanfield(
        [a, a' * a], Δ * a' * a, [a];
        rates = [κ], efficiencies = [η], order = 2,
    )
    serial = QC._closure(eqs.graph, QC._derive_frontier_serial; get_adjoints = false)
    parallel = QC._closure(eqs.graph, QC._derive_frontier_parallel; get_adjoints = false)
    assert_same_closure_graph(serial, parallel)
end
