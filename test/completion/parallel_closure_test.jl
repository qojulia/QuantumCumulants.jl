using QuantumCumulants
using Symbolics: @variables
using Test

const QC = QuantumCumulants
const SQA = QuantumCumulants.SecondQuantizedAlgebra

function _assert_same_graph(reference::QC.MomentGraph, candidate::QC.MomentGraph)
    rkeys = collect(keys(reference.nodes))
    ckeys = collect(keys(candidate.nodes))
    @test length(rkeys) == length(ckeys)
    @test all(i -> isequal(rkeys[i], ckeys[i]), eachindex(rkeys))
    for i in eachindex(rkeys)
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

function _ising_graph(; noise::Bool = false)
    N = 3
    h = ⊗([PauliSpace(Symbol(:spin, i)) for i in 1:N]...)
    σx(i) = Pauli(h, :σ, 1, i)
    σy(i) = Pauli(h, :σ, 2, i)
    σz(i) = Pauli(h, :σ, 3, i)
    σm(i) = (σx(i) - 1im * σy(i)) / 2
    @variables J hx γ η
    H = -J * sum(σz(i) * σz(i + 1) for i in 1:(N - 1)) - hx * sum(σx(i) for i in 1:N)
    if noise
        return meanfield(
            [σz(i) for i in 1:N], H, [σm(i) for i in 1:N];
            rates = fill(γ, N), efficiencies = fill(η, N), order = 2,
        ).graph
    end
    return meanfield(
        [σz(i) for i in 1:N], H, [σm(i) for i in 1:N];
        rates = fill(γ, N), order = 2,
    ).graph
end

function _indexed_collective_graph()
    h = NLevelSpace(:atom, 2)
    @variables N
    i = Index(h, :i, N, h)
    j = Index(h, :j, N, h)
    σ(x, y, k) = IndexedOperator(Transition(h, :σ, x, y), k)
    return meanfield(
        [σ(1, 2, i), σ(2, 2, i)],
        0 * Σ(σ(2, 2, i), i),
        [σ(1, 2, i)];
        rates = [DoubleIndexedVariable(:Γ, i, j)], order = 2,
    ).graph
end

_serial(g; kwargs...) = QC._closure(g, QC._derive_frontier_serial; kwargs...)
_threaded(g; kwargs...) = QC._closure(g, QC._derive_frontier_polyester; kwargs...)

@testset "parallel closure preserves deterministic hierarchy semantics" begin
    cases = (
        (; get_adjoints = true),
        (; get_adjoints = false),
        (; get_adjoints = false, foldable = _ -> false),
        (;
            get_adjoints = false,
            filter = avg -> length(SQA.operators(SQA.undo_average(avg))) <= 1,
        ),
    )

    for kwargs in cases
        reference = _serial(_ising_graph(); kwargs...)
        for _ in 1:5
            candidate = _threaded(_ising_graph(); kwargs...)
            _assert_same_graph(reference, candidate)
        end
    end
end

@testset "parallel closure preserves noise leaves" begin
    reference = _serial(_ising_graph(; noise = true); get_adjoints = false)
    @test any(nd -> nd.noise !== nothing, values(reference.nodes))
    for _ in 1:5
        candidate = _threaded(_ising_graph(; noise = true); get_adjoints = false)
        _assert_same_graph(reference, candidate)
    end
end

@testset "parallel closure preserves indexed collective-rate derivation" begin
    reference = _serial(_indexed_collective_graph(); get_adjoints = false)
    for _ in 1:10
        candidate = _threaded(_indexed_collective_graph(); get_adjoints = false)
        _assert_same_graph(reference, candidate)
    end
end

@testset "public closure selects the same hierarchy" begin
    reference = _serial(_ising_graph(); get_adjoints = false)
    candidate = QC.closure(_ising_graph(); get_adjoints = false)
    _assert_same_graph(reference, candidate)
end
