using QuantumCumulants
import Symbolics
using Symbolics: @variables
using Test

const QC = QuantumCumulants

@testset "ground reduction capability and serialized traversal" begin
    hp = PauliSpace(:spin)
    sx = Pauli(hp, :σ, 1)
    sy = Pauli(hp, :σ, 2)
    sz = Pauli(hp, :σ, 3)
    @variables Ω γ
    peqs = meanfield([sz], Ω * sx, [(sx - 1im * sy) / 2]; rates = [γ], order = 2)
    pk = first(keys(peqs.graph.nodes))
    pnode = QC.derive(pk, peqs.graph.sys, peqs.graph.ctx)
    @test !QC._may_need_ground_reduction(pnode.op_drift)
    raw = Symbolics.Num(
        QC.average_and_truncate(
            pnode.op_drift, peqs.graph.sys.order, peqs.graph.sys.mix_choice, peqs.graph.ctx,
        )
    )
    @test isequal(pnode.drift, raw)

    hn = NLevelSpace(:atom, 2)
    σ(i, j) = Transition(hn, :σ, i, j)
    H = Ω * (σ(1, 2) + σ(2, 1))
    neqs = meanfield([σ(2, 2)], H, [σ(1, 2)]; rates = [γ], order = 2)
    nk = first(keys(neqs.graph.nodes))
    nnode = QC.derive(nk, neqs.graph.sys, neqs.graph.ctx)
    @test QC._may_need_ground_reduction(nnode.op_drift)
    @test QC._may_need_ground_reduction(σ(1, 2))
    @test !QC._may_need_ground_reduction(sx)

    noise_eqs = meanfield(
        [σ(2, 2)], H, [σ(1, 2)]; rates = [γ], efficiencies = [1], order = 2,
    )
    @test first(values(noise_eqs.graph.nodes)).noise !== nothing

    ground = average(σ(1, 1))
    excited = average(σ(2, 2))
    shared = Symbolics.Num(ground + 2 * ground * excited + excited^2)
    tasks = [Threads.@spawn QC._reduce_ground_in_drift_threadsafe(shared) for _ in 1:16]
    vals = fetch.(tasks)
    expected = QC._reduce_ground_in_drift(shared)
    @test !isequal(expected, shared)
    @test all(x -> isequal(x, expected), vals)
end
