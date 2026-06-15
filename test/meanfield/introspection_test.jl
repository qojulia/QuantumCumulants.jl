using QuantumCumulants
using Symbolics: @variables
using ModelingToolkitBase: ModelingToolkitBase, unknowns, System
using Test

@testset "introspection: accessors agree with the derived view" begin
    h = FockSpace(:c)
    a = Destroy(h, :a)
    @variables Δ::Real κ::Real
    H = Δ * a' * a
    eqs = complete(meanfield([a, a' * a], H, [a]; rates = [κ], order = 2))

    @test states(eqs) === eqs.states
    @test operators(eqs) === eqs.operators

    # states[i] is the average of operators[i], by position.
    @test all(isequal(states(eqs)[i], average(operators(eqs)[i])) for i in eachindex(states(eqs)))

    # moments is the operator -> average correspondence in tracked order.
    m = moments(eqs)
    @test collect(keys(m)) == operators(eqs)
    @test all(isequal(m[op], average(op)) for op in operators(eqs))
end

@testset "introspection: moment_variable_map matches the System unknowns" begin
    h = FockSpace(:c)
    a = Destroy(h, :a)
    @variables Δ::Real κ::Real
    H = Δ * a' * a
    eqs = complete(meanfield([a, a' * a], H, [a]; rates = [κ], order = 2))

    mvm = moment_variable_map(eqs)
    @test all(isequal.(collect(keys(mvm)), states(eqs)))

    # Every mapped variable is an unknown of the built System, and the count matches.
    sys = System(eqs; name = :cav)
    u = Set(unknowns(sys))
    @test length(mvm) == length(u)
    @test all(v in u for v in values(mvm))
end

@testset "introspection: closure_report tracks open vs closed systems" begin
    h = FockSpace(:c)
    a = Destroy(h, :a)
    @variables Δ::Real κ::Real
    H = Δ * a' * a

    open_eqs = meanfield([a, a' * a], H, [a]; rates = [κ], order = 2)
    rep_open = closure_report(open_eqs)
    @test rep_open.closed == isempty(rep_open.missing)
    @test rep_open.n_states == length(states(open_eqs))
    @test sum(values(rep_open.by_order)) == rep_open.n_states

    closed_eqs = complete(open_eqs)
    rep_closed = closure_report(closed_eqs)
    @test rep_closed.closed
    @test isempty(rep_closed.missing)
    # closing only adds moments.
    @test rep_closed.n_states >= rep_open.n_states
end

@testset "introspection: noise_channels lists monitored channels" begin
    h = FockSpace(:c)
    a = Destroy(h, :a)
    @variables Δ::Real κ::Real η::Real
    H = Δ * a' * a

    det_eqs = meanfield([a, a' * a], H, [a]; rates = [κ], order = 2)
    @test isempty(noise_channels(det_eqs))

    # Two jumps, only the second monitored.
    b = Destroy(FockSpace(:c1) ⊗ FockSpace(:c2), :b, 2)
    hc2 = FockSpace(:c1) ⊗ FockSpace(:c2)
    a2 = Destroy(hc2, :a, 1)
    b2 = Destroy(hc2, :b, 2)
    @variables κa::Real κb::Real ηb::Real
    H2 = Δ * a2' * a2
    noise_eqs = meanfield(
        [a2, b2, a2' * a2, b2' * b2], H2, [a2, b2];
        rates = [κa, κb], efficiencies = [0, ηb], order = 2,
    )
    ch = noise_channels(noise_eqs)
    @test length(ch) == 1
    @test ch[1].index == 2
    @test isequal(ch[1].rate, κb)
    @test isequal(ch[1].efficiency, ηb)
end
