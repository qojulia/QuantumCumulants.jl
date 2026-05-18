# PENDING: port of master test/test_indexed_mixed_order.jl
#
# Status: blocked on `evaluate(eqs; limits = (N => k))` and indexed `scale`
# for `CorrelationFunction`. Both are CHANGELOG-deferred.
#
# To enable: implement evaluate + scale-on-correlation, then uncomment.

#=
using QuantumCumulants
using Symbolics: Symbolics, @variables
using SymbolicUtils
using Test

@testset "mixed_order_indexing: 1-atom-cavity collective closure" begin
    ha = NLevelSpace(:atom, 2)
    hf = FockSpace(:cavity)
    h = ha ⊗ hf

    a = Destroy(h, :a)
    σ(α, β, i) = IndexedOperator(Transition(h, :σ, α, β), i)

    @variables η::Real Γ::Real κ::Real Δc::Real N::Real ξ::Real ν::Real
    g(i) = IndexedVariable(:g, i)
    Δa(i) = IndexedVariable(:Δa, i)

    i = Index(h, :i, N, ha)
    j = Index(h, :j, N, ha)

    H = -Δc * a' * a - ∑(Δa(i) * σ(2, 2, i), i) +
        ∑(g(i) * (σ(2, 1, i) * a + σ(1, 2, i) * a'), i) +
        1im * η * (a' - a)

    J = [σ(1, 2, i), a, a' * a, σ(2, 2, i)]
    rates = [Γ, κ, ξ, ν]

    eqs = meanfield(a * σ(2, 2, j), H, J; rates = rates, order = [1, 2])
    complete!(eqs)
    @test eqs.order == [1, 2]
    @test length(eqs.equations) == 8

    evaled = evaluate(eqs; limits = (N => 3))
    @test length(evaled.equations) == 18
    @test isempty(find_missing(evaled; get_adjoints = false))
    @test evaled.order == [1, 2]

    SQA = QuantumCumulants.SecondQuantizedAlgebra
    for state in evaled.states
        aon = SQA.acts_on(state)
        term = SymbolicUtils.arguments(state)[1]
        len = term isa SQA.QMul ? length(term.args_nc) : 1
        if aon isa Vector
            @test len == 2
        elseif aon == 1
            @test len == 1
        else
            @test len <= 2
        end
    end

    corr = CorrelationFunction(a', a, eqs; steady_state = true)
    corr_sc = scale(corr)
    corr_ev = evaluate(corr; limits = Dict(N => 3))
    @test corr_ev.de.equations == evaluate(corr; limits = (N => 3)).de.equations
end
=#
