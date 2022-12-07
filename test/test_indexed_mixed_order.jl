using QuantumCumulants
using SymbolicUtils
using Test

@testset "mixed_order_indexing" begin


ha = NLevelSpace(:atom, 2)
hf = FockSpace(:cavity)
h = ha ⊗ hf

a = Destroy(h, :a)
σ(α,β,i) = IndexedOperator(Transition(h, :σ, α, β),i)

@cnumbers η Γ κ Δc N ξ ν
g(i) = IndexedVariable(:g, i)
Δa(i) = IndexedVariable(:Δa, i)

i = Index(h,:i,N,ha)
j = Index(h,:j,N,ha)

H = -Δc*a'a - ∑(Δa(i)*σ(2,2,i),i) + ∑(g(i)*(σ(2,1,i)*a + σ(1,2,i)*a'),i) + 1im*η*(a' - a)

J = [σ(1,2,i), a, a'a, σ(2,2,i)]
rates = [Γ, κ, ξ, ν]
order = [1, 2]

eqs = meanfield(a*σ(2,2,j),H,J;rates,order=order)
QuantumCumulants.complete!(eqs);

@test eqs.order == [1,2]
@test length(eqs) == 8

evaled = evaluate(eqs;limits=(N=>3))

@test length(evaled) == 18

missed = find_missing(evaled)
filter!(x->x ∉ conj.(evaled.states),missed)

@test isempty(missed)
@test evaled.order == [1,2]

for state in evaled.states
    aon = acts_on(state)
    term = arguments(state)[1]
    len = (term isa QuantumCumulants.QMul) ? length(term.args_nc) : 1
    if aon isa Vector
        @test len == 2
    elseif aon == 1
        @test len == 1
    else
        @test len <= 2
    end
end

for state in eqs.states
    aon = acts_on(state)
    term = arguments(state)[1]
    len = (term isa QuantumCumulants.QMul) ? length(term.args_nc) : 1
    if aon isa Vector
        @test len == 2
    elseif aon == 1
        @test len == 1
    else
        @test len <= 2
    end
end


end # testset
