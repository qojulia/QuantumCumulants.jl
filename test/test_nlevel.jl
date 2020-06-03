using Qumulants
using Test

@testset "nlevel" begin

# Symbolic levels
ha = NLevelSpace(:atom, (:g,:e))
@test ha.GS == :g

σ = Transition(ha, :σ, :g,:e)
@test σ' == Transition(ha, :σ, :e, :g)
σee = Transition(ha, :σ, :e,:e)
@test σee'==σee

σ_sym = Qumulants._to_symbolic(σ)
@test Qumulants._to_qumulants(σ_sym)==σ

ex = σ'*σ
@test simplify_operators(ex) == σee

ex = σ*σ'
σgg = Transition(ha, :σ,:g,:g)
@test simplify_operators(ex)==simplify_operators(σgg)==simplify_operators(one(σgg) - σee)


sz = σ'*σ - σ*σ'
@test simplify_operators(sz)==simplify_operators(2*σee - one(σgg))

# Integer levels
ha = NLevelSpace(:atom, 2)
@test ha.GS == 1

@test_throws AssertionError Transition(ha, :σ, :g,:e)
σ = Transition(ha, :σ, 1, 2)
@test σ' == Transition(ha, :σ, 2,1)
σee = Transition(ha, :σ, :2,2)
@test σee'==σee

ex = σ'*σ
@test simplify_operators(ex) == σee

ex = σ*σ'
σgg = Transition(ha, :σ,1,1)
@test simplify_operators(ex)==simplify_operators(σgg)==simplify_operators(one(σgg) - σee)


sz = σ'*σ - σ*σ'
@test simplify_operators(sz)==simplify_operators(2*σee - one(σgg))

# Product space
ha1 = NLevelSpace(:atom, (:g,:e))
ha2 = NLevelSpace(:atom, 2)
@test ha1 != ha2

hprod = ha1⊗ha2
σ1 = embed(hprod,Transition(ha1,:σ,:g,:e),1)
@test Qumulants._to_qumulants(Qumulants._to_symbolic(σ1)) == σ1
@test Qumulants._to_qumulants(Qumulants._to_symbolic(σ1'))==σ1'
σ2 = embed(hprod,Transition(ha2,:σ,1,2),2)
@test Qumulants._to_qumulants(Qumulants._to_symbolic(σ2)) == σ2
@test simplify_operators(σ1'*σ1)==embed(hprod,Transition(ha1,:σ,:e,:e),1)
@test simplify_operators(σ2*σ2')==simplify_operators(1 -embed(hprod,Transition(ha2,:σ,2,2),2))
@test simplify_operators(σ1*σ2)==σ1*σ2

@test_throws ErrorException Transition(hprod,:σ,:g,:e)
@test Transition(hprod,:σ,:g,:e,1)==σ1
@test Transition(hprod,:σ,1,2,2)==σ2
@test_throws AssertionError Transition(hprod,:σ,1,2,1)

end # testset
