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

ex = σ'*σ
@test qsimplify(ex) == σee

ex = σ*σ'
σgg = Transition(ha, :σ,:g,:g)
@test isequal(qsimplify(ex), qsimplify(σgg))
@test isequal(qsimplify(σgg), (one(σgg) + -1*σee))


sz = σ'*σ - σ*σ'
@test isequal(qsimplify(sz),qsimplify(2*σee - one(σgg)))

# Integer levels
ha = NLevelSpace(:atom, 2)
@test ha.GS == 1

@test_throws AssertionError Transition(ha, :σ, :g,:e)
σ = Transition(ha, :σ, 1, 2)
@test σ' == Transition(ha, :σ, 2,1)
σee = Transition(ha, :σ, :2,2)
@test σee'==σee

ex = σ'*σ
@test qsimplify(ex) == σee

ex = σ*σ'
σgg = Transition(ha, :σ,1,1)
@test isequal(qsimplify(ex), qsimplify(σgg))
@test isequal(qsimplify(σgg), (one(σgg) + -1*σee))

sz = σ'*σ - σ*σ'
@test isequal(qsimplify(sz),qsimplify(2*σee - one(σgg)))

# Product space
ha1 = NLevelSpace(:atom, (:g,:e))
ha2 = NLevelSpace(:atom, 2)
@test ha1 != ha2

hprod = ha1⊗ha2
σ1 = embed(hprod,Transition(ha1,:σ,:g,:e),1)
σ2 = embed(hprod,Transition(ha2,:σ,1,2),2)
@test qsimplify(σ1'*σ1)==embed(hprod,Transition(ha1,:σ,:e,:e),1)
@test isequal(qsimplify(σ2*σ2'), qsimplify(1 -embed(hprod,Transition(ha2,:σ,2,2),2)))
@test isequal(qsimplify(σ1*σ2), σ1*σ2)

@test_throws ErrorException Transition(hprod,:σ,:g,:e)
@test Transition(hprod,:σ,:g,:e,1)==σ1
@test Transition(hprod,:σ,1,2,2)==σ2
@test_throws AssertionError Transition(hprod,:σ,1,2,1)

end # testset
