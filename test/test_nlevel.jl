using Qumulants
using Test
using SymbolicUtils

ha = NLevelSpace(:atom, (:g,:e))
@test ha.GS == :g

σ = Transition(ha, :σ, :g,:e)
@test σ' == Transition(ha, :σ, :e, :g)
σee = Transition(ha, :σ, :e,:e)
@test σee'==σee

ex = σ'*σ
@test simplify_operators(ex) == σee

ex = σ*σ'
σgg = Transition(ha, :σ,:g,:g)
@test simplify_operators(ex)==simplify_operators(σgg)==simplify_operators(one(σgg) - σee)


sz = σ'*σ - σ*σ'
@test simplify_operators(sz)==simplify_operators(2*σee - one(σgg))
