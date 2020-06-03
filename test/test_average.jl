using Qumulants
using Test

@testset "average" begin

hf = FockSpace(:cavity)
ha = NLevelSpace(:atom,(:g,:e))
h = hf⊗ha

a = Destroy(h,:a)
σ = Transition(h,:σ,:g,:e)

@test average(2a)==2*average(a)
avg_sym = Qumulants._to_symbolic(average(a))
@test Qumulants._to_qumulants(avg_sym)==average(a)

@test average(2im*a'*σ)' == (-2im)*average(a*σ')

@test simplify_constants(average(σ)+average(σ))==average(2σ)
@test simplify_constants(average(σ)-average(σ))==0

ωc, ωa = parameters("ω_c ω_a")
@test average(ωc)==ωc
@test average(ωc*a)==ωc*average(a)

n = average(a'*a)
@test cumulant_expansion(n,2)==n
@test cumulant_expansion(n,1)==average(a')*average(a)

ex = average(a*σ)*average(a) - average(a)*average(a*σ)
# @test simplify_constants(ex)==0 TODO; bug in sorting

end # testset
