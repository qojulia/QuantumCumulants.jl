using Qumulants
using Test

@testset "average" begin

hf = FockSpace(:cavity)
ha = NLevelSpace(:atom,(:g,:e))
h = hf⊗ha

a_ = Destroy(hf,:a)
a = embed(h, a_, 1)
σ_ = Transition(ha,:σ,:g,:e)
σ = embed(h, σ_, 2)

@test average(2a)==2*average(a)
avg_sym = Qumulants._to_symbolic(average(a))
@test Qumulants._to_qumulants(avg_sym)==average(a)

@test simplify_constants(average(σ)+average(σ))==average(2σ)
@test simplify_constants(average(σ_)-average(σ_))==0

ωc, ωa = parameters("ω_c ω_a")
@test average(ωc)==ωc
@test average(ωc*a)==ωc*average(a)

n = average(a'*a)
@test cumulant_expansion(n,2)==n
@test cumulant_expansion(n,1)==average(a')*average(a)

end # testset
