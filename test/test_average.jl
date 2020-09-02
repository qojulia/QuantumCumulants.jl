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
@test average(2*(a+a)) == 2*(average(a) + average(a))
@test average(a^2) == average(a*a)
@test average(a^2,1) == average(a)^2
@test average((a'*a)^2) == average(a'*a*a'*a)
@test average(a'*a^2,2) == average(a')*average(a^2) + -2*average(a')*average(a)^2 + 2*average(a)*average(a'*a)

@test simplify_constants(average(σ)+average(σ))==average(2σ)
@test simplify_constants(average(σ)-average(σ))==0

ωc, ωa = parameters("ω_c ω_a")
@test average(ωc)==ωc
@test average(ωc*a)==ωc*average(a)
@test average(ωc*(a+a')) == ωc*(average(a) + average(a'))

n = average(a'*a)
@test cumulant_expansion(n,2)==n
@test cumulant_expansion(n,1)==average(a')*average(a)

ex = average(a*σ)*average(a) - average(a)*average(a*σ)
@test simplify_constants(ex)==0

# Test cumulants
hs = FockSpace[]
for i=1:4
    push!(hs, FockSpace(Symbol(:fock, i)))
end
h = ⊗(hs...)

a = Destroy(h,:a,1)
b = Destroy(h,:b,2)
c = Destroy(h,:c,3)
d = Destroy(h,:d,4)

@test cumulant(a*b)==average(a*b)+ -1*average(a)*average(b)
@test cumulant(a*b,1) == average(a*b)
@test iszero(simplify_constants(expand(cumulant(a*b*c) - (average(a*b*c) +
            2*average(a)*average(b)*average(c) - average(a)*average(b*c) -
            average(b)*average(a*c) - average(c)*average(a*b)))))
@test simplify_constants(expand(average(a*b*c*d) - cumulant(a*b*c*d))) == cumulant_expansion(average(a*b*c*d),3)

end # testset
