using Qumulants
using Test

@testset "average" begin

hf = FockSpace(:cavity)
ha = NLevelSpace(:atom,(:g,:e))
h = hf⊗ha

a = Destroy(h,:a)
σ = Transition(h,:σ,:g,:e)

@test isequal(average(2a),2*average(a))
avg_sym = Qumulants._to_symbolic(average(a))
@test isequal(Qumulants._to_qumulants(avg_sym),average(a))

@test isequal(average(2im*a'*σ)' , (-2im)*average(a*σ'))

@test isequal(simplify_constants(average(σ)+average(σ)),average(2σ))
@test simplify_constants(average(σ)-average(σ))==0

ωc, ωa = parameters("ω_c ω_a")
@test isequal(average(ωc), ωc)
@test isequal(average(ωc*a), ωc*average(a))

n = average(a'*a)
@test isequal(cumulant_expansion(n,2), n)
@test isequal(cumulant_expansion(n,1), average(a')*average(a))

ex = average(a*σ)*average(a) - average(a)*average(a*σ)
@test isequal(simplify_constants(ex), 0)

end # testset
