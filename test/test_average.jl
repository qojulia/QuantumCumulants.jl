using QuantumCumulants
using SymbolicUtils
using Test

@testset "average" begin

hf = FockSpace(:cavity)
ha = NLevelSpace(:atom,(:g,:e))
h = hf⊗ha

a = Destroy(h,:a)
σ = Transition(h,:σ,:g,:e)

@test isequal(average(2a), 2*average(a))

@test isequal(QuantumCumulants._conj(average(2im*a'*σ)), (-2im)*average(a*σ'))
@test isequal(average(2*(a+a)), 2*(average(a) + average(a)))
@test isequal(average(a^2), average(a*a))
@test isequal(average(a^2,1), average(a)^2)
@test isequal(average((a'*a)^2), average(a'*a*a'*a))
@test isequal(average(a'*a^2,2), average(a')*average(a^2) + -2*average(a')*average(a)^2 + 2*average(a)*average(a'*a))

@test isequal(simplify(average(σ)+average(σ)), average(2σ))
@test iszero(simplify(average(σ)-average(σ)))

ωc, ωa = cnumbers("ω_c ω_a")
@test isequal(average(ωc),ωc)
@test isequal(average(ωc*a),ωc*average(a))
@test isequal(average(ωc*(a+a')) , ωc*average(a) + ωc*average(a'))

n = average(a'*a)
@test isequal(cumulant_expansion(n,2), n)
@test isequal(cumulant_expansion(n,1), average(a')*average(a))

@test iszero(average(a*σ)*average(a) - average(a)*average(a*σ))

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

@test isequal(cumulant(a*b),simplify(average(a*b)+ -1*average(a)*average(b)))
@test isequal(cumulant(a*b,1), average(a*b))
@test iszero(simplify(cumulant(a*b*c) - (average(a*b*c) +
            2*average(a)*average(b)*average(c) - average(a)*average(b*c) -
            average(b)*average(a*c) - average(c)*average(a*b))))
@test isequal(simplify(average(a*b*c*d) - cumulant(a*b*c*d)), cumulant_expansion(average(a*b*c*d),3))

# cumulant expansion test with function of moments
# https://en.wikipedia.org/wiki/Cumulant#First_several_cumulants_as_functions_of_the_moments
hf1 = FockSpace(:c1)
hf2 = FockSpace(:c2)
hf3 = FockSpace(:c3)
hf4 = FockSpace(:c4)
hf_all = hf1⊗hf2⊗hf3⊗hf4

a1 = Destroy(hf_all,:a1,1)
a2 = Destroy(hf_all,:a2,2)
a3 = Destroy(hf_all,:a3,3)
a4 = Destroy(hf_all,:a4,4)

Δa1 = a1 - average(a1)
Δa2 = a2 - average(a2)
Δa3 = a3 - average(a3)
Δa4 = a4 - average(a4)

@test isequal(average(Δa1*Δa2), average(a1*a2) - average(a1)*average(a2))
@test isequal(average(Δa1*Δa2*Δa3), cumulant(a1*a2*a3))
@test iszero(average(Δa1*Δa2*Δa3)-cumulant(a1*a2*a3))
@test iszero(expand(cumulant(a1^4) - ( average(Δa1^4) - 3*(average(Δa1^2))^2 ) ))
@test iszero(expand(cumulant(a1*a2*a3*a4) - ( average(Δa1*Δa2*Δa3*Δa4) - ( average(Δa1*Δa2)*average(Δa3*Δa4) + average(Δa1*Δa3)*average(Δa2*Δa4) + average(Δa1*Δa4)*average(Δa3*Δa2) ) ) ))

end # testset
