using Qumulants
using Test

@testset "mixed-order" begin

# Hilbertspace
hf = FockSpace(:cavity)
ha = NLevelSpace(:atom, 2)
h = hf⊗ha

# Operators
a = Destroy(h,:a,1)
σ(i,j) = Transition(h,:σ,i,j)

@test isequal(average(a'*a, [2,1]), average(a'*a))
@test isequal(average(a'*σ(1,2), [2,1]; mix_choice=minimum), average(a')*average(σ(1,2)))

he = heisenberg([a'*a,σ(2,2)],a'*a + σ(2,2) + a'*σ(1,2) + a*σ(2,1))
he_avg1 = average(he,2)
he_avg2 = average(he,[2,1])
he_avg3 = average(he,[2,1];mix_choice=minimum)

@test isequal(he_avg1, he_avg2)
@test !isequal(he_avg1, he_avg3)

he = heisenberg(a'*σ(1,2), a'*a + σ(2,2) + a*σ(2,1))
@test_throws ErrorException average(he,[2,1];mix_choice=minimum)

end # testset
