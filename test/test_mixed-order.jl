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

# N-atom laser
# Parameters
N = 2 #number of atoms
κ, g, Γ23, Γ13, Γ12, Ω, Δc, Δ3 = params("κ g Γ_{23} Γ_{13} Γ_{12} Ω Δ_c Δ_3")

# Hilbertspace
hf = FockSpace(:cavity)
ha = ⊗([NLevelSpace(Symbol(:atom,i),3) for i=1:N]...)
h = hf ⊗ ha

# Operators
a = Destroy(h,:a)
σ(i,j,k) = Transition(h,Symbol("σ_{$k}"),i,j,k+1)

# Hamiltonian
H = -Δc*a'a + sum(g*(a'*σ(1,2,i) + a*σ(2,1,i)) for i=1:N) + sum(Ω*(σ(3,1,i) + σ(1,3,i)) for i=1:N) - sum(Δ3*σ(3,3,i) for i=1:N)

# Jumps
J = [a;[σ(1,2,i) for i=1:N];[σ(1,3,i) for i=1:N];[σ(2,3,i) for i=1:N]]

# Rates
rates = [κ;[Γ12 for i=1:N];[Γ13 for i=1:N];[Γ23 for i=1:N]]

# list of operators
ops = [a'a, σ(2,2,1), σ(3,3,1)]

he = heisenberg(ops,H,J; rates=rates)
he_avg_ = average(he,[2,1,1]) #second order average

he_avg = complete(he_avg_;order=[2,1,1])
@test isempty(find_missing(he_avg))
@test isempty(findall(x -> (aon = acts_on(x); 2 in aon && 3 in aon), he_avg.lhs))

he_avg_ = average(he,[2,1,1];mix_choice=minimum)
he_avg = complete(he_avg_;order=[2,1,1],mix_choice=minimum)
@test isempty(find_missing(he_avg))
@test isempty(findall(x -> (aon = acts_on(x); 2 in aon && 3 in aon), he_avg.lhs))
@test isempty(findall(x -> (aon = acts_on(x); 1 in aon && 3 in aon), he_avg.lhs))
@test isempty(findall(x -> (aon = acts_on(x); 1 in aon && 2 in aon), he_avg.lhs))

end # testset
