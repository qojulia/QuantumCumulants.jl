using Qumulants
using Test

δ(i,j) = Qumulants.KroneckerDelta(i,j)

@testset "indexing" begin

# Test single mode
a = Destroy(:a)
i = Index(:i,1,:N)
@test i ∈ Qumulants.IndexOrder
i = Index(:i,1,:N)
@test length(findall(isequal(i),Qumulants.IndexOrder))==1
j = Index(:j,1,:N)
@test length(findall(x->(x==i||x==j),Qumulants.IndexOrder))==2

# Test with symbolic index
@test a[i] == a[i]
@test a[i] != a[j]
@test simplify_operators(a[i] + a[i]) == 2*a[i]
@test simplify_operators(a[i] + a[j]) == a[i]+a[j]
@test a[i]' == getindex(a', i)
@test one(a[i]) == one(a)
@test Qumulants.apply_comms(a[i]*a[i]') == a[i]'*a[i] + one(a[i]) == Qumulants.expression2index(a'*a,i) + one(a)

@test a[i]*a[j] == a[j]*a[i]
@test a[j]*a[i]'*a[j]'*a[i] == a[i]'*a[i]*a[j]*a[j]'
@test a[1]*a[j] == a[j]*a[1]

@test simplify_operators(a[i]*a[j] + a[j]*a[i]) == 2*a[i]*a[j] == 2*a[j]*a[i]
@test simplify_operators(a[i]*a[j] - a[j]*a[i]) == zero(a[j])

# Test with numeric index
@test a[1] == a[1]
@test a[1] != a[2]
@test simplify_operators(a[1] + a[1]) == 2*a[1]
# @test simplify_operators(a[1] + a[2]) == a[1]+a[2]
@test a[1]' == getindex(a', 1)
@test one(a[1]) == one(a)[1]
@test Qumulants.apply_comms(a[2]*a[2]') == a[2]'*a[2] + one(a[2]) == Qumulants.expression2index(a'*a + one(a),2)

# Test transitions
σ(i,j) = Transition(:σ,i,j,2)

@test σ(1,2)[i] == σ(1,2)[i]
@test σ(1,2)[i] != σ(1,2)[j]
@test simplify_operators(σ(1,2)[i]*σ(1,2)[i]) == Zero()[i]
@test simplify_operators(σ(2,1)[i]*σ(1,2)[i]) == σ(2,2)[i]
@test simplify_operators(σ(1,2)[i]*σ(2,1)[i]) == Qumulants.expression2index(Identity() - σ(2,2), i)

ex = σ(2,1)[i]*σ(1,2)[j]
@test simplify_operators(ex) == σ(2,2)[i]*δ(i,j) + (1-δ(i,j))*ex
@test σ(2,1)[i]*σ(1,2)[j] == copy(ex)

# Test laser with i-th atom
a = Destroy(:a) ⊗ Identity()
s = Identity() ⊗ σ(1,2)[i]
sps = simplify_operators(s'*s)
sz = simplify_operators(s'*s - s*s')
ωc, ωa, g, γ, κ, ν = (0.1, 3, 0.5, 0.25, 1.0, 4)
H = ωc*a'*a + ωa*s'*s + g*(a'*s + a*s')
J = [sqrt(κ)*a,sqrt(γ)*s,sqrt(ν)*s']

ops = [a'*a, sps, a'*s]
eqs = heisenberg(ops,H,J)

rhs1 = simplify_operators(-κ*a'*a -1.0im*g*(a'*s - a*s'))
rhs2 = simplify_operators(-(γ + ν)*s'*s + ν*one(s) + 1.0im*g*(a'*s - a*s'))
rhs3 = simplify_operators(1.0im*(ωc - ωa)*a'*s - 0.5*(κ + ν + γ)*a'*s + 1.0im*g*sz*a'*a + 1.0im*g*s'*s)
eqs_test = DifferentialEquationSet(ops,simplify_operators.([rhs1,rhs2,rhs3]))
@test eqs == eqs_test

ex = Identity()⊗(σ(2,1)[i]*σ(1,2)[j])
simplify_operators(ex)

end # testset
