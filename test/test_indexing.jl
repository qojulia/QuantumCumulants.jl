using Qumulants
using Test

i = Index(:i,1,:n)
j = Index(:j,1,:n)
k = Index(:k,1,:n)

h = FockSpace(:fock)
a = Destroy(h,:a)


@test isequal(simplify_operators(a[i]*a[j]'), (i==j) + a[j]'*a[i])
@test isequal(simplify_operators(a[i]*a[i]'), true+a[i]'*a[i])


# Nlevel
h = NLevelSpace(:nlevel,2)
σ(i,j) = Transition(h,:σ,i,j)
n = (i!=j)*Qumulants.nip(σ(2,1)[i],σ(1,2)[j])


@test simplify_operators(σ(2,1)[i]*σ(1,2)[i]) == σ(2,2)[i]
@test simplify_operators(σ(2,1)[i]*σ(1,2)[j]) == n + (i==j)*σ(2,2)[i]
@test simplify_operators(n*σ(2,1)[i])==simplify_operators(n*σ(1,2)[j])==0
@test simplify_operators(n*σ(1,2)[i])==(i!=j)*Qumulants.nip(σ(2,2)[i],σ(1,2)[j])

@test simplify_operators(n*σ(2,1)[i])==simplify_operators(n*σ(1,2)[j])==0
@test simplify_operators(n*σ(1,2)[i])==(i!=j)*Qumulants.nip(σ(2,2)[i],σ(1,2)[j])



n = Qumulants.nip(σ(2,1)[i],σ(3,1)[j])

ex = Qumulants._to_symbolic(σ(2,1)[i]*n)
ex2 = Qumulants._to_symbolic(σ(2,))
