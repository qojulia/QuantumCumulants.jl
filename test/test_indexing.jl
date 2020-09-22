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
n2 = (j!=k)*Qumulants.nip(σ(2,1)[j],σ(1,2)[k])

@test simplify_operators(σ(2,1)[i]*σ(1,2)[i]) == σ(2,2)[i]
@test simplify_operators(σ(2,1)[i]*σ(1,2)[j]) == n + (i==j)*σ(2,2)[i]
@test simplify_operators(n*σ(2,1)[i])==simplify_operators(n*σ(1,2)[j])==0
@test simplify_operators(n*σ(1,2)[i])==(i!=j)*Qumulants.nip(σ(1,2)[j],σ(2,2)[i])

@test simplify_operators(n*σ(2,1)[i])==simplify_operators(n*σ(1,2)[j])==0
@test simplify_operators(n*σ(1,2)[i])==(i!=j)*Qumulants.nip(σ(1,2)[j],σ(2,2)[i])

hf = FockSpace(:field)
ha = NLevelSpace(:atom,2)
h = hf⊗ha
a = Destroy(h,:a)
σ(i,j) = Transition(h,:σ,i,j)
@parameters g
H = Sum(g[i]*(a'*σ(1,2)[i] + a*σ(2,1)[i]), i)

n = (j!=k)*Qumulants.nip(σ(2,1)[k],σ(1,2)[j])
ops = [a,σ(1,2)[j],σ(2,2)[j],a'*σ(1,2)[j],n]
he = heisenberg(ops,H)

H1 = g*(a'*σ(1,2) + a*σ(2,1))
ops1 = [a,σ(1,2),σ(2,2),a'*σ(1,2)]
he1 = heisenberg(ops1, H1)
