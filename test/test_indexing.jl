using Qumulants
using Test

hf = FockSpace(:fock)

n = Parameter{Int}(:n)
i = Index(:i, 1, n)
j = Index(:j, 1, n)
a(i) = Destroy(hf,:a; index=i)

@test simplify_operators(a(i)+a(i))==2*a(i)
@test simplify_operators(a(i) + a(j)) == simplify_operators(a(j) + a(i))
@test simplify_operators(a(i)*a(j)) == simplify_operators(a(j)*a(i))

@test simplify_operators(a(i)*a(i)') == 1 + a(i)'*a(i)
@test simplify_operators(a(i)*a(j)') == (i==j) + a(j)'*a(i)
