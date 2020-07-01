using Qumulants
using Test

hf = FockSpace(:fock)

n = Parameter{Int}(:n)
i = Index(:i, 1, n)
j = Index(:j, 1, n)
k = Index(:k, 1, n)
a(i) = Destroy(hf,:a; index=i)

@test simplify_operators(a(i)+a(i))==2*a(i)
@test simplify_operators(a(i) + a(j)) == simplify_operators(a(j) + a(i)) == a(i)+a(j)
@test simplify_operators(a(i)*a(j)) == simplify_operators(a(j)*a(i)) == a(i)*a(j)

@test simplify_operators(a(i)*a(i)') == 1 + a(i)'*a(i)
@test simplify_operators(a(i)*a(j)') == (i==j) + a(j)'*a(i)

ex = 2*(i!=j)*(k!=j)*a(i)
ex_sym = Qumulants._to_symbolic(ex)
@test isequal(Qumulants.find_neq_inds(ex), Qumulants.find_neq_inds(ex_sym))
neq_inds = Qumulants.find_neq_inds(ex_sym.arguments)
@test (i!=j) in neq_inds
@test (k!=j) in neq_inds
@test !(j!=k in neq_inds)


ha = NLevelSpace(:atom, (:g,:e))
σ(i,j,k) = Transition(ha,:σ,i,j;index=k)
σ(i) = σ(:g,:e,i)

ex = σ(i)*σ(i)
@test iszero(simplify_operators(ex))
ex = σ(i)'*σ(i)
@test isequal(simplify_operators(ex), σ(:e,:e,i))
ex = σ(i)*σ(j)
ex_ = Qumulants.neq_inds_prod(ex.arguments, [i!=j])
@test isequal(simplify_operators(ex), ex_)
ex = σ(i)'*σ(j)
