using Qumulants
using Test

i = Index(:i,1,:n)
j = Index(:j,1,:n)

h = FockSpace(:fock)
ai = IndexedDestroy(h,:a,1,i)
aj = IndexedDestroy(h,:a,1,j)

@test isequal(simplify_operators(ai*aj'), (i==j) + aj'*ai)
@test isequal(simplify_operators(ai*ai'), true+ai'*ai)
