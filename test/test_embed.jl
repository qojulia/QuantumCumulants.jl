using Qumulants
using SymbolicUtils
SymbolicUtils.show_simplified[]=false
using Test

# Test
hf = FockSpace(:c)

a = Destroy(hf,:a)

h2 = hf⊗hf
a_embed = embed(h2,a,2)
a_sym = Qumulants._to_symbolic(a_embed)
@test Qumulants._to_qumulants(a_sym)==a_embed

A(i) = embed(h2,a,i)

@test A(1)*A(2)'==simplify_operators(A(1)*A(2)')
