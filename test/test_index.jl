using Qumulants
using SymbolicUtils
using Test
SymbolicUtils.show_simplified[] = false

N = 10
I = IndexSet(:I, 1:N)
i = Index(2, I)
@test i==I[2]

@parameters p

p[i]

simplify_constants(p[i] + p[i])

p[i,i]

[p[j] for j in I]
