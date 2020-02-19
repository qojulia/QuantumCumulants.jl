using Qumulants
using SymPy

a = Destroy(:a)⊗Identity()
σ(i,j,k) = Identity()⊗Transition(:σ,i,j,2)[k]
i = Index(:i,1,:N)
j = Index(:j,1,:N)

g = symbols("g",real=true)
Δ = symbols("Δ",real=true)
Γ = symbols("Γ",real=true,positive=true)
κ = symbols("κ",real=true,positive=true);

# Hamiltonian

#cavity
Hc = -Δ*(a'*a)
#atom-cavity
Hac = g[i]*(a'*σ(1,2,i) + a*σ(2,1,i))

H = Hc + Sum(Hac,i)

# Jumps
J = [a,σ(1,2,i)]
rates = [κ, Γ[i]];

# Derive q. Langevin eq. (without noise)
ops = [a,σ(1,2,j), a'*a, a'*σ(1,2,j), σ(2,2,j)] #operators for which the equations should be derived
me = master(ops,H,J;rates=rates)
me_avg = average(me,2)

isym = sympify(i)
tmp = SymPy.sympy.Sum(g[i], (isym, isym.__pyobject__.lower, isym.__pyobject__.upper))


k = Index(:k, 1, :N; neq=[i])
average(σ(2,1,i)*σ(1,2,k))
