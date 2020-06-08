using Qumulants
using SymbolicUtils

# Hilbert space
N = 2
hf = FockSpace(:field)
ha = NLevelSpace(:atom,(:g,:e))
# hN = ncopy(ha, N)
h_tot = hf⊗ha

# Operators
a = Destroy(h_tot,:a)
σ(i,j,k) = k isa Index ? Transition(h_tot, :σ, i, j; index=k) : Transition(h_tot, :σ, i, j; index=I[k])
σ(k) = σ(:g,:e,k)

# Indices
I = IndexSet(:I, 1, N)

H = +((σ(i)'*σ(i) for i=I)...) + +((σ(i) + σ(i)' for i=I)...)
he = heisenberg([σ(I[1]), σ(:e,:e,I[1])], H)

# Parameters
@parameters Δ Ω
H = +((Δ[i]*σ(i)'*σ(i) for i=I)...) + +((Ω[i]*(σ(i) + σ(i)') for i=I)...)
he = heisenberg([σ(I[1]), σ(:e,:e,I[1])], H)

# Many-atom laser
@parameters Δ g κ γ ν

H = +((Δ[i]*σ(i)'*σ(i) for i=I)...) + +((g[i]*(a'*σ(i) + a*σ(i)') for i=I)...)
J = [a;[σ(i) for i=I];[σ(i)' for i=I]]
rates = [κ;[γ[i] for i=I];[ν[i] for i=I]]

# he = heisenberg([a'*a,a'*σ(I[1]),σ(:e,:e,I[1])],H,J;rates=rates)
he = heisenberg(σ(1)'*σ(2),H)

simplify_operators(σ(1)'*σ(2)*H)
ex = Qumulants._to_symbolic(σ(1)'*σ(2)*H)
ex2 = simplify(ex;rules=Qumulants.SIMPLIFY_COMMUTATOR_RULES,fixpoint=false)
for i=1:2
    global ex2 = simplify(ex2;rules=Qumulants.SIMPLIFY_COMMUTATOR_RULES,fixpoint=false)
end
ex3 = SymbolicUtils.flatten_term(*, ex2.arguments)

ex4 = Qumulants._to_symbolic(g[I[2]]*a*σ(I[1])*(1 + (-1*σ(:e,:e,I[2]))))
