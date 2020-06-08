using Qumulants
using SymbolicUtils

# Hilbert space
N = 4
hf = FockSpace(:field)
ha = NLevelSpace(:atom,(:g,:e))
# hN = ncopy(ha, N)
h_tot = hf⊗ha

# Operators
a = Destroy(h_tot,:a)
σ(i,j,k) = Transition(h_tot, :σ, i, j; index=k)
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

he = heisenberg([a'*a,a'*σ(I[1]),σ(:e,:e,I[1])],H,J;rates=rates)
