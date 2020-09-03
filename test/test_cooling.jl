using Qumulants
using Test


# Hilbert space
hm = FockSpace(:motion)
ha = NLevelSpace(:atom,(:g,:e))
h = hm⊗ha

# Operators
σ(i,j) = Transition(h,:σ,i,j)
a = Destroy(h,:a)
x = (a+a')
p = -im*(a-a')

# Parameters
@parameters Δ Ω k γ m

# Hamiltonian
H = -Δ*σ(:e,:e) + Ω*cos(k*x)*(σ(:g,:e) + σ(:e,:g)) + p^2/(2m)
J = [σ(:g,:e)]
rates = [γ]

# Operators
ops = [a,σ(:g,:e),σ(:e,:e)]
he = heisenberg(ops,H,J;rates=rates)
