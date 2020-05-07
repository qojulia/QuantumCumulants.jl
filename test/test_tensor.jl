using Qumulants
using Test
using SymbolicUtils

hf = FockSpace(:cavity)
ha = NLevelSpace(:atom,(:g,:e))

a_ = Destroy(hf,:a)
a = a_⊗one(ha)
σ_ = Transition(ha,:σ,:g,:e)
σ = one(hf)⊗σ_

ex = a^2
@test simplify_operators(ex) == a_^2 ⊗ one(ha)

# JC
(ωc,ωa,g) = (1.1341,0.4321,2.15013)
H = ωc*a'*a + ωa*σ'*σ + g*(a'*σ + σ'*a)

commutator(a,b) = simplify_operators(a*b - b*a)
da = commutator(1.0im*H,a)
ds = commutator(1.0im*H,σ)
