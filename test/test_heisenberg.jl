using Qumulants
using Test
using SymbolicUtils

hf = FockSpace(:cavity)
ha = NLevelSpace(:atom,(:g,:e))

a_ = Destroy(hf,:a)
a = a_⊗one(ha)
σ_ = Transition(ha,:σ,:g,:e)
σ = one(hf)⊗σ_

# JC
(ωc,ωa,g) = (1.1341,0.4321,2.15013)
H = ωc*a'*a + ωa*σ'*σ + g*(a'*σ + σ'*a)

da = commutator(1.0im*H,a)
@test da == -1.0im*ωc*a + (-1.0im*g)*σ
ds = commutator(1.0im*H,σ)
@test ds == (-1.0im*g)*a + (2.0im*g)*a_⊗Transition(ha,:σ,:e,:e) + (-1.0im*ωa)*σ

he = heisenberg([a,σ],H)
@test he.rhs[1] == da
@test he.rhs[2] == ds

# Lossy JC
κ,γ = 3.333,0.1313131313
J = [a,σ]
he_diss = heisenberg([a,σ,σ'*σ],H,J;rates=[κ,γ])

@test he_diss.rhs[1] == (-1.0im*ωc - 0.5κ)*a + (-1.0im*g)*σ
@test he_diss.rhs[2] == (-1.0im*g)*a + (2.0im*g)*a_⊗Transition(ha,:σ,:e,:e) + (-1.0im*ωa - 0.5γ)*σ
@test he_diss.rhs[3] == (1.0im*g*a_')⊗σ_ + (-1.0im*g*a_)⊗σ_' + (-γ*one(a_))⊗Transition(ha,:σ,:e,:e)

# Single-atom laser
ν = 3.44444444
J = [a,σ,σ']
he_laser = heisenberg([a'*a,σ'*σ,a'*σ],H,J;rates=[κ,γ,ν])
