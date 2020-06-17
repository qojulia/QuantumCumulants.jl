using Qumulants
using Test

@testset "heisenberg" begin

hf = FockSpace(:cavity)
ha = NLevelSpace(:atom,(:g,:e))
h = hf⊗ha

a_ = Destroy(hf,:a)
a = embed(h,a_,1)
σ_ = Transition(ha,:σ,:g,:e)
σ = embed(h,σ_,2)
σee = embed(h,Transition(ha,:σ,:e,:e),2)

# JC
(ωc,ωa,g) = (1.1341,0.4321,2.15013)
H = ωc*a'*a + ωa*σ'*σ + g*(a'*σ + σ'*a)

da = commutator(1.0im*H,a)
@test iszero(simplify_operators(da - ( -1.0im*ωc*a + (-1.0im*g)*σ )))
ds = commutator(1.0im*H,σ)
@test iszero(simplify_operators(ds - ( (-1.0im*g)*a + (-1.0im*ωa)*σ + (2.0im*g)*a*σee )))
dn = commutator(1.0im*H,a'*a)
@test iszero(simplify_operators(dn - ( (-1.0im*g)*a'*σ + (1.0im*g)*a*σ' )))

he = heisenberg([a,σ,a'*a],H)
@test iszero(simplify_operators(he.rhs[1] - da))
@test iszero(simplify_operators(he.rhs[2] - ds))
@test iszero(simplify_operators(he.rhs[3] - dn))

# Lossy JC
κ,γ = 3.333,0.1313131313
J = [a,σ]
he_diss = heisenberg([a,σ,σ'*σ],H,J;rates=[κ,γ])

@test iszero(simplify_operators(he_diss.rhs[1] - ( (-1.0im*ωc - 0.5κ)*a + (-1.0im*g)*σ )))
@test iszero(simplify_operators(he_diss.rhs[2] - ( (-1.0im*g)*a + (-1.0im*ωa - 0.5γ)*σ + (2.0im*g)*a*σee )))
@test iszero(simplify_operators(he_diss.rhs[3] - ( (-γ)*σee + (1.0im*g)*a'*σ + (-1.0im*g)*a*σ' )))

# Single-atom laser
ν = 3.44444444
J = [a,σ,σ']
he_laser = heisenberg([a'*a,σ'*σ,a'*σ],H,J;rates=[κ,γ,ν])

ex = (-κ)*a'*a + (-1.0im*g)*a'*σ + (1.0im*g)*a*σ'
@test iszero(simplify_operators(he_laser.rhs[1] - ex))
ex = ν + (-ν - γ)*σee + (1.0im*g)*a'*σ + (-1.0im*g)*a*σ'
@test iszero(simplify_operators(he_laser.rhs[2] - ex))
ex = (1.0im*g)*σee + (-1.0im*g)*a'*a + (1.0im*(ωc - ωa) - 0.5*(κ + γ + ν))*a'*σ + (2.0im*g)*a'*a*σee
@test iszero(simplify_operators(he_laser.rhs[3] - ex))

end # testset
