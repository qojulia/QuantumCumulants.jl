using QuantumCumulants
using SymbolicUtils
using Test

@testset "meanfield" begin

hf = FockSpace(:cavity)
ha = NLevelSpace(:atom,(:g,:e))
h = hf⊗ha

a = Destroy(h,:a)
σ = Transition(h,:σ,:g,:e)
σee = Transition(h,:σ,:e,:e,2)

# JC
(ωc,ωa,g) = (1.1341,0.4321,2.15013)
H = ωc*a'*a + ωa*σ'*σ + g*(a'*σ + σ'*a)

da = commutator(im*H,a)
@test iszero(simplify(da - (-im*ωc*a + (-im*g)*σ)))
ds = commutator(im*H,σ)
@test iszero(simplify(ds - ((-im*g)*a + (-im*ωa)*σ + (2im*g)*a*σee)))
dn = commutator(im*H,a'*a)
@test iszero(simplify(dn - ((-im*g)*a'*σ + (im*g)*a*σ')))

he = meanfield([a,σ,a'*a],H)
@test iszero(simplify(he.operator_equations[1].rhs - (da)))
@test iszero(simplify(he.operator_equations[2].rhs - (ds)))
@test iszero(simplify(he.operator_equations[3].rhs - (dn)))

# Lossy JC
κ,γ = 3.333,0.1313131313
J = [a,σ]
he_diss = meanfield([a,σ,σ'*σ],H,J;rates=[κ,γ])

@test iszero(simplify(he_diss.operator_equations[1].rhs - ((-im*ωc - 0.5κ)*a + (-im*g)*σ)))
@test iszero(simplify(he_diss.operator_equations[2].rhs - ((-im*g)*a + (-im*ωa - 0.5γ)*σ + (2im*g)*a*σee)))
@test iszero(simplify(he_diss.operator_equations[3].rhs - ((-γ)*σee + (im*g)*a'*σ + (-im*g)*a*σ')))

# Single-atom laser
ν = 3.44444444
J = [a,σ,σ']
he_laser = meanfield([a'*a,σ'*σ,a'*σ],H,J;rates=[κ,γ,ν])

@test iszero(simplify(he_laser.operator_equations[1].rhs - ((-κ)*a'*a + (-im*g)*a'*σ + (im*g)*a*σ')))
@test iszero(simplify(he_laser.operator_equations[2].rhs - (ν + (-ν - γ)*σee + (im*g)*a'*σ + (-im*g)*a*σ')))
ex = (im*g)*σee + (-im*g)*a'*a + (im*(ωc - ωa) - 0.5*(κ + γ + ν))*a'*σ + (2im*g)*a'*a*σee
@test iszero(simplify(he_laser.operator_equations[3].rhs - ex))

end # testset
