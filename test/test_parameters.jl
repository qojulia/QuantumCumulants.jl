using Qumulants
using Test

@testset "parameters" begin

hf = FockSpace(:cavity)
ha = NLevelSpace(:atom,(:g,:e))
h = hf⊗ha

a_ = Destroy(hf,:a)

@parameters p q
@test p*q isa NumberTerm
@test p*q*a_ isa OperatorTerm

a = embed(h,a_,1)
σ_ = Transition(ha,:σ,:g,:e)
σ = embed(h,σ_,2)
σee = embed(h,Transition(ha,:σ,:e,:e),2)

# JC
@parameters ωc ωa g
H = ωc*a'*a + ωa*σ'*σ + g*(a'*σ + σ'*a)

da = commutator(1.0im*H,a)
@test da == simplify_operators((-1.0im*g)*σ + (-1.0im*ωc)*a)
ds = commutator(1.0im*H,σ)
@test ds == simplify_operators((-1.0im*g)*a + (-1.0im*ωa)*σ + (2.0im*g)*a*σee)
dn = commutator(1.0im*H,a'*a)
@test dn == simplify_operators((-1.0im*g)*a'*σ + (1.0im*g)*a*σ')

he = heisenberg([a,σ,a'*a],H)
@test he.rhs[1] == da
@test he.rhs[2] == ds
@test he.rhs[3] == dn

# Lossy JC
@parameters κ γ
J = [a,σ]
he_diss = heisenberg([a,σ,σ'*σ],H,J;rates=[κ,γ])

@test he_diss.rhs[1] == simplify_operators((-1.0im*ωc - 0.5κ)*a + (-1.0im*g)*σ)
@test he_diss.rhs[2] == simplify_operators((-1.0im*g)*a + (-1.0im*ωa - 0.5γ)*σ + (2.0im*g)*a*σee)
@test he_diss.rhs[3] == simplify_operators((-γ)*σee + (1.0im*g)*a'*σ + (-1.0im*g)*a*σ')

# Single-atom laser
@parameters ν
J = [a,σ,σ']
he_laser = heisenberg([a'*a,σ'*σ,a'*σ],H,J;rates=[κ,γ,ν])

@test he_laser.rhs[1] == simplify_operators((-κ)*a'*a + (-1.0im*g)*a'*σ + (1.0im*g)*a*σ')
@test he_laser.rhs[2] == simplify_operators((1.0im*g)*a'*σ + (-1.0im*g)*a*σ' + ν + (-ν - γ)*σee)
@test he_laser.rhs[3] == simplify_operators((1.0im*g)*σee + (-1.0im*g)*a'*a + (1.0im*(ωc - ωa) - 0.5*(κ + γ + ν))*a'*σ + (2.0im*g)*a'*a*σee)

end # testset
