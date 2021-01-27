using Qumulants
using Test

@testset "parameters" begin

hf = FockSpace(:cavity)
ha = NLevelSpace(:atom,(:g,:e))
h = hf⊗ha

a_ = Destroy(hf,:a)

@params p q
# @test p*q isa NumberTerm
@test p*q*a_ isa OperatorTerm

a = embed(h,a_,1)
σ_ = Transition(ha,:σ,:g,:e)
σ = embed(h,σ_,2)
σee = embed(h,Transition(ha,:σ,:e,:e),2)

# JC
@params ωc ωa g
H = ωc*a'*a + ωa*σ'*σ + g*(a'*σ + σ'*a)

da = commutator(1.0im*H,a)
@test isequal(da, simplify_operators((-1.0im*g)*σ + (-1.0im*ωc)*a))
ds = commutator(1.0im*H,σ)
@test isequal(ds, simplify_operators((-1.0im*g)*a + (-1.0im*ωa)*σ + (2.0im*g)*a*σee))
dn = commutator(1.0im*H,a'*a)
@test isequal(dn, simplify_operators((-1.0im*g)*a'*σ + (1.0im*g)*a*σ'))

he = heisenberg([a,σ,a'*a],H)
@test isequal(he.rhs[1], da)
@test isequal(he.rhs[2], ds)
@test isequal(he.rhs[3], dn)

# Lossy JC
@params κ γ
J = [a,σ]
he_diss = heisenberg([a,σ,σ'*σ],H,J;rates=[κ,γ])

@test iszero(simplify_operators(he_diss.rhs[1] - ((-1.0im*ωc - 0.5κ)*a + (-1.0im*g)*σ)))
@test iszero(simplify_operators(he_diss.rhs[2] - ((-1.0im*g)*a + (-1.0im*ωa - 0.5γ)*σ + (2.0im*g)*a*σee)))
@test isequal(he_diss.rhs[3], simplify_operators((-γ)*σee + (1.0im*g)*a'*σ + (-1.0im*g)*a*σ'))

# Single-atom laser
@params ν
J = [a,σ,σ']
he_laser = heisenberg([a'*a,σ'*σ,a'*σ],H,J;rates=[κ,γ,ν])

@test isequal(he_laser.rhs[1], simplify_operators((-κ)*a'*a + (-1.0im*g)*a'*σ + (1.0im*g)*a*σ'))
@test isequal(he_laser.rhs[2], simplify_operators(ν + (-ν - γ)*σee + (1.0im*g)*a'*σ + (-1.0im*g)*a*σ'))
@test iszero(simplify_operators(he_laser.rhs[3] - ((1.0im*g)*σee + (-1.0im*g)*a'*a + (1.0im*ωc + -1.0im*ωa - 0.5*(κ + γ + ν))*a'*σ + (2.0im*g)*a'*a*σee)))

# Test sorting of longer term
@params λ Γ N
b = Destroy(h, :b)
yy = 1.0im*σ*b*(N-1)*λ*Γ
@test isequal(simplify_operators(yy), OperatorTerm(*, [-1.0im,Γ,λ,b,σ]) + OperatorTerm(*, [1.0im,N,Γ,λ,b,σ]))

end # testset
