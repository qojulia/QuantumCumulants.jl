using QuantumCumulants
using SymbolicUtils
using Test

@testset "parameters" begin

hf = FockSpace(:cavity)
ha = NLevelSpace(:atom,(:g,:e))
h = hf⊗ha

a = Destroy(h,:a)

@cnumbers p q
# @test p*q isa NumberTerm
@test p*q*a isa QuantumCumulants.QTerm

σ = Transition(h,:σ,:g,:e)
σee = Transition(h,:σ,:e,:e)

# JC
@cnumbers ωc ωa g
H = ωc*a'*a + ωa*σ'*σ + g*(a'*σ + σ'*a)

da = commutator(1im*H,a)
@test iszero(simplify(da - ((-1im*g)*σ + (-1im*ωc)*a)))
ds = commutator(1im*H,σ)
@test iszero(simplify(ds - ((-1im*g)*a + (-1im*ωa)*σ + (2im*g)*a*σee)))
dn = commutator(1im*H,a'*a)
@test iszero(simplify(dn - ((-1im*g)*a'*σ + (1im*g)*a*σ')))

he = heisenberg([a,σ,a'*a],H)
@test iszero(simplify(he.operator_equations[1].rhs - da))
@test iszero(simplify(he.operator_equations[2].rhs - ds))
@test iszero(simplify(he.operator_equations[3].rhs - dn))

# Lossy JC
@cnumbers κ γ
J = [a,σ]
he_diss = heisenberg([a,σ,σ'*σ],H,J;rates=[κ,γ])

@test iszero(simplify(he_diss.operator_equations[1].rhs - ((-1im*ωc - 0.5κ)*a + (-1im*g)*σ); polynorm=true))
@test iszero(simplify(he_diss.operator_equations[2].rhs - ((-1im*g)*a + (-1im*ωa - 0.5γ)*σ + (2im*g)*a*σee); polynorm=true))
@test iszero(simplify(he_diss.operator_equations[3].rhs - ((-γ)*σee + (1im*g)*a'*σ + (-1im*g)*a*σ')))

# Single-atom laser
@cnumbers ν
J = [a,σ,σ']
he_laser = heisenberg([a'*a,σ'*σ,a'*σ],H,J;rates=[κ,γ,ν])

@test iszero(simplify(he_laser.operator_equations[1].rhs - ((-κ)*a'*a + (-1im*g)*a'*σ + (1im*g)*a*σ')))
@test iszero(simplify(he_laser.operator_equations[2].rhs - (ν + (-ν - γ)*σee + (1im*g)*a'*σ + (-1im*g)*a*σ'); polynorm=true))
@test iszero(simplify(he_laser.operator_equations[3].rhs - ((1im*g)*σee + (-1im*g)*a'*a + (1im*ωc + -1im*ωa - 0.5*(κ + γ + ν))*a'*σ + (2im*g)*a'*a*σee); polynorm=true))

# Test sorting of longer term
@cnumbers λ Γ N
b = Destroy(h, :b)
yy = 1im*σ*b*(N-1)*λ*Γ
using SymbolicUtils: Term
@test isequal(simplify(yy), QuantumCumulants.QMul(im*Γ*λ*(N-1), [b, σ]))

# Test numbers in commutator
@cnumbers Δ
h = NLevelSpace(:a, 2)
s(i,j) = Transition(h, :s, i, j)
H = Δ*s(2,2) - Δ*s(1,1)
H = simplify(H)
@test iszero(simplify(heisenberg(s(1,2), H).operator_equations[1].rhs  + 2im*Δ*s(1,2)))

end # testset
