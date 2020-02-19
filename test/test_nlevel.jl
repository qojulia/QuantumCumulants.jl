using Qumulants
using Test

@testset "nlevel" begin

# Test single two-level atom
σ(i,j) = Transition(:σ,i,j,1:2;GS=2)
σ12 = σ(1,2)
tmp = σ12*σ12
@test simplify_operators(tmp) == zero(σ12)
@test simplify_operators(σ12*σ12') == σ(1,1)
@test σ12*σ12' == σ(1,1)

σ(i,j) = Transition(:σ,i,j,(:g,:e))
σge = σ(:g,:e)
tmp = σge*σge
@test simplify_operators(tmp) == zero(σge)
@test simplify_operators(σge'*σge) == σ(:e,:e)

# Rabi oscillation model
σ(i,j) = Transition(:σ,i,j,1:2;GS=2)
Δ,Ω,γ = (0.33,0.2,1.0)
H = Δ*σ(1,1) + Ω*(σ(1,2) + σ(2,1))
J = [sqrt(γ)*σ(2,1)]

σ21 = σ(2,1)
dσ21 = simplify_operators(1.0im*(H*σ(2,1) - σ(2,1)*H))
dσ21_sim = Qumulants.apply_comms(dσ21)
sz = simplify_operators(σ(1,1) - σ(2,2))
@test sz == 2*σ(1,1) - one(σ(1,1))

@test dσ21 == dσ21_sim == heisenberg(σ21,H).rhs == simplify_operators(-1.0im*Δ*σ21 + 1.0im*Ω*sz)

dσ11 = heisenberg(σ(1,1),H)
@test dσ11.rhs == simplify_operators(1.0im*Ω*(-σ21' + σ21))

dσ21_qle = simplify_operators(1.0im*(H*σ(2,1) - σ(2,1)*H) + sum(j'*σ21*j - 0.5*(j'*j*σ21 + σ21*j'*j) for j=J))
dσ21_qle_sim = Qumulants.apply_comms(dσ21_qle)
@test dσ21_qle_sim == heisenberg(σ21,H,J).rhs == simplify_operators(-1.0im*Δ*σ21 + 1.0im*Ω*sz - 0.5γ*σ21)
dσ11_qle = heisenberg(σ(1,1),H,J)
@test dσ11_qle.rhs == simplify_operators(1.0im*Ω*(-σ21' + σ21) - γ*σ(1,1))

# Test single-atom laser
a = Destroy(:a) ⊗ Identity()
s = Identity() ⊗ σ(2,1)
sz = simplify_operators(2*s'*s - one(s))
ωc, ωa, g, γ, κ, ν = (0.1, 3, 0.5, 0.25, 1.0, 4)
H = ωc*a'*a + ωa*s'*s + g*(a'*s + a*s')
J = [sqrt(κ)*a,sqrt(γ)*s,sqrt(ν)*s']

dada = heisenberg(a'*a, H, J)
dsps = heisenberg(s'*s, H, J)
dads = heisenberg(a'*s, H, J)

@test dada.rhs == simplify_operators(-κ*a'*a -1.0im*g*(a'*s - a*s'))
@test dsps.rhs == simplify_operators(-(γ + ν)*s'*s + ν*one(s) + 1.0im*g*(a'*s - a*s'))
@test dads.rhs == simplify_operators(1.0im*(ωc - ωa)*a'*s - 0.5*(κ + ν + γ)*a'*s + 1.0im*g*sz*a'*a + 1.0im*g*s'*s)

end # testset
