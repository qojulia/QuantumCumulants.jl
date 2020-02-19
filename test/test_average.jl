using Qumulants
using SymPy
using Test

@testset "average" begin

# Test single mode
a = Destroy(:a)
a_avg = average(a)

a2 = a*a
a2_avg = average(a2)
@test isa(a2_avg,SymPy.Sym)
a2_mf = average(a2,1)
@test a2_mf == a_avg*a_avg

ad = a'
ad_avg = average(ad)
ad_avg2 = replace_adjoints(ad_avg, [a],1)[1]
@test ad_avg2 == a_avg'

@test average(one(a)) == 1

# Test nlevel averaging
σ(i,j) = Transition(:σ,i,j,2;GS=2)
@test average(σ(1,2)) == replace_adjoints(average(σ(2,1)), [σ(1,2)], 1)[1]'
@test average(σ(1,1)).__pyobject__.is_real
@test average(σ(2,2)) == average(one(σ(2,2)) - σ(1,1)) # Why does this fail when using `1 - average(σ(1,1))` ?
tmp = average(one(σ(2,2)) - σ(1,1))
tmp2 = 1-average(σ(1,1))

# Test single-atom laser
σ(i,j) = Transition(:σ,i,j,2;GS=1)
a = Destroy(:a) ⊗ Identity()
s = Identity() ⊗ σ(1,2)
sz = simplify_operators(2*s'*s - one(s))
sps = simplify_operators(s'*s)
ωc, ωa, g, γ, κ, ν = (0.1, 3, 0.5, 0.25, 1.0, 4)
H = ωc*a'*a + ωa*s'*s + g*(a'*s + a*s')
J = [sqrt(κ)*a,sqrt(γ)*s,sqrt(ν)*s']

dada = heisenberg(a'*a, H, J)
dsps = heisenberg(s'*s, H, J)
dads = heisenberg(a'*s, H, J)

@test average(dada.rhs) == -κ*average(a'*a) - 1.0im*g*average(a'*s) + 1.0im*g*average(a*s')
@test average(dsps.rhs) == -(γ + ν)*average(sps) + ν + 1.0im*g*average(a'*s) - 1.0im*g*average(a*s')
@test average(dads.rhs,2) == -(1.0im*(ωa - ωc) + 0.5κ + 0.5ν + 0.5γ)*average(a'*s) + 1.0im*g*average(sps) + 2.0im*g*average(sps*a'*a,2) - 1.0im*g*average(a'*a)

end # testset
