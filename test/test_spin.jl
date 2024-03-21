using QuantumCumulants
using ModelingToolkit
using OrdinaryDiffEq
using Test

@testset "spin" begin

hs1 = SpinSpace(:Spin1)
hs2 = SpinSpace(:Spin2)
h = hs1 ⊗ hs2

s(axis) = Sigma(hs1, :σ, axis) # axis ∈ [1,2,3] → [x,y,z]
@test isequal(s(1)*s(2),1im*s(3))
@test !isequal(s(1)*s(2),1im*s(2))
@test isequal(s(1)*s(3),-1im*s(2))
@test isequal(s(3)*s(3),1)
@test isequal(s(3)*s(1),1im*s(2))
@test isequal(s(1)*s(2)*s(3),1im)

σ(i, axis) = Sigma(h,Symbol(:σ_,i), axis, i)
σx(i) = σ(i, 1)
σy(i) = σ(i, 2)
σz(i) = σ(i, 3)

sx(i) = σ(i, :x)
sy(i) = σ(i, :y)
sz(i) = σ(i, :z)
sx(1) == σx(1)

@test isequal(σx(2)*σx(1), σx(1)*σx(2))
@test isequal(σy(2)*σx(1), σx(1)*σy(2))
@test isequal(σz(2)*σz(1), σz(1)*σz(2))

@cnumbers J
Δi(i) = cnumber(Symbol(:Δ_,i))
H = Δi(1)*σz(1) + Δi(2)*σz(2) + J*σx(1)*σx(2)

ops = [σz(1)]
eqs = meanfield(ops, H)
eqs_c = complete(eqs, order=2)
@test length(eqs_c) == 6

ops2 = [σx(1), σy(1), σz(1), σx(2), σy(2), σz(2)]
eqs2 = meanfield(ops2, H)
eqs2_c = complete(eqs2, order=2)
@test length(eqs2_c) == 14

# Time evolution driven Dicke model
@cnumbers Δ g κ η
hf = FockSpace(:cavity)
ha1 = NLevelSpace(:atom1,2)
ha2 = NLevelSpace(:atom2,2)
h = hf ⊗ ha1 ⊗ ha2
a = Destroy(h,:a)
s1(i,j) = Transition(h,:s1,i,j,2)
s2(i,j) = Transition(h,:s2,i,j,3)
H = Δ*a'*a + g*(a' + a)*(s1(2,1) + s1(1,2) + s2(2,1) + s2(1,2)) + η*(a' + a)
J = [a]
rates = [κ]
#
eq = meanfield(a'a,H,J;rates=rates,order=2)
eqs = complete(eq)
ps = [Δ, g, κ, η]
@named sys = ODESystem(eqs)
u0 = zeros(ComplexF64, length(eqs))
p0 = [0.5, 1.0, 1.25, 0.85]
prob = ODEProblem(sys,u0,(0.0,0.5),ps.=>p0)
sol = solve(prob,RK4())

s1y = get_solution(sol, 1im*(s1(1,2) - s1(2,1)))[end]
s2x = get_solution(sol, (s2(2,1) + s2(1,2)))[end]
s1z = get_solution(sol, (s1(2,2) - s1(1,1)))[end]
n = sol[a'a][end]
s1zs2z = get_solution(sol, (s1(2,2) - s1(1,1))*(s2(2,2) - s2(1,1)))[end]
s1xs2y = get_solution(sol, (s1(2,1) + s1(1,2))*1im*(s2(1,2) - s2(2,1)))[end]

### Spin description
hs1 = SpinSpace(:spin1)
hs2 = SpinSpace(:spin2)
h_ = hf ⊗ hs1 ⊗ hs2
a2 = Destroy(h_,:a2)
σ1(axis) = Sigma(h_,:σ1,axis,2)
σ2(axis) = Sigma(h_,:σ2,axis,3)
H = Δ*a2'*a2 + g*(a2' + a2)*(σ1(1) + σ2(1)) + η*(a2' + a2)
J = [a2]
rates = [κ]
#
eq_ = meanfield([σ1(3), σ2(3), σ1(3)*σ2(3)],H,J;rates=rates,order=2)
eqs_ = complete(eq_)
eqs_.states
ps = [Δ, g, κ, η]
@named sys_ = ODESystem(eqs_)
u0_ = zeros(ComplexF64, length(eqs_))
u0_[1] = u0_[2] = -1
u0_[3] = 1
p0 = [0.5, 1.0, 1.25, 0.85]
prob_ = ODEProblem(sys_,u0_,(0.0,0.5),ps.=>p0)
sol_ = solve(prob_,RK4())

s1y_ = get_solution(sol_, σ1(2))[end]
s2x_ = sol_[σ2(1)][end]
s1z_ = get_solution(sol_, σ1(3))[end]
n_ = sol_[a2'a2][end]
s1zs2z_ = sol_[σ1(3)*σ2(3)][end]
s1xs2y_ = sol_[σ1(1)*σ2(2)][end]

@test isapprox(s1y, s1y_; atol=1e-5)
@test isapprox(s2x, s2x_; atol=1e-5)
@test isapprox(s1z, s1z_; atol=1e-5)
@test isapprox(s1zs2z, s1zs2z_; atol=1e-5)
@test isapprox(n, n_; atol=1e-5)
@test isapprox(s1zs2z, s1zs2z_; atol=1e-5)
@test isapprox(s1xs2y, s1xs2y_; atol=1e-5)

end # testset
