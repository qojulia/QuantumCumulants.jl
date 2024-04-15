using QuantumCumulants
using ModelingToolkit
using OrdinaryDiffEq
using QuantumOpticsBase
using Test

@testset "spin" begin

hs1 = PauliSpace(:Spin1)
hs2 = PauliSpace(:Spin2)
h = hs1 ⊗ hs2

s(axis) = Pauli(hs1, :σ, axis) # axis ∈ [1,2,3] → [x,y,z]
@test isequal(s(1)*s(2),1im*s(3))
@test !isequal(s(1)*s(2),1im*s(2))
@test isequal(s(1)*s(3),-1im*s(2))
@test isequal(s(3)*s(3),1)
@test isequal(s(3)*s(1),1im*s(2))
@test isequal(s(1)*s(2)*s(3),1im)

σ(i, axis) = Pauli(h,Symbol(:σ_,i), axis, i)
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
@cnumbers Δ_ g κ η
hf = FockSpace(:cavity)
ha1 = NLevelSpace(:atom1,2)
ha2 = NLevelSpace(:atom2,2)
h = hf ⊗ ha1 ⊗ ha2
a = Destroy(h,:a)
s1(i,j) = Transition(h,:s1,i,j,2)
s2(i,j) = Transition(h,:s2,i,j,3)
H = Δ_*a'*a + g*(a' + a)*(s1(2,1) + s1(1,2) + s2(2,1) + s2(1,2)) + η*(a' + a)
J = [a]
rates = [κ]
#
eq = meanfield(a'a,H,J;rates=rates,order=2)
eqs = complete(eq)
ps = [Δ_, g, κ, η]
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
hs1 = PauliSpace(:spin1)
hs2 = PauliSpace(:spin2)
h_ = hf ⊗ hs1 ⊗ hs2
a2 = Destroy(h_,:a2)
σ1(axis) = Pauli(h_,:σ1,axis,2)
σ2(axis) = Pauli(h_,:σ2,axis,3)
H = Δ_*a2'*a2 + g*(a2' + a2)*(σ1(1) + σ2(1)) + η*(a2' + a2)
J = [a2]
rates = [κ]
#
eq_ = meanfield([σ1(3), σ2(3), σ1(3)*σ2(3)],H,J;rates=rates,order=2)
eqs_ = complete(eq_)
eqs_.states
ps = [Δ_, g, κ, η]
@named sys_ = ODESystem(eqs_)
u0_ = zeros(ComplexF64, length(eqs_))
u0_[1] = u0_[2] = -1
u0_[3] = 1

# initial state: numeric conversion
b_field = FockBasis(4)
bp1 = SpinBasis(1/2) 
bp2 = SpinBasis(1/2)
b = b_field ⊗ bp1 ⊗ bp2 
ψf = fockstate(b_field,0)
ψ1 = spindown(bp1)
ψ2 = spindown(bp2)
ψ_ = ψf ⊗ ψ1 ⊗ ψ2
u0_pauli = initial_values(eqs_, ψ_)
@test all(u0_test .≈ u0_cs1)


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

### collective spin ###
hcs1 = SpinSpace(:Spin1)
hcs2 = SpinSpace(:Spin2)
h = hcs1 ⊗ hcs2

S(axis) = Spin(hcs1, :S, axis) # axis ∈ [1,2,3] → [x,y,z]

S(1)==S(:x)
S(2)==S(:Y)
S(:z)==S(:Z)
S(1)≠S(2)

isequal(simplify(S(:x)*S(:y) - S(:y)*S(:x)), 1im*S(:z))
isequal(simplify(S(:x)*S(:z) - S(:z)*S(:x)), -1im*S(:y))
isequal(simplify(S(:y)*S(:z) - S(:z)*S(:y)), 1im*S(:1))

isequal(S(:x)*S(:y), 1*S(:x)*S(:y))
!isequal(S(:x)*S(:y), S(:x)*S(:z))
isequal(S(1)*S(1), (S(1))^2)
isequal(simplify(S(1) + S(2)), simplify(S(2) + S(1)))

# error
isequal(2*S(1)*S(2)*S(3), S(1)*S(2)*S(3)*2)

S(i, axis) = Spin(h,Symbol(:S_,i), axis, i)
Sx(i) = S(i, 1)
Sy(i) = S(i, 2)
Sz(i) = S(i, 3)
Sm(i) = Sx(i) - 1im*Sy(i)
Sp(i) = Sx(i) + 1im*Sy(i)

isequal(Sx(2)*Sx(1), Sx(1)*Sx(2))
isequal(Sy(2)*Sx(1), Sx(1)*Sy(2))
isequal(Sz(2)*Sz(1), Sz(1)*Sz(2))
isequal(simplify(Sy(1)Sz(2)Sx(1)), Sx(1)Sy(1)Sz(2) - 1im*Sz(1)*Sz(2))


### simple time evolution 
@cnumbers δcs Ωcs gcs Γcs
Hcs1 = δcs/2*Sz(1)
Jcs1 = [Sm(1)]
Rcs1 = [Γcs]

ops_cs1 = [Sx(1), Sy(1), Sz(1)]
eqs_cs1 = meanfield(ops_cs1,Hcs1,Jcs1;rates=Rcs1,order=2)
eqs_cs1_c = complete(eqs_cs1)
@named sys_cs1 = ODESystem(eqs_cs1_c);

eqs_cs1_c.states
u0_cs1 = zeros(ComplexF64, length(eqs_cs1_c))
# full exciation 
Ncs1 = 20
Ncs1_ = Ncs1/2
Ncs2 = 8
Ncs2_ = Ncs2/2

u0_cs1[3] = Ncs1_ # z
u0_cs1[6] = u0_cs1[7] = Ncs1/4 # xx, yy
u0_cs1[8] = 1im*Ncs1/4 # xy
u0_cs1[9] = Ncs1_*Ncs1_ # zz

# initial state: numeric conversion
bs1 = SpinBasis(Ncs1_) 
bs2 = SpinBasis(Ncs2_) # random
b = bs1 ⊗ bs2 # used space is prodcut space, but second space is not used
ψ1 = spinup(bs1)
ψ2 = spinup(bs2)
ψ = ψ1 ⊗ ψ2
u0_test = initial_values(eqs_cs1_c, ψ)
@test all(u0_test .≈ u0_cs1)

ps_cs1 = [δcs, Γcs]
p0_cs1 = [0, 1]
prob_cs1 = ODEProblem(sys_cs1,u0_cs1,(0.0, 0.2), ps_cs1.=>p0_cs1)
sol_cs1 = solve(prob_cs1,Tsit5(),abstol=1e-8,reltol=1e-8)

@test sol_cs1[Sz(1)][1] == 10.0
@test real.(sol_cs1[Sz(1)][end]) < 0
@test abs(imag.(sol_cs1[Sz(1)][end])) < 0.01
