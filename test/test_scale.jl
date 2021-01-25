using Qumulants
using OrdinaryDiffEq
using Test

@testset "scaling" begin

order = 2
N_c = 1 #number of clusters
@parameters Δc κ Γ2 Γ3 Γ23 η ν3 ν2
Δ2 = [Parameter(Symbol(:Δ2_, i)) for i=1:N_c]
Δ3 = [Parameter(Symbol(:Δ3_, i)) for i=1:N_c]
Ω3 = [Parameter(Symbol(:Ω3_, i)) for i=1:N_c]
g = [Parameter(Symbol(:g_, i)) for i=1:N_c]
N = [Parameter(Symbol(:N_, i)) for i=1:N_c]

# Define hilbert space
hf = FockSpace(:cavity)
ha = [NLevelSpace(Symbol(:atoms, j),3) for j=1:N_c]
ha_c = [ClusterSpace(ha[j],N[j],order) for j=1:N_c]
h = ⊗(hf, ha_c...)
# Define the fundamental operators
a = Destroy(h,:a,1)
S(i,j,c) = Transition(h,Symbol(:σ, c),i, j, 1+c) #c=cluster

# Hamiltonian
H = Δc*a'*a + sum(Δ2[c]*sum(S(2,2,c)) for c=1:N_c) + sum(Δ3[c]*sum(S(3,3,c)) for c=1:N_c) +
    sum(Ω3[c]*(sum(S(3,1,c)) + sum(S(1,3,c))) for c=1:N_c) + sum(g[c]*(a'*sum(S(1,2,c)) + a*sum(S(2,1,c))) for c=1:N_c)
H = simplify_operators(H)
# Collapse operators
J = [a,[S(1,2,c) for c=1:N_c]..., [S(1,3,c) for c=1:N_c]..., [S(2,3,c) for c=1:N_c]..., [S(3,3,c) for c=1:N_c]..., [S(2,2,c) for c=1:N_c]..., a'a]
rates = [κ,[Γ2 for i=1:N_c]...,[Γ3 for i=1:N_c]...,[Γ23 for i=1:N_c]...,[ν3 for i=1:N_c]...,[ν2 for i=1:N_c]...,η]

# Derive equation for average photon number
ops = [a'a, S(2,2,1)[1]]
he_ops = heisenberg(ops,H,J;rates=rates)
# Custom filter function -- include only phase-invaraint terms
ϕ(x) = 0
ϕ(x::Destroy) = -1
ϕ(x::Create) = 1
function ϕ(t::Transition)
    if (t.i==1 && t.j==2) || (t.i==3 && t.j==2)
        -1
    elseif (t.i==2 && t.j==1) || (t.i==2 && t.j==3)
        1
    else
        0
    end
end
ϕ(avg::Average) = ϕ(avg.operator)
function ϕ(t::OperatorTerm)
    @assert t.f === (*)
    p = 0
    for arg in t.arguments
        p += ϕ(arg)
    end
    return p
end
phase_invariant(x) = iszero(ϕ(x))

he_scale = complete(he_ops;filter_func=phase_invariant, order=order, multithread=true)
@test length(he_scale) == 9

ps = [Δc; κ; Γ2; Γ3; Γ23; η; ν3; ν2; Δ2; Δ3; Ω3; g; N]
meta_f = build_ode(he_scale, ps)
f = Meta.eval(meta_f)
u0 = zeros(ComplexF64, length(he_scale))

N0 = 1000/N_c
N_ = N0*[1.0 for c=1:N_c]

p0 = [ones(length(ps)-1); N_]
prob1 = ODEProblem(f,u0,(0.0,1.0),p0)
sol1 = solve(prob1, Tsit5(), reltol=1e-12, abstol=1e-12)
sol1.u[end]
# spectrum
corr = CorrelationFunction(a', a, he_scale; filter_func=phase_invariant, steady_state=true)
@test length(corr.de) == 3
meta_f_corr = build_ode(corr, ps)
f_corr = Meta.eval(meta_f_corr)
s = Spectrum(corr,ps)

##################
### 2 clusters ###
##################

order = 2
N_c = 2 #number of clusters
@parameters Δc κ Γ2 Γ3 Γ23 η ν3 ν2
Δ2 = [Parameter(Symbol(:Δ2_, i)) for i=1:N_c]
Δ3 = [Parameter(Symbol(:Δ3_, i)) for i=1:N_c]
Ω3 = [Parameter(Symbol(:Ω3_, i)) for i=1:N_c]
g = [Parameter(Symbol(:g_, i)) for i=1:N_c]
N = [Parameter(Symbol(:N_, i)) for i=1:N_c]

# Define hilbert space
hf = FockSpace(:cavity)
ha = [NLevelSpace(Symbol(:atoms, j),3) for j=1:N_c]
ha_c = [ClusterSpace(ha[j],N[j],order) for j=1:N_c]
h = ⊗(hf, ha_c...)
# Define the fundamental operators
a = Destroy(h,:a,1)
S(i,j,c) = Transition(h,Symbol(:σ, c),i, j, 1+c) #c=cluster

# Hamiltonian
H = Δc*a'*a + sum(Δ2[c]*sum(S(2,2,c)) for c=1:N_c) + sum(Δ3[c]*sum(S(3,3,c)) for c=1:N_c) +
    sum(Ω3[c]*(sum(S(3,1,c)) + sum(S(1,3,c))) for c=1:N_c) + sum(g[c]*(a'*sum(S(1,2,c)) + a*sum(S(2,1,c))) for c=1:N_c)
H = simplify_operators(H)
# Collapse operators
J = [a,[S(1,2,c) for c=1:N_c]..., [S(1,3,c) for c=1:N_c]..., [S(2,3,c) for c=1:N_c]..., [S(3,3,c) for c=1:N_c]..., [S(2,2,c) for c=1:N_c]..., a'a]
rates = [κ,[Γ2 for i=1:N_c]...,[Γ3 for i=1:N_c]...,[Γ23 for i=1:N_c]...,[ν3 for i=1:N_c]...,[ν2 for i=1:N_c]...,η]

# Derive equation for average photon number
ops = [a'a, S(2,2,1)[1]]
he_ops = heisenberg(ops,H,J;rates=rates)
# Custom filter function -- include only phase-invaraint terms

he_scale = complete(he_ops;filter_func=phase_invariant, order=order, multithread=true)
@test length(he_scale) == 21

ps = [Δc; κ; Γ2; Γ3; Γ23; η; ν3; ν2; Δ2; Δ3; Ω3; g; N]
meta_f = build_ode(he_scale, ps)
f = Meta.eval(meta_f)
u0 = zeros(ComplexF64, length(he_scale))

N0 = 1000/N_c
N_ = N0*[1.0 for c=1:N_c]
ps
p0 = [ones(length(ps)-N_c); N_]
prob2 = ODEProblem(f,u0,(0.0,1.0),p0)
sol2 = solve(prob2, Tsit5(), reltol=1e-12, abstol=1e-12)

# spectrum
corr = CorrelationFunction(a', a, he_scale; filter_func=phase_invariant, steady_state=true)
@test length(corr.de) == 5
meta_f_corr = build_ode(corr, ps)
f_corr = Meta.eval(meta_f_corr)
s = Spectrum(corr,ps)

@test isapprox(abs.(sol2.u[end][1]), abs.(sol1.u[end][1]), rtol=1e-6)
@test isapprox(abs.(sol2.u[end][2]), abs.(sol1.u[end][2]), rtol=1e-6)

################
### Holstein ###
################

M = 3 # Order
# Prameters
@parameters λ ν Γ η Δ γ N
# Hilbert space
h_in = NLevelSpace(:internal, 2)
hv = FockSpace(:vib)
hv_c = ClusterSpace(hv,N,M)
h = ⊗(h_in, hv_c)
# Operators
σ(i,j) = Transition(h, :σ, i, j)
b = Destroy(h, :b)
# Hamiltonian
H0 = Δ*σ(2,2) + ν*sum(b'b)
H_holstein = -1*λ*sum((b' + b))*σ(2,2)
Hl = η*(σ(1,2) + σ(2,1))
H = H0 + H_holstein + Hl
# Jumps
J = [σ(1,2), b]
rates = [γ,Γ]

# Equations
ops = b[1]
he0 = heisenberg(ops, H, J; rates=rates)
he = complete(he0, order=M, multithread=true)
@test length(he) == 27
ps = (Δ,η,γ,λ,ν,Γ,N)
# Generate function
f = generate_ode(he,ps;check_bounds=true)
p0 = [ones(length(ps)-1)..., 4]

u0 = zeros(ComplexF64,length(he))
prob1 = ODEProblem(f,u0,(0.0,1.0),p0)
sol1 = solve(prob1,Tsit5(),abstol=1e-12,reltol=1e-12)

bdb1 = get_solution(b[1]'b[1], sol1, he)[end]
σ22_1 = get_solution(σ(2,2), sol1, he)[end]
σ12_1 = get_solution(σ(1,2), sol1, he)[end]


########################
### explicit 4 atoms ###
########################

n = 4 # Number of vibration modes
h_in = NLevelSpace(:internal,2)
hv = [FockSpace(Symbol(:vib, i)) for i=1:n]
h = ⊗(h_in, hv...)
# Operators
σ(i,j) = Transition(h, :σ, i, j)
bb(k) = Destroy(h, Symbol(:b_, k), k+1)
# Hamiltonian
H0 = Δ*σ(2,2) + ν*sum(bb(k)'*bb(k) for k=1:n)
H_holstein = - λ*sum((bb(k)' + bb(k)) for k=1:n)*σ(2,2)
Hl = η*(σ(1,2) + σ(2,1))
H = H0 + H_holstein + Hl
# Jumps
J = [σ(1,2); [bb(k) for k=1:n]]
rates = [γ;[Γ for i=1:n]]
# Equations
he_in = average(heisenberg(σ(2,2), H, J; rates=rates), M)
he = complete(he_in; order=M, multithread=true);
@test length(he) == 154
ps = (Δ,η,γ,λ,ν,Γ)
# Generate function
f = generate_ode(he,ps;check_bounds=true)
# Numerical parameters
p0 = ones(length(ps))
u0 = zeros(ComplexF64,length(he))
prob2 = ODEProblem(f,u0,(0.0,1.0),p0)
sol2 = solve(prob2,Tsit5(),abstol=1e-12,reltol=1e-12);

bdb2 = get_solution(bb(1)'bb(1), sol2, he)[end]
σ22_2 = get_solution(σ(2,2), sol2, he)[end]
σ12_2 = get_solution(σ(1,2), sol2, he)[end]

@test isapprox(abs.(bdb2), abs.(bdb1), rtol=1e-6)
@test isapprox(abs.(σ22_2), abs.(σ22_1), rtol=1e-6)
@test isapprox(abs.(σ12_2), abs.(σ12_1), rtol=1e-6)

end #testset
