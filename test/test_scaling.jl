using QuantumCumulants
using Test
using OrdinaryDiffEq
using ModelingToolkit

@testset "scaling" begin

order = 2
N_c = 1 #number of clusters
@cnumbers Δc κ Γ2 Γ3 Γ23 η ν3 ν2
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
@qnumbers a::Destroy(h,1)
S(i,j,c) = Transition(h,Symbol(:σ, c),i, j, 1+c) #c=cluster

# Hamiltonian
H = Δc*a'*a + sum(Δ2[c]*sum(S(2,2,c)) for c=1:N_c) + sum(Δ3[c]*sum(S(3,3,c)) for c=1:N_c) +
    sum(Ω3[c]*(sum(S(3,1,c)) + sum(S(1,3,c))) for c=1:N_c) + sum(g[c]*(a'*sum(S(1,2,c)) + a*sum(S(2,1,c))) for c=1:N_c)
# Collapse operators
J = [a,[S(1,2,c) for c=1:N_c]..., [S(1,3,c) for c=1:N_c]..., [S(2,3,c) for c=1:N_c]..., [S(3,3,c) for c=1:N_c]..., [S(2,2,c) for c=1:N_c]..., a'a]
rates = [κ,[Γ2 for i=1:N_c]...,[Γ3 for i=1:N_c]...,[Γ23 for i=1:N_c]...,[ν3 for i=1:N_c]...,[ν2 for i=1:N_c]...,η]

# Derive equation for average photon number
ops = [a'a, S(2,2,1)[1], a'*S(1,2,1)[1]]
he = meanfield(ops,H,J;rates=rates,order=2)

# Custom filter function -- include only phase-invariant terms
ϕ(x) = 0
ϕ(x::Destroy) = -1
ϕ(x::Create) = 1
function ϕ(t::Transition)
    if (t.i==1 && t.j==2) || (t.i==3 && t.j==2) || (t.i==:g && t.j==:e)
        -1
    elseif (t.i==2 && t.j==1) || (t.i==2 && t.j==3) || (t.i==:e && t.j==:g)
        1
    else
        0
    end
end
ϕ(avg::Average) = ϕ(avg.arguments[1])
function ϕ(t::QuantumCumulants.QMul)
    p = 0
    for arg in t.args_nc
        p += ϕ(arg)
    end
    return p
end
phase_invariant(x) = iszero(ϕ(x))

he_avg = cumulant_expansion(he,2)
he_scale = complete(he_avg;filter_func=phase_invariant, order=order, multithread=false)
@test length(he_scale) == 9

ps = [Δc; κ; Γ2; Γ3; Γ23; η; ν3; ν2; Δ2; Δ3; Ω3; g; N]
sys = ODESystem(he_scale)
u0 = zeros(ComplexF64, length(he_scale))

N0 = 1000/N_c
N_ = N0*[1.0 for c=1:N_c]

p0 = ps .=> [ones(length(ps)-1); N_]
prob1 = ODEProblem(sys,u0,(0.0,1.0),p0)
sol1 = solve(prob1, Tsit5(), reltol=1e-12, abstol=1e-12)

@test sol1.u[end][1] ≈ 0.0758608728203

avg = average(a*S(2,1,1)[2])
# TODO: @test sol1[avg], sol1, he_scale) == map(conj, get_solution(average(a'*S(1,2,1)[1]), sol1, he_scale)) == conj.(getindex.(sol1.u, 3))

## Two-level laser
M = 2
@cnumbers N
hf = FockSpace(:cavity)
ha = NLevelSpace(:atom, (:g,:e))
hc = ClusterSpace(ha, N, M)
h = tensor(hf, hc)

@qnumbers a::Destroy(h)
σ(i,j) = Transition(h,:σ,i,j)

@cnumbers Δ g κ γ ν

H = Δ*a'*a + g*sum(a'*σ(:g,:e)[i] + a*σ(:e,:g)[i] for i=1:M)
J = [a;[σ(:g,:e)[i] for i=1:M];[σ(:e,:g)[i] for i=1:M]]
rates = [κ; [γ for i=1:M]; [ν for i=1:M]]

he = meanfield(a'*a, H, J; rates=rates, order=2)

# Complete
he_scaled = complete(he;filter_func=phase_invariant)

names = he_scaled.names
avg = average(σ(:e,:g)[1]*σ(:e,:e)[2])
@test isequal(average(σ(:e,:e)[1]*σ(:e,:g)[2]), QuantumCumulants.substitute_redundants(avg,[QuantumCumulants.ClusterAon(2,1),QuantumCumulants.ClusterAon(2,2)],names))

@test QuantumCumulants.lt_reference_order(σ(:e,:g)[1],σ(:g,:e)[2])
@test !QuantumCumulants.lt_reference_order(σ(:g,:e)[1],σ(:e,:g)[2])

he_avg = cumulant_expansion(he_scaled,2)
@test isempty(find_missing(he_avg))

ps = (Δ, g, γ, κ, ν, N)
sys = ODESystem(he_avg)
p0 = ps .=> (0, 1.5, 0.25, 1, 4, 7)
u0 = zeros(ComplexF64, length(he_scaled))
prob = ODEProblem(sys, u0, (0.0, 50.0), p0)
sol = solve(prob, RK4(), abstol=1e-10, reltol=1e-10)

@test sol.u[end][1] ≈ 12.601868534

# Spectrum
corr = CorrelationFunction(a',a,he_avg;steady_state=true,filter_func=phase_invariant)
Spec = Spectrum(corr,ps)
s = Spec(range(-π, π, length=301), sol.u[end], getindex.(p0,2))
@test all(s .>= 0.0)

## Some abstract tests
M = 4
@cnumbers N
hc = FockSpace(:cavity)
hvib = FockSpace(:mode)
hcluster = ClusterSpace(hvib,N,M)
h = ⊗(hc, hcluster)
a = Destroy(h,:a,1)
b = Destroy(h,:b,2)

names = [:a,[Symbol(:b_,i) for i=1:M]]
scale_aons = [QuantumCumulants.ClusterAon(2,i) for i=1:M]

avg = average(a*b[1]*b[2]*b[2])
avg_sub = QuantumCumulants.substitute_redundants(avg, scale_aons, names)
@test isequal(average(a*b[1]*b[1]*b[2]), avg_sub)

avg = average(b[1]'*b[2]'*b[2])
avg_sub = QuantumCumulants.substitute_redundants(avg, scale_aons, names)
@test isequal(avg_sub, average(b[1]'*b[1]*b[2]'))

avg = average(b[1]*b[1]*b[2]'*b[2])
avg_sub = QuantumCumulants.substitute_redundants(avg, scale_aons, names)
@test isequal(avg_sub, average(b[1]'*b[1]*b[2]*b[2]))

avg = average(b[1]*b[2]'*b[3]*b[4]')
avg_sub = QuantumCumulants.substitute_redundants(avg, scale_aons, names)
@test isequal(avg_sub, average(b[1]'*b[2]'*b[3]*b[4]))

# Test Holstein
M = 2
hc = FockSpace(:cavity)
hvib = FockSpace(:mode)
hcluster = ClusterSpace(hvib,N,M)
h = ⊗(hc, hcluster)
@cnumbers G Δ κ γ Ω
a = Destroy(h,:a,1)
b = Destroy(h,:b,2)

H = Δ*a'*a + G*sum(b[i] + b[i]' for i=1:M)*a'*a + Ω*(a+a')
J = [a,b]
rates = [κ,γ]
ops = [a,a'*a,a*a,b[1],a*b[1],a'*b[1],b[1]'*b[1],b[1]*b[1],b[1]'*b[2],b[1]*b[2]]
he = meanfield(ops,H,J;rates=rates)

he_avg = cumulant_expansion(he,2)
@test isempty(find_missing(he_avg))

ps = (G,Δ,κ,γ,Ω,N)
sys = ODESystem(he_avg)

u0 = zeros(ComplexF64, length(he_avg))
p0 = ps .=> ones(length(ps))
prob = ODEProblem(sys, u0, (0.0,1.0), p0)
sol = solve(prob, Tsit5())

# Test molecule
M = 2 # Order
# Prameters
@cnumbers λ ν Γ η Δ γ N
# Hilbert space
h_in = NLevelSpace(:internal, 2)
hv = FockSpace(:vib)
hc = ClusterSpace(hv, N, M)
h = ⊗(h_in, hc)
# Operators
σ(i,j) = Transition(h, :σ, i, j)
b = Destroy(h, :b)
# Hamiltonian
H0 = Δ*σ(2,2) + ν*sum(b_'*b_ for b_ in b)
H_holstein = -1*λ*sum((b_' + b_ for b_ in b))*σ(2,2)
Hl = η*(σ(1,2) + σ(2,1))
H = H0 + H_holstein + Hl
# Jumps
J = [σ(1,2), b]
rates = [γ,Γ]
# Equations
ops = [σ(2,2),σ(1,2),b[1],σ(1,2)*b[1],σ(2,1)*b[1],σ(2,2)*b[1],b[1]'*b[1],b[1]*b[1],b[1]'*b[2],b[1]*b[2]]
he = meanfield(ops, H, J; rates=rates)
he_avg = cumulant_expansion(he,2)
@test isempty(find_missing(he_avg))
ps = (Δ,η,γ,λ,ν,Γ,N)
# Generate function
sys = ODESystem(he_avg)
p0 = ps .=> [ones(length(ps)-1); 4]
u0 = zeros(ComplexF64,length(he_avg))
prob1 = ODEProblem(sys,u0,(0.0,1.0),p0)
sol1 = solve(prob1,Tsit5(),abstol=1e-12,reltol=1e-12)
bdb1 = sol1[b[1]'*b[1]][end]
σ22_1 = sol1[σ(2,2)][end]
σ12_1 = sol1[σ(1,2)][end]


## Two clusters
N_c = 2
N = cnumbers([Symbol(:N_, i) for i=1:N_c]...)
M = 2
hf = FockSpace(:cavity)
ha = NLevelSpace(:atom, (:g,:e))
hc = [ClusterSpace(ha, N[i], M) for i=1:N_c]
h = tensor(hf, hc...)

@qnumbers a::Destroy(h)
σ(i,j,c) = Transition(h,Symbol(:σ_, c),i,j,c+1)

@cnumbers κ
ν = cnumbers([Symbol(:ν_, c) for c=1:N_c]...)
γ = cnumbers([Symbol(:γ_, c) for c=1:N_c]...)
Δ = cnumbers([Symbol(:Δ_, c) for c=1:N_c]...)
g = cnumbers([Symbol(:g_, c) for c=1:N_c]...)

H = sum(Δ[c]*σ(:e,:e,c)[k] for c=1:N_c, k=1:M) + sum(g[c]*(a'*σ(:g,:e,c)[i] + a*σ(:e,:g,c)[i]) for i=1:M, c=1:N_c)
J = [a;[σ(:g,:e,c) for c=1:N_c];[σ(:e,:g,c) for c=1:N_c]]
rates = [κ,γ...,ν...]

ops = [a'*a]
he = meanfield(ops, H, J; rates=rates, order=2)

# Scale
he_scaled = complete(he; filter_func=phase_invariant)

@test isempty(find_missing(he_scaled))

ps = (κ, Δ..., g..., γ..., ν..., N...)
sys = ODESystem(he_scaled)
if N_c==2
    p0 = ps .=> (1, [0 for i=1:N_c]..., [1.5 for i=1:N_c]..., [0.25 for i=1:N_c]..., [4 for i=1:N_c]..., 4, 3)
elseif N_c==3
    p0 = ps .=> (1, [0 for i=1:N_c]..., [1.5 for i=1:N_c]..., [0.25 for i=1:N_c]..., [4 for i=1:N_c]..., 2, 3, 2)
end
u0 = zeros(ComplexF64, length(he_scaled))
prob = ODEProblem(sys, u0, (0.0, 50.0), p0)
sol = solve(prob, RK4(), abstol=1e-10, reltol=1e-10)

@test sol.u[end][1] ≈ 12.601868534

### 6-level laser

order = 2
# Define parameters
@cnumbers Δc κ g Γ12 Γ13 Γ24 Γ34 Γ54 Γ16 Γ64 Δ3 Δ4 Δ5 Δ6 Ω13 Ω34 Ω54 Ω64 η ν12 ν13 ν34 ν54 ν64 N

# Define hilbert space
hf = FockSpace(:cavity)
ha = NLevelSpace(:atom,6)
ha_c = ClusterSpace(ha,N,order)
h = ⊗(hf, ha_c)
# Define the fundamental operators
a = Destroy(h,:a,1)
σ(i,j) = Transition(h,:σ,i,j)
@test length(σ(6,6)) == order

# Hamiltonian
H = Δc*a'a + Δ3*sum(σ(3,3)) + Δ4*sum(σ(4,4)) + Δ5*sum(σ(5,5)) + Δ6*sum(σ(6,6)) +
    Ω13*(sum(σ(3,1)) + sum(σ(1,3))) + Ω34*(sum(σ(3,4)) + sum(σ(4,3))) + Ω54*(sum(σ(5,4)) + sum(σ(4,5))) + Ω64*(sum(σ(6,4)) + sum(σ(4,6))) +
    g*(a'*sum(σ(1,2)) + a*sum(σ(2,1)))
# Collapse operators
J = [a, σ(1,2), σ(1,3), σ(2,4), σ(3,4), σ(5,4), σ(6,4), σ(1,6), σ(2,2), σ(3,3)+σ(4,4)+σ(5,5)+σ(6,6), σ(4,4)+σ(5,5)+σ(6,6), σ(5,5), σ(6,6), a'a]
rates = [κ, Γ12, Γ13, Γ24, Γ34, Γ54, Γ64, Γ16, ν12, ν13, ν34, ν54, ν64, η]

# Derive equation for average photon number
ops = [a'a, σ(2,2)[1], σ(3,3)[1], σ(4,4)[1], σ(5,5)[1], σ(6,6)[1]]
eqs_ops = meanfield(ops,H,J;rates=rates, order=order, multithread=true)

@test length(eqs_ops) == length(ops)


### 4th order 2-level laser
order = 4
# Define parameters
@cnumbers δA ΩA wA νA γ κ δc NA
# Define hilbert space
hf = FockSpace(:cavity)
haA = NLevelSpace(:atomA,2)
heA = ClusterSpace(haA, NA, order) #atom ensemble 1
h = ⊗(hf, heA)
a = Destroy(h,:a,1)
σA(i,j) = Transition(h,:σA,i,j,2)
# Hamiltonian
H = δc*a'a + δA*sum(σA(2,2)) + ΩA/2*(a'sum(σA(1,2)) + a*sum(σA(2,1)))
# dissipative processes
J = [a, σA(1,2), σA(2,1), σA(2,2)]
rates = [κ, γ, wA, γ, νA]
# Derive equation for average photon number
ops = [a'a, σA(2,2)[1]]
he_ops = meanfield(ops,H,J;rates=rates, multithread=true, order=order)

he_scale = complete(he_ops; filter_func=phase_invariant, order=order, multithread=true)
@test length(he_scale) == 18
@test isempty(find_missing(he_scale))

### 4th order synchronization

order = 4
# Define parameters
@cnumbers δA δB ΩA ΩB wA wB νA νB γ κ δc NA NB
# Define hilbert space
hf = FockSpace(:cavity)
haA = NLevelSpace(:atomA,2)
haB = NLevelSpace(:atomB,2)
heA = ClusterSpace(haA, NA, order) #atom ensemble 1
heB = ClusterSpace(haB, NB, order) #atom ensemble 2
h = ⊗(hf, heA, heB)
# Define the fundamental operators
a = Destroy(h,:a,1)
σA(i,j) = Transition(h,:σA,i,j,2)
σB(i,j) = Transition(h,:σB,i,j,3)
# Hamiltonian
H = δc*a'a + δA*sum(σA(2,2)) + δB*sum(σB(2,2)) +
    ΩA/2*(a'sum(σA(1,2)) + a*sum(σA(2,1))) + ΩB/2*(a'sum(σB(1,2)) + a*sum(σB(2,1)))
# dissipative processes
J = [a, σA(1,2), σA(2,1), σB(1,2), σB(2,1), σA(2,2), σB(2,2)]
rates = [κ, γ, wA, γ, wB, νA, νB]
# Derive equation for average photon number
ops = [a'a, σA(2,2)[1], σB(2,2)[1]]
he_ops = meanfield(ops,H,J;rates=rates, multithread=true, order=order)

he_scale = complete(he_ops; filter_func=phase_invariant, order=order, multithread=true)
@test length(he_scale) == 66
@test isempty(find_missing(he_scale))

sys = ODESystem(he_scale)
u0 = zeros(ComplexF64, length(he_scale))
ps = [δA, δB, ΩA, ΩB, wA, wB, νA, νB, γ, κ, δc, NA, NB]
p0 = ps.=>[1.0 + i/20 for i=1:length(ps)]
prob = ODEProblem(sys,u0,(0.0,1.0),p0)
sol = solve(prob,RK4())

uend = copy(sol.u[end])
@test length(unique(uend)) == 66
uend_f = filter(x->imag(x) != 0, uend)
@test 2*length(uend_f) == length(unique([uend_f; adjoint.(uend_f)]))

end # testset
