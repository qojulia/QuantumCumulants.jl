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

@test round.(abs.(sol2.u[end][1]); digits=6) == round.(abs.(sol1.u[end][1]); digits=6)
@test round.(abs.(sol2.u[end][2]); digits=6) == round.(abs.(sol1.u[end][2]); digits=6)


################
### Holstein ###
################



end #testset
