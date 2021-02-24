using Qumulants
using Test
using OrdinaryDiffEq

M = 2
hf = FockSpace(:cavity)
ha = [NLevelSpace(Symbol(:atom, i), (:g,:e)) for i=1:M]
h = tensor(hf, ha...)

@qnumbers a::Destroy(h)
σ(i,j,k) = Transition(h,Symbol(:σ_,k),i,j,k+1)

@cnumbers Δ g κ γ ν

H = Δ*a'*a + g*sum(a'*σ(:g,:e,i) + a*σ(:e,:g,i) for i=1:M)
J = [a;[σ(:g,:e,i) for i=1:M];[σ(:e,:g,i) for i=1:M]]
rates = [κ; [γ for i=1:M]; [ν for i=1:M]]

ops = [a'*a, a'*σ(:g,:e,1), σ(:e,:e,1), σ(:e,:g,1)*σ(:g,:e,2)]

he = heisenberg(ops, H, J; rates=rates)

ϕ(x) = 0
ϕ(::Destroy) = -1
ϕ(::Create) = 1
function ϕ(t::Transition)
    if (t.i==:e && t.j==:g)
        1
    elseif (t.i==:g && t.j==:e)
        -1
    else
        0
    end
end
ϕ(avg::Average) = ϕ(avg.arguments[1])
function ϕ(t::QTerm)
    @assert t.f === (*)
    p = 0
    for arg in t.arguments
        p += ϕ(arg)
    end
    return p
end
phase_invariant(x) = iszero(ϕ(x))

# Scale
@cnumbers N
he_scaled = scale(he, [2,3], N)

avg = average(σ(:e,:g,1)*σ(:e,:e,2))
@test isequal(average(σ(:e,:e,1)*σ(:e,:g,2)), Qumulants.substitute_redundants(avg,[2,3],[:a,:σ_1,:σ_2]))

@test Qumulants.lt_reference_order(σ(:e,:g,1),σ(:g,:e,2))
@test !Qumulants.lt_reference_order(σ(:g,:e,1),σ(:e,:g,2))

he_avg = average(he_scaled,2)
missed = find_missing(he_avg)
filter!(!phase_invariant, missed)
@test isequal(missed, find_missing(he_avg))

subs = Dict(missed .=> 0)
he_nophase = substitute(he_avg, subs)
@test isempty(find_missing(he_nophase))

ps = (Δ, g, γ, κ, ν, N)
f = generate_ode(he_nophase, ps)
p0 = (0, 1.5, 0.25, 1, 4, 7)
u0 = zeros(ComplexF64, length(he_scaled))
prob = ODEProblem(f, u0, (0.0, 50.0), p0)
sol = solve(prob, RK4(), abstol=1e-10, reltol=1e-10)

@test sol.u[end][1] ≈ 12.601868534

# Some abstract tests
M = 4
hc = FockSpace(:cavity)
hvib = [FockSpace(Symbol(:mode, i)) for i=1:M]
h = ⊗(hc, hvib...)
a = Destroy(h,:a,1)
b = [Destroy(h,Symbol(:b_,i),i+1) for i=1:M]

names = [:a;[Symbol(:b_,i) for i=1:M]]

avg = average(a*b[1]*b[2]*b[2])
Qumulants.substitute_redundants(avg, [2,3], names)

avg = average(b[1]'*b[2]'*b[2])
Qumulants.substitute_redundants(avg, [2,3], names)

avg = average(b[1]*b[1]*b[2]'*b[2])
Qumulants.substitute_redundants(avg, [2,3], names)

avg = average(b[1]*b[2]'*b[3]*b[4]')
Qumulants.substitute_redundants(avg, [2:5;], names)

# Test Holstein
M = 2
hc = FockSpace(:cavity)
hvib = [FockSpace(Symbol(:mode, i)) for i=1:M]
h = ⊗(hc, hvib...)
@cnumbers G Δ κ γ Ω
a = Destroy(h,:a,1)
b = [Destroy(h,Symbol(:b_,i),i+1) for i=1:M]

avg = average(a*b[1]*b[2]*b[2])
Qumulants.substitute_redundants(avg, [2,3], names)

avg = average(b[1]'*b[2]'*b[2])
Qumulants.substitute_redundants(avg, [2,3], names)

avg = average(b[1]*b[1]*b[2]'*b[2])
Qumulants.substitute_redundants(avg, [2,3], names)

H = Δ*a'*a + G*sum(b[i] + b[i]' for i=1:M)*a'*a + Ω*(a+a')
J = [a;b]
rates = [κ;[γ for i=1:M]]
ops = [a,a'*a,a*a,b[1],a*b[1],a'*b[1],b[1]'*b[1],b[1]*b[1],b[1]'*b[2],b[1]*b[2]]
he = heisenberg(ops,H,J;rates=rates)

@cnumbers N
he_scaled = scale(he,[2,3],N;simplify=false)

he_avg = average(he_scaled,2)
@test isempty(find_missing(he_avg))

ps = (G,Δ,κ,γ,Ω,N)
f = generate_ode(he_avg, ps)

u0 = zeros(ComplexF64, length(he_avg))
p0 = ones(length(ps))
prob = ODEProblem(f, u0, (0.0,1.0), p0)
sol = solve(prob, Tsit5())

# Test molecule
M = 2 # Order
# Prameters
@cnumbers λ ν Γ η Δ γ N
# Hilbert space
h_in = NLevelSpace(:internal, 2)
hv = [FockSpace(Symbol(:vib, i)) for i=1:M]
h = ⊗(h_in, hv...)
# Operators
σ(i,j) = Transition(h, :σ, i, j)
b = [Destroy(h, Symbol(:b_, i), i+1) for i=1:M]
# Hamiltonian
H0 = Δ*σ(2,2) + ν*sum(b_'*b_ for b_ in b)
H_holstein = -1*λ*sum((b_' + b_ for b_ in b))*σ(2,2)
Hl = η*(σ(1,2) + σ(2,1))
H = H0 + H_holstein + Hl
# Jumps
J = [σ(1,2); b]
rates = [γ;[Γ for i=1:M]]
# Equations
ops = [σ(2,2),σ(1,2),b[1],σ(1,2)*b[1],σ(2,1)*b[1],σ(2,2)*b[1],b[1]'*b[1],b[1]*b[1],b[1]'*b[2],b[1]*b[2]]
he = heisenberg(ops, H, J; rates=rates)
he_scaled = scale(he,[2,3],N)
he_avg = average(he_scaled,2)
@test isempty(find_missing(he_avg))
ps = (Δ,η,γ,λ,ν,Γ,N)
# Generate function
f = generate_ode(he_avg,ps;check_bounds=true)
p0 = [ones(length(ps)-1); 4]
u0 = zeros(ComplexF64,length(he_avg))
prob1 = ODEProblem(f,u0,(0.0,1.0),p0)
sol1 = solve(prob1,Tsit5(),abstol=1e-12,reltol=1e-12)
bdb1 = get_solution(b[1]'b[1], sol1, he)[end]
σ22_1 = get_solution(σ(2,2), sol1, he)[end]
σ12_1 = get_solution(σ(1,2), sol1, he)[end]


## Two clusters
N_c = 2
M = 2
hf = FockSpace(:cavity)
ha = [NLevelSpace(Symbol(:atom, i, j), (:g,:e)) for j=1:N_c, i=1:M]
h = tensor(hf, ha...)

@qnumbers a::Destroy(h)
σ(i,j,c) = [Transition(h,Symbol(:σ_,c, :_, k, :_),i,j,k+1+M*(c-1)) for k=1:M]

@cnumbers κ
ν = cnumbers([Symbol(:ν_, c) for c=1:N_c]...)
γ = cnumbers([Symbol(:γ_, c) for c=1:N_c]...)
Δ = cnumbers([Symbol(:Δ_, c) for c=1:N_c]...)
g = cnumbers([Symbol(:g_, c) for c=1:N_c]...)

H = sum(Δ[c]*σ(:e,:e,c)[k] for c=1:N_c, k=1:M) + sum(g[c]*(a'*σ(:g,:e,c)[i] + a*σ(:e,:g,c)[i]) for i=1:M, c=1:N_c)
H = qsimplify(H)
J = QNumber[a]
rates = [κ]
for c=1:N_c
    append!(J, σ(:g,:e,c))
    append!(rates, [γ[c] for i=1:M])
end
for c=1:N_c
    append!(J, σ(:e,:g,c))
    append!(rates, [ν[c] for i=1:M])
end

ops = [a'*a; [a'*σ(:g,:e,c)[1] for c=1:N_c]; [σ(:e,:e,c)[1] for c=1:N_c]; [σ(:e,:g,c)[1]*σ(:g,:e,c)[2] for c=1:N_c]]
for i=1:N_c
    for j = i+1:N_c
        push!(ops, σ(:e,:g,i)[1]*σ(:g,:e,j)[1])
    end
end

he = heisenberg(ops, H, J; rates=rates)

ϕ(x) = 0
ϕ(::Destroy) = -1
ϕ(::Create) = 1
function ϕ(t::Transition)
    if (t.i==:e && t.j==:g)
        1
    elseif (t.i==:g && t.j==:e)
        -1
    else
        0
    end
end
ϕ(avg::Average) = ϕ(avg.arguments[1])
function ϕ(t::QTerm)
    @assert t.f === (*)
    p = 0
    for arg in t.arguments
        p += ϕ(arg)
    end
    return p
end
phase_invariant(x) = iszero(ϕ(x))

# Scale
N = [cnumbers([Symbol(:N_, i) for i=1:N_c]...)...]
scale_aons = [[2+k for k=(i-1)*M:M*i-1] for i=1:N_c]
he_scaled = scale(he, scale_aons, N)

he_avg = average(he_scaled,2)
missed = find_missing(he_avg)
filter!(!phase_invariant, missed)
@test isequal(missed, find_missing(he_avg))

subs = Dict(missed .=> 0)
he_nophase = substitute(he_avg, subs)
@test isempty(find_missing(he_nophase))

ps = (κ, Δ..., g..., γ..., ν..., N...)
f = generate_ode(he_nophase, ps)
if N_c==2
    p0 = (1, [0 for i=1:N_c]..., [1.5 for i=1:N_c]..., [0.25 for i=1:N_c]..., [4 for i=1:N_c]..., 4, 3)
elseif N_c==3
    p0 = (1, [0 for i=1:N_c]..., [1.5 for i=1:N_c]..., [0.25 for i=1:N_c]..., [4 for i=1:N_c]..., 2, 3, 2)
end
u0 = zeros(ComplexF64, length(he_scaled))
prob = ODEProblem(f, u0, (0.0, 50.0), p0)
sol = solve(prob, RK4(), abstol=1e-10, reltol=1e-10)

@test sol.u[end][1] ≈ 12.601868534

# order = 2
# N_c = 2 #number of clusters
# @cnumbers Δc κ Γ2 Γ3 Γ23 η ν3 ν2
# Δ2 = cnumbers([Symbol(:Δ2_, i) for i=1:N_c]...)
# Δ3 = cnumbers([Symbol(:Δ3_, i) for i=1:N_c]...)
# Ω3 = cnumbers([Symbol(:Ω3_, i) for i=1:N_c]...)
# g =cnumbers( [Symbol(:g_, i) for i=1:N_c]...)
# N =cnumbers( [Symbol(:N_, i) for i=1:N_c]...)
# # Define hilbert space
# hf = FockSpace(:cavity)
# ha = [NLevelSpace(Symbol(:atoms, j, k),3) for k=1:order, j=1:N_c]
# h = ⊗(hf, ha...)
# # Define the fundamental operators
# a = Destroy(h,:a,1)
# S(i,j,c) = [Transition(h,Symbol(:σ_, c, :_, k, :_),i, j, 1+k+(order)*(c-1)) for k=1:order] #c=cluster
# # Hamiltonian
# H = Δc*a'*a + sum(Δ2[c]*sum(S(2,2,c)) for c=1:N_c) + sum(Δ3[c]*sum(S(3,3,c)) for c=1:N_c) +
#     sum(Ω3[c]*(sum(S(3,1,c)) + sum(S(1,3,c))) for c=1:N_c) + sum(g[c]*(a'*sum(S(1,2,c)) + a*sum(S(2,1,c))) for c=1:N_c)
# H = qsimplify(H)
# # Collapse operators
# J = [a,[S(1,2,c)[k] for c=1:N_c, k=1:order]..., [S(1,3,c)[k] for c=1:N_c, k=1:order]..., [S(2,3,c)[k] for c=1:N_c, k=1:order]..., [S(3,3,c)[k] for c=1:N_c, k=1:order]..., [S(2,2,c)[k] for c=1:N_c, k=1:order]..., a'a]
# rates = [κ,[Γ2 for i=1:order*N_c]...,[Γ3 for i=1:order*N_c]...,[Γ23 for i=1:order*N_c]...,[ν3 for i=1:order*N_c]...,[ν2 for i=1:order*N_c]...,η]
# # Derive equation for average photon number
# ops = [a'*a; [S(2,2,c)[1] for c=1:N_c]; [a'*S(1,2,c)[1] for c=1:N_c]; [S(2,1,c)[1]*S(1,2,c)[2] for c=1:N_c];
#         [S(2,2,c)[1]*S(2,2,c)[2] for c=1:N_c]]
# for i=1:N_c
#     for j=i+1:N_c
#         push!(ops, S(2,1,i)[1]*S(1,2,j)[1])
#         push!(ops, S(2,2,i)[1]*S(2,2,j)[1])
#     end
# end
# he_ops = heisenberg(ops,H,J;rates=rates)
# scale_aons = [[2,3],[4,5]]
# N = [N...]
# he_scale = scale(he_ops, scale_aons, N)
# he_avg = average(he_scale, 2)
# # Custom filter function -- include only phase-invaraint terms
# # he_scale = complete(he_ops;filter_func=phase_invariant, order=order, multithread=true)
#
# ϕ(x) = 0
# ϕ(x::Destroy) = -1
# ϕ(x::Create) = 1
# function ϕ(t::Transition)
#     if (t.i==1 && t.j==2) || (t.i==3 && t.j==2)
#         -1
#     elseif (t.i==2 && t.j==1) || (t.i==2 && t.j==3)
#         1
#     else
#         0
#     end
# end
# ϕ(avg::Average) = ϕ(avg.arguments[1])
# function ϕ(t::QTerm)
#     @assert t.f === (*)
#     p = 0
#     for arg in t.arguments
#         p += ϕ(arg)
#     end
#     return p
# end
# phase_invariant(x) = iszero(ϕ(x))
#
# missed = find_missing(he_avg)
# subs = Dict(filter(!phase_invariant, missed) .=> 0)
# he_avg = qsimplify(substitute(he_avg, subs))
# missed = find_missing(he_avg)
# filter!(phase_invariant, missed)
# for i=1:length(missed)
#     missed[i] = Qumulants.substitute_redundants(missed[i], scale_aons, names)
# end
# missed = unique_ops(missed)
#
# @test length(he_scale) == 21
# ps = [Δc; κ; Γ2; Γ3; Γ23; η; ν3; ν2; Δ2; Δ3; Ω3; g; N]
# meta_f = build_ode(he_scale, ps)
# f = Meta.eval(meta_f)
# u0 = zeros(ComplexF64, length(he_scale))
# N0 = 1000/N_c
# N_ = N0*[1.0 for c=1:N_c]
# ps
# p0 = [ones(length(ps)-N_c); N_]
# prob2 = ODEProblem(f,u0,(0.0,1.0),p0)
# sol2 = solve(prob2, Tsit5(), reltol=1e-12, abstol=1e-12)
#
# names = Qumulants.get_names(he_ops)
#
# avg = average(a'*S(1,2,2)[2])
# Qumulants.is_redundant_aon(avg, [4,5])
# Qumulants.substitute_redundants(avg, [4,5], names)
