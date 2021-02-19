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
sol = solve(prob, RK4())

@test sol.u[end][1] ≈ 12.6018748

# Test Holstein
M = 2
hc = FockSpace(:cavity)
hvib = [FockSpace(Symbol(:mode, i)) for i=1:M]
h = ⊗(hc, hvib...)
@cnumbers G Δ κ γ Ω
a = Destroy(h,:a,1)
b = [Destroy(h,Symbol(:b_,i),i+1) for i=1:M]

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
