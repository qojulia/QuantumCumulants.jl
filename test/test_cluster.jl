using Qumulants

M = 2

# Define parameters
@parameters Δ g γ κ ν N

# Define hilbert space
hf = FockSpace(:cavity)
ha = NLevelSpace(Symbol(:atoms),(:g,:e))
hc = ClusterSpace(ha,N,M)
h = hf⊗hc

a = Destroy(h,:a)
S(i,j) = Transition(h,:σ,i, j)

H = Δ*a'*a + g*(a'*sum(S(:g,:e)) + a*sum(S(:e,:g)))

# Collapse operators
J = [a;S(:g,:e);S(:e,:g)]
rates = [κ;[γ for i=1:M];[ν for i=1:M]]

# Derive equation for average photon number
he_n = heisenberg(a'*a,H,J;rates=rates)
he_avg = average(he_n, M)

phase_invariant(x) = iszero(phase(x))
phase(avg::Average) = phase(avg.operator)
phase(op::Destroy) = -1
phase(op::Create) = 1
function phase(t::Transition)
    lvls = Qumulants.levels(t.hilbert, t.aon)
    i = findfirst(isequal(t.i), lvls)
    j = findfirst(isequal(t.j), lvls)
    if i < j
        -1
    elseif i > j
        1
    else
        0
    end
end
phase(op::OperatorTerm{<:typeof(*)}) = sum(phase(arg) for arg in op.arguments)

he_nophase = complete(he_avg;order=M,filter_func=phase_invariant)

ps = (Δ, g, γ, κ, ν, N)
meta_f = build_ode(he_nophase, ps)
f = Meta.eval(meta_f)

using OrdinaryDiffEq
u0 = zeros(ComplexF64, length(he_nophase))
p0 = (0, 1.5, 0.25, 1, 4, 2)
prob = ODEProblem(f, u0, (0.0, 50.0), p0)
sol = solve(prob, RK4())


# Spectrum
corr = CorrelationFunction(a', a, he_nophase; filter_func=phase_invariant, steady_state=true)
s = Spectrum(corr, ps)

using PyPlot; pygui(true)
w_ls = range(-2pi,2pi,length=501)
plot(w_ls, s(w_ls, sol.u[end], p0))


### Test harmonic oscillator cluster
using Qumulants

M = 2 # Order

# Prameters
@parameters λ ν Γ η Δ γ N

# Hilbert space
h_in = NLevelSpace(:internal, (:g,:e))
hv = ClusterSpace(FockSpace(:vib), N, M)
h = h_in ⊗ hv

# Operators
b = Destroy(h, :b)
σ(i,j) = Transition(h, :σ, i, j)


# Hamiltonian
H0 = Δ*σ(:e,:e) + N*λ^2*ν*σ(:e,:e) + ν*(b'*b)
H_holstein = -1*λ*ν*(sum(b') + sum(b))*σ(:e,:e)
Hl = η*(σ(:g,:e) + σ(:e,:g))
H = H0 + H_holstein + Hl

H_ = simplify_operators(H)

# Jumps
J = [σ(:g,:e); b]
rates = [γ;[Γ for i=1:M]]

he_in = average(heisenberg(σ(:e,:e), H, J; rates=rates), M)
he = complete(he_in; order=M)

ps = (λ, ν, Γ, η, Δ, γ, N)
meta_f = build_ode(he, ps)
f = Meta.eval(meta_f)

u0 = zeros(ComplexF64, length(he))
p0 = [ones(length(ps)-1); 3] # TODO
prob = ODEProblem(f, u0, (0.0,20.0), p0)
sol = solve(prob,RK4())

dt = 1e-3
tspan = (0.0,dt)
prob = ODEProblem(f, u0, tspan, p0)
sol = solve(prob,RK4(),adaptive=false,dt=dt)
