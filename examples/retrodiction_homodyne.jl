# # Prediction and Retrodiction

# In this example, we implement the prediction and retrodiction of a continuously monitored Gaussian system as described in [J. Zhang and K. Mølmer, Phys. Rev. A 96, 062131 (2017)](10.1103/PhysRevA.96.062131). The system is a harmonic oscillator with frequency $\Omega$, where the decay $\sqrt{\Gamma} a$ is continuously monitored with an efficiency $\eta$.  
# We start by loading the needed packages and specifying the model.

using QuantumCumulants
using ModelingToolkit
using OrdinaryDiffEq
using StochasticDiffEq
using Plots
using Random # hide

## Hilbert space, operators and parameter
h = PhaseSpace(:motion)
@qnumbers x::Position(h, 1)
p = Momentum(h, :p, 1)

# Since `SymbolicUtils` has problems with fractions we avoid terms like $1/\sqrt{2}$ and $1/m$ by defining $s = 1/\sqrt{2}$ and $m = 1$.

@rnumbers Ω Γ η s
m = 1

@syms t::Real
@register_symbolic Ydot(t) # measurement current
@register_symbolic Ydot_rev(t) # time-reverse measurement current

a = (x + 1im*p)*s
nothing # hide

## Hamiltonian
H = p^2/(2m) + 0.5m*Ω^2*x^2

## Lindblad
J = [a]
R = [Γ]
eff = [η]

# We derive the stochastic differential equation for the measurement backaction and solve the dynamics. 

ops = [x, p, x*x, x*p, p*p]
eqs = meanfield(ops, H, J; rates = R, efficiencies = eff, order = 2)

# We want to note here that the time evolution of the covariance matrix would in principle be deterministic (no noise term). However, this is not the case for the set of equations we derive her, which makes them not optimal to describe the system. Nevertheless, the numerics is stable and the results are still correct. 

## numerical parameter
Ω_ = 1
Γ_ = 1/6
η_ = 0.5

Tend = 3.0/Γ_
dt = Tend/5e4 # time step
T_saveat = [0:dt:Tend;]

ps = [Ω, Γ, η, s]
pn = [Ω_, Γ_, η_, 1/√2]

## initial state
x0 = 5.0;
p0 = 0.0
u0 = [x0, p0, 5+x0*x0, im/2+x0*p0, 5+p0*p0] # σx = σpp = 5
nothing # hide

## define numerical SDE
Random.seed!(2) # hide
@named sys_fw = SDESystem(eqs)
dict_fw = merge(Dict(unknowns(sys_fw) .=> u0), Dict(ps .=> pn))
noise = StochasticDiffEq.RealWienerProcess(0.0, 0.0; save_everystep = true)
prob_fw = SDEProblem(sys_fw, dict_fw, (0, Tend); noise = noise)
nothing # hide

#

sol_fw = solve(prob_fw, EM(); dt = dt, saveat = T_saveat)
nothing # hide

pp = plot(sol_fw.t * Γ_, real(sol_fw[x]), xlabel = "Γt", ylabel = "⟨x⟩", legend = false)

# The above trajectory is now used to extract a measurement record, which serves as our "experimental" result. In the following, we calculate the "experimental" measurement current from the noisy trajectory to simulate the system again with this measurement (Kalman filter).

t_W = noise.t # noise time steps
W_W = noise.W # noise values

## calculate measurement record 
N = length(t_W)-1
dY_W = zeros(N)
dYdt_W = zeros(N)

sol_x = real(sol_fw[x])
for k = 1:N
    dW = W_W[k+1] - W_W[k]
    # cd_p_c = (sol_x[k] + sol_x[k+1])/2*√(2)*√(Γ_)
    cd_p_c = sol_x[k]*√(2)*√(Γ_)

    dYdt_W[k] = dW/dt + √η_*cd_p_c
    dY_W[k] = dW + √η_*cd_p_c*dt
end
dYdt_data = [dYdt_W; dYdt_W[end]]
dYdt(t) = dYdt_data[Int(round(t/dt))+1] # measurement current

# We modify the noiseless equations to add the (deterministic) term from the measurement current. 

## deterministic Kalman filtering equations
eqs_det = meanfield(ops, H, J; rates = R, order = 2)
function f_measure(lhs, rhs)
    term_ = √(η*Γ)*(lhs*a + a'*lhs - average(a+a')*lhs)*(Ydot(t) - √(η*Γ)*average(a+a'))
    term = cumulant_expansion(average(term_), 2)
    return rhs + term
end
eqs_kal = modify_equations(eqs_det, f_measure)

#

Ydot(t) = dYdt(t)
@named sys_kal = System(eqs_kal)
prob_kal = ODEProblem(sys_kal, dict_fw, (0, Tend))
nothing # hide

#

sol_kal = solve(prob_kal, Euler(); dt = dt, saveat = T_saveat)

nothing # hide

#

x_fw = real(sol_kal[x])
xx_fw = real(sol_kal[x*x])
σx_fw = [xx_fw[k] - x_fw[k]*x_fw[k] for k = 1:(N+1)]
sqr_σx_fw = sqrt.(σx_fw)

p_fw = real(sol_kal[p])
pp_fw = real(sol_kal[p*p])
σp_fw = [pp_fw[k] - p_fw[k]*p_fw[k] for k = 1:(N+1)]
sqr_σp_fw = sqrt.(σp_fw)

#

p1 = plot(
    t_W*Γ_,
    x_fw;
    size = (500, 250),
    label = "",
    xlabel = "Γt",
    ylabel = "⟨x⟩",
    ylim = (-10, 10),
    grid = true,
)
plot!(t_W .* Γ_, x_fw; ribbon = sqr_σx_fw, fillalpha = 0.2, label = "")

# We can see that the same dynamics as before is obtained.

# The same measurement current is now used to propagate backward in time to obtain the retrodicted, more precise, state of the system. This procedure can also be described with the so-called past quantum state, and is similar to the smoothing procedure in classical measurement theory, where past and future observations are used to estimate the state of a system. To this end, we drive the backward propagation equation, where we again modify the usual equations. We need to adjust the recycling term of the Lindblad operator and include the measurement backaction. 

eqs_back_kal = meanfield(ops, -H, J; rates = R, order = 2)
eqs_back_kal_c = complete(eqs_back_kal)

## adjust recycling term for the back propagation
function f_back_lind(lhs, rhs)
    term_1 = Γ*a*lhs*a' - Γ*a'lhs*a # adapt recycling term
    term_2 = -Γ*(average(a*a') - average(a'a))*lhs
    term = cumulant_expansion(average(term_1 + term_2), 2)
    return rhs + term
end
eqs_back_kal_c1 = modify_equations(eqs_back_kal_c, f_back_lind)

## include measurement backaction 
function f_back_kal(lhs, rhs)
    term_ = √(η*Γ)*(lhs*a' + a*lhs - average(a+a')*lhs)*(Ydot_rev(t) - √(η*Γ)*average(a+a'))
    term = cumulant_expansion(average(term_), 2)
    return rhs + term
end
eqs_back_kal_c2 = modify_equations(eqs_back_kal_c1, f_back_kal)

#

# deterministic smoothing equations
Ydot_rev(t) = dYdt(Tend-t)
@named sys_back_kal = System(eqs_back_kal_c2)
u0_bw = [0.0+0im, 0, 100, 0, 100]
dict_bw = merge(Dict(unknowns(sys_back_kal) .=> u0_bw), Dict(ps .=> pn))
prob_bw_kal = ODEProblem(sys_back_kal, dict_bw, (0, Tend))
nothing # hide

#

sol_bw = solve(prob_bw_kal, Euler(); dt = dt, saveat = T_saveat)
nothing # hide

#

x_bw = reverse(real(sol_bw[x]))
xx_bw = reverse(real(sol_bw[x*x]))
σx_bw = [xx_bw[k] - x_bw[k]*x_bw[k] for k = 1:(N+1)]
sqr_σx_bw = sqrt.(σx_bw)

p_bw = reverse(real(sol_bw[p]))
pp_bw = reverse(real(sol_bw[p*p]))
σp_bw = [pp_bw[k] - p_bw[k]*p_bw[k] for k = 1:(N+1)]
sqr_σp_bw = sqrt.(σp_bw)
nothing # hide

#

p2 = plot(
    t_W*Γ_,
    x_bw;
    size = (500, 250),
    label = "",
    xlabel = "Γt",
    ylabel = "⟨x⟩",
    ylim = (-10, 10),
    grid = true,
)
plot!(t_W*Γ_, x_bw; ribbon = sqr_σx_bw, fillalpha = 0.2, label = "")


# For a Gaussian system we can combine the mean values and covariances of the forward and backward equations in the following way to determine the retrodicted expressions: # TODO

x_pq = (x_fw .* σx_bw + x_bw .* σx_fw) ./ (σx_fw .+ σx_bw)
σx_pq = [1/(1/σx_fw[i] + 1/σx_bw[i]) for i = 1:length(σx_bw)]
sqr_σx_pq = sqrt.(σx_pq)
nothing # hide

#

p3 = plot(
    t_W*Γ_,
    x_pq;
    size = (500, 250),
    label = "",
    xlabel = "Γt",
    ylabel = "⟨x⟩",
    ylim = (-10, 10),
    grid = true,
)
plot!(t_W*Γ_, x_pq; ribbon = sqr_σx_pq, fillalpha = 0.2, label = "")

# Note that the uncertainty of the retrodicted position is reduced. 

# ## Package versions

# These results were obtained using the following versions:

using InteractiveUtils
versioninfo()

using Pkg
Pkg.status(
    ["QuantumCumulants", "OrdinaryDiffEq", "ModelingToolkit", "StochasticDiffEq", "Plots"],
    mode = PKGMODE_MANIFEST,
)
