# # Ramsey Spectroscopy

# A typical problem in quantum optics is the interrogation of an atom with an external driving field. In this brief example we apply [Ramsey interferometry](https://en.wikipedia.org/wiki/Ramsey_interferometry) on a single two-level atom. The Hamiltonian is

# $H = - \Delta \sigma^{22} +  \Omega(t) (\sigma^{21} + \sigma^{12}),$

# with $\Delta = \omega_l - \omega_a$ the laser-atom detuning, and $\Omega(t)$ a time-dependent driving field. Additionally we include atomic decay and dephasing with rates $\Gamma$ and $\nu$.

using QuantumCumulants
using OrdinaryDiffEq, ModelingToolkitBase
using Plots

@variables Δ Ω Γ ν
@register_symbolic f(t)

h = NLevelSpace(:atom, 2)
σ(i, j) = Transition(h, :σ, i, j)

# Seed a meanfield call to grab an MTK-aware iv `t`, then build the
# time-dependent Hamiltonian against it.
eqs_seed = meanfield([σ(2, 2), σ(1, 2)], -Δ * σ(2, 2), [σ(1, 2), σ(2, 2)]; rates = [Γ, ν])
t = eqs_seed.iv
H = -Δ * σ(2, 2) + f(t) * Ω / 2 * (σ(1, 2) + σ(2, 1))

eqs = meanfield([σ(2, 2), σ(1, 2)], H, [σ(1, 2), σ(2, 2)]; rates = [Γ, ν], iv = t)
complete!(eqs)

sys = to_system(eqs; name = :sys)
sys_c = mtkcompile(sys)

Γ_ = 1.0
Ω_ = 500Γ_
Δ_ = 0Γ_
ν_ = 0.2Γ_

tp = π / (2Ω_)
tf = 1 / (20Γ_)

function f(t)
    if t < tp || (t > tp + tf && t < 2tp + tf)
        return 1.0
    else
        return 0.0
    end
end

ps = Dict(Γ => Γ_, Ω => Ω_, Δ => Δ_, ν => ν_)
u0 = initial_values(eqs)

prob = ODEProblem(sys_c, merge(u0, ps), (0.0, 2tp + tf))
sol = solve(prob, Tsit5(), maxiters = 1.0e7)

ts = sol.t
s22 = real.(get_solution(sol, σ(2, 2), eqs).(ts))
plot(ts, s22, xlabel = "tΓ", ylabel = "⟨σ22⟩", legend = false, size = (600, 300))

# Scanning over the detuning produces the well-known Ramsey fringes.

Δ_ls = [-2000:20:2000;] .* Γ_
s22_ls = zeros(length(Δ_ls))

for i in eachindex(Δ_ls)
    ps_i = Dict(Γ => Γ_, Ω => Ω_, Δ => Δ_ls[i], ν => ν_)
    prob_i = ODEProblem(sys_c, merge(u0, ps_i), (0.0, 2tp + tf))
    sol_i = solve(prob_i, Tsit5(); reltol = 1.0e-8, abstol = 1.0e-8)
    s22_ls[i] = real(get_solution(sol_i, σ(2, 2), eqs)(sol_i.t[end]))
end

plot(Δ_ls, s22_ls, xlabel = "Δ/Γ", ylabel = "⟨σ22⟩", legend = false, size = (600, 300))
