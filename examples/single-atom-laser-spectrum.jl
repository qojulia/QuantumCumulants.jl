# # Single-atom Laser (meanfield setup)

# A single-atom laser system: cavity mode plus a two-level atom with incoherent
# pumping. The Hamiltonian is
#
# $H = \Delta a^\dagger a + g(a^\dagger \sigma^{ge} + a\sigma^{eg}).$
#
# We construct the equations of motion at second order. The original example
# also computed the laser spectrum via `CorrelationFunction` + `Spectrum`; in v1.0
# the spectrum machinery is being re-implemented and only the meanfield/MTK
# boundary is exercised here.

using QuantumCumulants
using OrdinaryDiffEq, ModelingToolkitBase
using Plots

@variables Δ g γ κ ν

hf = FockSpace(:cavity)
ha = NLevelSpace(:atom, (:g, :e))
h = hf ⊗ ha

a = Destroy(h, :a)
s = Transition(h, :σ, :g, :e)

H = Δ * a' * a + g * (a' * s + a * s')
J = [a, s, s']
rates = [κ, γ, ν]

eqs = meanfield(a' * a, H, J; rates = rates, order = 2)
complete!(eqs)

sys = to_system(eqs; name = :laser)
sys_c = mtkcompile(sys)

u0 = initial_values(eqs)
p0 = Dict(Δ => 1.0, g => 1.5, γ => 0.25, κ => 1.0, ν => 4.0)
prob = ODEProblem(sys_c, merge(u0, p0), (0.0, 10.0))
sol = solve(prob, Tsit5())

ts = sol.t
n_t = real.(get_solution(sol, a' * a, eqs).(ts))
plot(ts, n_t; xlabel = "t", ylabel = "⟨a⁺a⟩")
