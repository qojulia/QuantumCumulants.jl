# # Mollow Triplet (meanfield setup)

# Resonance spectrum of a single, coherently driven atom (Mollow 1969). The
# Hamiltonian is
# $H = -\Delta\sigma^{ee} + \Omega(\sigma^{ge} + \sigma^{eg})$
# with collapse operator $\sigma^{ge}$ at rate $\gamma$. The original example
# computed the laser spectrum via `Spectrum(c, ps)` (Laplace-domain solve); in
# v1.0 the spectrum machinery is being re-implemented and only the steady-state
# integration is shown here.

using QuantumCumulants
using ModelingToolkitBase, OrdinaryDiffEq
using Plots

h = NLevelSpace(:atom, (:g, :e))

@variables Δ Ω γ
σ(i, j) = Transition(h, :σ, i, j)
H = Δ * σ(:e, :e) + Ω * (σ(:g, :e) + σ(:e, :g))
J = [σ(:g, :e)]

eqs = meanfield([σ(:e, :g), σ(:e, :e)], H, J; rates = [γ])
complete!(eqs)

sys = to_system(eqs; name = :sys)
sys_c = mtkcompile(sys)

p0 = Dict(Δ => 0.0, Ω => 2.0, γ => 1.0)
u0 = initial_values(eqs)
prob = ODEProblem(sys_c, merge(u0, p0), (0.0, 20.0))
sol = solve(prob, Tsit5())

ts = sol.t
σee = real.(get_solution(sol, σ(:e, :e), eqs).(ts))
plot(ts, σee; xlabel = "γt", ylabel = "Excited-state population")
