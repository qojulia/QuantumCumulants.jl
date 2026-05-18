# # Laser with Filter Cavities (meanfield setup)

# Single-atom laser coupled to a frequency-comb of filter cavities. The full
# example also computed the cavity-photon-number spectrum via `Spectrum`; in
# v1.0 the spectrum machinery is still being ported and only the meanfield
# setup is exercised here.

using QuantumCumulants
using OrdinaryDiffEq, ModelingToolkitBase
using Plots

M = 3 # number of filter cavities

@variables N Δ g γ κ ν gf δ κf

hc = FockSpace(:cavity)
hf = FockSpace(:filter)
ha = NLevelSpace(:atom, (:g, :e))
h = hc ⊗ hf ⊗ ha

@qnumbers a::Destroy(h, 1)
b(k) = IndexedOperator(Destroy(h, :b, 2), Index(h, Symbol(:i, k), M, hf))
σ(α, β) = Transition(h, :σ, α, β)

H = Δ * a' * a + g * (a' * σ(:g, :e) + a * σ(:e, :g))

J = [a, σ(:g, :e), σ(:e, :g)]
rates = [κ, γ, ν]

eqs = meanfield(a' * a, H, J; rates = rates, order = 2)
complete!(eqs)

sys = to_system(eqs; name = :filter_cavity)
sys_c = mtkcompile(sys)

u0 = initial_values(eqs)
p0 = Dict(Δ => 1.0, g => 1.5, γ => 0.25, κ => 1.0, ν => 4.0)
prob = ODEProblem(sys_c, merge(u0, p0), (0.0, 10.0))
sol = solve(prob, Tsit5())

ts = sol.t
n_t = real.(get_solution(sol, a' * a, eqs).(ts))
plot(ts, n_t; xlabel = "t", ylabel = "⟨a⁺a⟩")
