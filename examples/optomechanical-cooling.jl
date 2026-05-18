# # Optomechanical Cooling

# In this example, we implement a cooling scheme based on radiation pressure coupling of light to a mechanical oscillator inside an optical cavity. The model is based on [C. Genes et. al., Phys. Rev. A 77, 033804 (2008)](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.77.033804), with Hamiltonian

# $H = -\hbar\Delta a^\dagger a + \hbar\omega_m b^\dagger b + \hbar Ga^\dagger a \left(b + b^\dagger\right) + \hbar E \left(a + a^\dagger\right),$

using QuantumCumulants
using OrdinaryDiffEq, ModelingToolkitBase
using Plots

hc = FockSpace(:cavity)
hm = FockSpace(:motion)
h = hc ⊗ hm

@qnumbers a::Destroy(h, 1) b::Destroy(h, 2)

@variables Δ ωm E G κ

H = -Δ * a' * a + ωm * b' * b + G * a' * a * (b + b') + E * (a + a')

J = [a]
rates = [κ]

ops = [a' * a, b' * b]
eqs = meanfield(ops, H, J; rates = rates, order = 2, simplify = false)
complete!(eqs; simplify = false)

sys = to_system(eqs; name = :optomech)
sys_c = mtkcompile(sys)

u0 = initial_values(
    eqs;
    defaults = Dict(average(b' * b) => 4.0e6 + 0im)
)

p0 = Dict(Δ => -10.0, ωm => 1.0, E => 200.0, G => 0.0125, κ => 20.0)
prob = ODEProblem(sys_c, merge(u0, p0), (0.0, 60000.0))
sol = solve(prob, Tsit5())

ts = sol.t
phonons = real.(get_solution(sol, b' * b, eqs).(ts))
T = 7.5e-5 .* phonons
photons = real.(get_solution(sol, a' * a, eqs).(ts))

p1 = plot(ts, T, ylabel = "T in K", legend = false)
p2 = plot(ts, photons, xlabel = "t⋅ωm", ylabel = "⟨a⁺a⟩", legend = false)
plot(p1, p2, layout = (2, 1), size = (650, 400))
