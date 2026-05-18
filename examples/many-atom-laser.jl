# # Many-atom laser

# Second-order laser system consisting of $N$ three-level atoms coupled to a single
# mode cavity. An auxiliary state $|3\rangle$ is coherently pumped to invert the lasing
# transition $|1\rangle \leftrightarrow |2\rangle$. Hamiltonian:

# ```math
# H = -\Delta_c a^\dagger a + \sum_{i=1}^N \left[\Delta_3 \sigma_i^{33} - g(a^\dagger \sigma_i^{12} + a\sigma_i^{21}) - \Omega(\sigma_i^{31} + \sigma_i^{13})\right].
# ```

using QuantumCumulants
using ModelingToolkitBase, OrdinaryDiffEq
using Plots

N = 2 # number of atoms
@variables κ g Γ23 Γ13 Γ12 Ω Δc Δ3

hf = FockSpace(:cavity)
ha = ⊗([NLevelSpace(Symbol(:atom, i), 3) for i = 1:N]...)
h = hf ⊗ ha

a = Destroy(h, :a)
σ(i, j, k) = Transition(h, Symbol("σ_", k), i, j, k + 1)
nothing # hide

H =
    -Δc * a' * a +
    sum(g * (a' * σ(1, 2, i) + a * σ(2, 1, i)) for i = 1:N) +
    sum(Ω * (σ(3, 1, i) + σ(1, 3, i)) for i = 1:N) -
    sum(Δ3 * σ(3, 3, i) for i = 1:N)

J = [a; [σ(1, 2, i) for i = 1:N]; [σ(1, 3, i) for i = 1:N]; [σ(2, 3, i) for i = 1:N]]
rates = [κ; [Γ12 for i = 1:N]; [Γ13 for i = 1:N]; [Γ23 for i = 1:N]]
nothing # hide

ops = [a' * a, σ(2, 2, 1), σ(3, 3, 1)]
eqs = meanfield(ops, H, J; rates = rates, order = 2, simplify = false)
complete!(eqs; simplify = false)

sys = to_system(eqs; name = :laser)
sys_c = mtkcompile(sys)

u0 = initial_values(eqs)

Γ12n = 1.0
Γ23n = 20Γ12n
Γ13n = 2Γ12n
Ωn = 5Γ13n
gn = 2Γ12n
Δcn = 0.0
Δ3n = 0.0
κn = 0.5Γ12n

p0 = Dict(g => gn, Γ23 => Γ23n, Γ13 => Γ13n, Γ12 => Γ12n,
          Ω => Ωn, Δc => Δcn, Δ3 => Δ3n, κ => κn)
tend = 10.0 / κn

prob = ODEProblem(sys_c, merge(u0, p0), (0.0, tend))
sol = solve(prob, Tsit5(), reltol = 1e-8, abstol = 1e-8)
nothing # hide

n_t = real.(get_solution(sol, a' * a, eqs).(sol.t))
σ22 = real.(get_solution(sol, σ(2, 2, 1), eqs).(sol.t))
σ33 = real.(get_solution(sol, σ(3, 3, 1), eqs).(sol.t))
σ22m11_t = 2 .* σ22 .+ σ33 .- 1

p1 = plot(sol.t, n_t, xlabel = "tΓ₁₂", ylabel = "⟨a⁺a⟩", legend = false)
p2 = plot(sol.t, σ22m11_t, xlabel = "tΓ₁₂", ylabel = "⟨σ22⟩ - ⟨σ11⟩", legend = false)
plot(p1, p2, layout = (1, 2), size = (800, 300))
