# # Noisy excitation transport

# Energy transport in a one-dimensional chain of atoms, where only the first atom is
# driven. Position-dependent dipolar interactions move excitations from neighbor to
# neighbor.

using QuantumCumulants
using ModelingToolkitBase, OrdinaryDiffEq
using Plots

N = 5 # Hilbert space for N atoms
h = ⊗([NLevelSpace(Symbol(:atom, i), (:g, :e)) for i = 1:N]...)

σ(i, j, k) = Transition(h, Symbol(:σ_, k), i, j, k)

@variables Ω γ Δ J0
x = [first(@variables $(Symbol("x_$i"))) for i = 1:N]
J(xᵢ, xⱼ) = J0 / abs(xᵢ - xⱼ)^3

H =
    -Δ * sum(σ(:e, :e, k) for k = 1:N) +
    Ω * (σ(:e, :g, 1) + σ(:g, :e, 1)) +
    sum(
        J(x[k], x[k + 1]) *
        (σ(:e, :g, k) * σ(:g, :e, k + 1) + σ(:g, :e, k) * σ(:e, :g, k + 1)) for
        k = 1:(N - 1)
    )

c_ops = [σ(:g, :e, k) for k = 1:N]

eqs = meanfield(σ(:g, :e, 1), H, c_ops;
                rates = [γ for i = 1:N], order = 2, simplify = false)
complete!(eqs; simplify = false)

sys = to_system(eqs; name = :sys)
sys_c = mtkcompile(sys)

d = 0.75
x0 = [d * (k - 1) for k = 1:N]
p = merge(
    Dict(γ => 1.0, Δ => 0.0, Ω => 2.0, J0 => 1.25),
    Dict(x .=> x0),
)

u0 = initial_values(eqs)
prob = ODEProblem(sys_c, merge(u0, p), (0.0, 15.0))

sol = solve(prob, Tsit5())

pop1 = real.(get_solution(sol, σ(:e, :e, 1), eqs).(sol.t))
popN = real.(get_solution(sol, σ(:e, :e, N), eqs).(sol.t))
graph = plot(sol.t, pop1, label = "Driven atom", xlabel = "γt",
             ylabel = "Excited state population")
plot!(graph, sol.t, popN, label = "End of chain")
