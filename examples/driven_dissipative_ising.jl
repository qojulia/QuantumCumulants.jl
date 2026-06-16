# # Driven-dissipative transverse-field Ising chain

# In this example we use **Pauli operators** to study a one-dimensional transverse-field Ising
# chain that is also subject to dissipation. Pauli operators are the natural language for
# two-level systems (qubits, spin-1/2), and because they are Hermitian the resulting cumulant
# equations carry real coefficients.

# The Hamiltonian of a chain of $N$ spins reads

# ```math
# H = -J \sum_{i=1}^{N-1} \sigma^z_i \sigma^z_{i+1} - h_x \sum_{i=1}^{N} \sigma^x_i,
# ```

# where $J$ is the nearest-neighbour Ising coupling and $h_x$ is the strength of a transverse
# field. On top of the coherent dynamics, each spin decays towards its ground state with rate
# $\gamma$ through the collapse operator

# ```math
# \sigma^-_i = \tfrac{1}{2}\left(\sigma^x_i - i\,\sigma^y_i\right) = |g\rangle\langle e|_i.
# ```

# The competition between the entangling Ising term, the transverse drive and the local
# dissipation builds up correlations between the spins. A plain mean-field treatment
# (first-order cumulants) factorises all of these correlations and gives qualitatively wrong
# dynamics. The cumulant expansion is a *systematically improvable* hierarchy, though: by going
# to higher order we keep ever larger groups of correlated spins and converge towards the exact
# result. We will see this explicitly by comparing orders 1, 2 and 3 against the exact master
# equation for a chain small enough to solve directly.

# As always, we start by loading the packages and setting up the operators.

using QuantumCumulants
using ModelingToolkitBase
using ModelingToolkitBase: unknowns
using OrdinaryDiffEqLowOrderRK
using QuantumOptics
using Plots

N = 6 # number of spins in the chain

# Each site lives in its own `PauliSpace`; the full Hilbert space is their tensor product. The
# last argument of `Pauli` selects the site, and the `axis` argument selects the component
# ($1 = x$, $2 = y$, $3 = z$).

h = ⊗([PauliSpace(Symbol(:spin, i)) for i in 1:N]...)
σx(i) = Pauli(h, :σ, 1, i)
σy(i) = Pauli(h, :σ, 2, i)
σz(i) = Pauli(h, :σ, 3, i)
σm(i) = (σx(i) - 1im * σy(i)) / 2 # lowering operator |g><e| on site i

@variables J hx γ

H = -J * sum(σz(i) * σz(i + 1) for i in 1:(N - 1)) - hx * sum(σx(i) for i in 1:N)
c_ops = [σm(i) for i in 1:N]
nothing # hide

# We derive the equations of motion for the magnetisation $\langle\sigma^z_i\rangle$ of every
# spin and complete the set at a given cumulant order. `order = 1` reproduces the bare
# mean-field equations; higher orders progressively add two-spin, three-spin, ... correlations.

function magnetization_eqs(order)
    eqs = meanfield([σz(i) for i in 1:N], H, c_ops; rates = [γ for i in 1:N], order = order)
    complete!(eqs)
    return eqs
end

eqs = [magnetization_eqs(order) for order in 1:3]
nothing # hide

# The number of equations grows quickly with the cumulant order, since we track all correlations
# up to that many spins.

[length(e) for e in eqs]

# To solve numerically we build a `System`, choose the physical parameters, and set the initial
# state. We start with every spin pointing up, $|e e \dots e\rangle$. The initial averages are
# obtained directly from this product state with `initial_values(eqs, ψ₀)`, which converts the
# [QuantumOptics.jl](https://qojulia.org) state into the correct moments.

p = Dict(J => 1.0, hx => 1.0, γ => 0.2)
ψ0 = tensor([spinup(SpinBasis(1 // 2)) for _ in 1:N]...)

function solve_cumulant(eqs; tend = 20.0)
    sys = mtkcompile(System(eqs; name = :sys))
    u0 = initial_values(eqs, ψ0)
    dict = parameter_map(sys, merge(Dict(unknowns(sys) .=> u0), p))
    # `eval_expression = true` compiles the right-hand side via `eval` instead of a
    # RuntimeGeneratedFunction, which is noticeably faster to compile for the larger
    # higher-order systems; `build_initializeprob = false` skips the (here redundant)
    # initialization problem since `initial_values` already provides every state.
    prob = ODEProblem(sys, dict, (0.0, tend); eval_expression = true, build_initializeprob = false)
    sol = solve(prob, RK4(), saveat = 0.1)
    mz = [real.(get_solution(sol, σz(i), eqs).(sol.t)) for i in 1:N]
    return sol.t, mz
end

sols = [solve_cumulant(e) for e in eqs]
nothing # hide

# For a chain of only six spins we can also integrate the full Lindblad master equation
# directly in [QuantumOptics.jl](https://qojulia.org) and use it as the ground truth. We map
# the symbolic operators onto numeric ones acting on a tensor product of spin-1/2 bases with
# [`to_numeric`](@ref).

b = tensor([SpinBasis(1 // 2) for _ in 1:N]...)
Hn = to_numeric(-p[J] * sum(σz(i) * σz(i + 1) for i in 1:(N - 1)) - p[hx] * sum(σx(i) for i in 1:N), b)
Jn = [to_numeric(σm(i), b) for i in 1:N]

tout, ρt = timeevolution.master(collect(0:0.1:20.0), ψ0, Hn, Jn; rates = [p[γ] for _ in 1:N])
mz_exact = [real.(expect.(Ref(to_numeric(σz(i), b)), ρt)) for i in 1:N]
nothing # hide

# Let us compare the magnetisation of the central spin across cumulant orders.

c = N ÷ 2
colors = [:firebrick, :steelblue, :seagreen]
plt = plot(
    tout,
    mz_exact[c];
    label = "master equation",
    lw = 3,
    color = :black,
    xlabel = "γt",
    ylabel = "⟨σz⟩ (central spin)",
    legend = :topright,
)
for order in 1:3
    t, mz = sols[order]
    plot!(plt, t, mz[c]; label = "cumulants (order $order)", lw = 2, color = colors[order])
end
plt

# The story is the convergence of the hierarchy. Mean-field (order 1) shows large, persistent
# oscillations: by factorising the spin-spin correlations it badly misjudges how the Ising
# coupling damps the collective dynamics. Order 2 already removes most of the spurious
# oscillation, and order 3 sits essentially on top of the exact master equation. This is exactly
# the regime where higher cumulant orders pay off: strong coupling ($J \sim h_x$) with only weak
# dissipation ($\gamma \ll J$) builds up correlations that low orders cannot capture.

# Finally, we visualise the full space-time dynamics of the chain from the third-order solution.
# Each spin starts up (bright) and relaxes towards the steady state, with the transverse field
# and Ising coupling imprinting a transient pattern along the chain.

t3, mz3 = sols[3]
heatmap(
    t3,
    1:N,
    reduce(hcat, mz3)';
    xlabel = "γt",
    ylabel = "site i",
    colorbar_title = "⟨σz_i⟩",
    color = :balance,
    clims = (-1, 1),
)
