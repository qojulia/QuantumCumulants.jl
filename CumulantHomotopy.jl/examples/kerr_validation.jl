# # Phase 0 validation: parametric Kerr oscillator
#
# This notebook demonstrates the full Phase 0 concept of CumulantHomotopy.jl on a
# parametric Kerr oscillator (SPEC §11.2): generate the stationary cumulant hierarchy at
# several truncation orders, export each to a minimal real polynomial system, enumerate all
# steady states with HomotopyContinuation.jl, and compare the physical branch against the
# exact stationary moments of the full Liouvillian.

using QuantumCumulants
using Symbolics
using CumulantHomotopy
import QuantumCumulants as QC
using QuantumToolbox

# ## The model
#
# A single mode with detuning `Δ`, Kerr nonlinearity `U`, two-photon (parametric) drive `G`,
# and photon loss at rate `κ`:
#
# ```math
# H = Δ\, a^\dagger a + \tfrac{U}{2}\, a^\dagger a^\dagger a a + \tfrac{G}{2} (a a + a^\dagger a^\dagger), \qquad J = \{a\}.
# ```

@variables Δ::Real U::Real κ::Real G::Real
h = FockSpace(:cavity)
@qnumbers a::Destroy(h)
H = Δ * a' * a + U * a' * a' * a * a /2 + G * (a*a + a'*a') /2
J = [a]

# A deeply quantum working point: below the parametric threshold (`G < κ/2`) with strong Kerr
# (`U ≈ κ`) and sub-single-photon occupation. Here mean field predicts an empty cavity, so the
# finite steady-state occupation is entirely a quantum-fluctuation effect that only the higher
# cumulant orders can recover.

parameters = Dict(Δ => 0.0, U => 1.0, κ => 1.0, G => 0.3)

# ## Exact reference (full Liouvillian)
#
# The exact stationary state is the steady state of the full Liouvillian in a truncated Fock
# space, using the *same* Hamiltonian and jump operator. QuantumToolbox.jl's `steadystate`
# solves it directly (sparse). This is the ground truth the hierarchy is compared against.

function exact_stationary(; Δ, U, κ, G, N = 60)
    a = destroy(N)
   H = -Δ  * a' * a + U * (a'^2 * a^2) / 2 - G * (a' * a' + a * a) / 2
    ρss = steadystate(H, [sqrt(κ) * a])
    return (; n = real(expect(a' * a, ρss)), a = expect(a, ρss))
end
exact = exact_stationary(; Δ = 0.0, U = 1.0, κ = 1.0, G = 0.3)
@show exact.n exact.a;

# ## Truncation orders to study
#
# Start at order 1 (pure mean field) and go up to order 6.

orders = 1:5
n_operator = QC.average(a' * a)
a_operator = QC.average(a)

# Photon number of a solution. Orders ≥ 2 track `⟨a'a⟩` directly; order 1 (mean field) does
# not, so fall back to the factorization `⟨a'a⟩ ≈ ⟨a'⟩⟨a⟩ = |⟨a⟩|²`.
function photon_number(sol)
    haskey(sol, n_operator) && return real(sol[n_operator])
    return abs2(sol[a_operator])
end

# ## Generate and close the hierarchies
#
# One completed `MeanfieldEquations` per order, kept for inspection.

hierarchies = Dict(
    order => complete(meanfield(a, H, J; rates = [κ], order = order))
    for order in orders)

# ## Export a real polynomial system per order
#
# Each hierarchy is realified to its own `StationaryPolynomialSystem` (SPEC §3.2) and saved
# separately, so it can be inspected or re-solved later.

systems = Dict(order => realify(hierarchies[order]) for order in orders)

systems[2]

# ## Enumerate steady states per order
#
# Solve each system with HomotopyContinuation.jl. `all_roots` keeps every complex root;
# `physical` keeps the real roots (the physical steady states).

# all_roots = Dict(
#     order => stationary_states(systems[order], parameters; only_physical = false, show_progress = false)
#     for order in orders)
physical = Dict(
    order => stationary_states(systems[order], parameters; show_progress = false)
    for order in orders)

# ## Compare the physical branch against the exact stationary state
#
# For each order, pick the physical root closest to the exact photon number and tabulate the
# error.

begin
println("order | #real vars | #physical | photon number n | |n - n_exact|")
println("------+------------+-----------+-----------------+--------------")
for order in orders
    ns = [photon_number(sol) for sol in physical[order]]
    n_best = isempty(ns) ? NaN : ns[argmin(abs.(ns .- exact.n))]
    println(
        "  ", order,
        "   |     ", lpad(length(systems[order].variables), 2),
        "     |     ", lpad(length(physical[order]), 2),
        "    |    ", rpad(round(n_best; digits = 6), 10),
        "   |  ", round(abs(n_best - exact.n); digits = 6),
    )
end
end

# ## Result
#
# The exact occupation is `n ≈ 0.087` at a point *below* the classical parametric threshold:
# mean field (order 1) predicts an empty cavity (`n = 0`), so the entire steady-state
# occupation is a quantum effect. The cumulant hierarchy recovers it as the order grows:
#
# | order | best `n` | `|n − n_exact|` |
# |-------|----------|-----------------|
# | 1     | 0.000    | 8.7e-2          |
# | 2     | 0.068    | 1.8e-2          |
# | 3     | 0.068    | 1.8e-2          |
# | 4     | 0.083    | 4.3e-3          |
# | 5     | 0.083    | 4.3e-3          |
#
# Two structural features are visible. First, a **parity effect** (SPEC §6.4): the model has a
# `Z₂` symmetry (`a → -a`), so an odd order adds no new information over the even order below
# it, and convergence runs through the even orders `2 → 4 → …`. Second, **root proliferation**:
# the number of real roots grows steeply with order (1, 5, 5, 29, 71) and includes clearly
# unphysical ones (negative occupations). Here the physical branch is chosen by proximity to
# the exact answer; doing that selection *without* the exact answer, through physical
# diagnostics (positivity, held-out equations, order-lift defects) and continuation, is the
# subject of the later roadmap phases.
#
# This closes the Phase 0 loop: generate → realify → enumerate with HomotopyContinuation.jl →
# compare against exact stationary moments.
