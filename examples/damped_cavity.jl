# # Damped Cavity
#
# This is a minimal example demonstrating the v1 workflow for `QuantumCumulants.jl`.
# We derive the equation of motion for the amplitude `⟨a⟩` of a damped cavity, build
# an MTKBase `System`, and integrate it numerically.
#
# The Hamiltonian is `H = ω a^† a`, and the cavity loses photons at rate `κ` via the
# collapse operator `a`. The analytical solution is `⟨a⟩(t) = ⟨a⟩(0) exp((-iω - κ/2) t)`.

using QuantumCumulants
using ModelingToolkitBase
using OrdinaryDiffEq
using Plots

# Hilbert space and operators
hc = FockSpace(:cavity)
@qnumbers a::Destroy(hc)
@variables ω κ

# Build the meanfield equations
H = ω * a' * a
eqs = meanfield([a], H, [a]; rates = [κ])
complete!(eqs)

# Convert to an MTKBase System and compile
sys = to_system(eqs; name = :cavity)
sys_c = mtkcompile(sys)

# Set initial conditions and parameters, solve
u0 = initial_values(eqs; defaults = Dict(average(a) => 1.0 + 0.0im))
p  = Dict(ω => 2.0, κ => 0.5)

prob = ODEProblem(sys_c, merge(u0, p), (0.0, 5.0))
sol  = solve(prob, Tsit5(); reltol = 1e-10, abstol = 1e-12)

# Plot vs analytical solution
ts   = range(0.0, 5.0; length = 200)
num  = real.(get_solution(sol, a, eqs).(ts))
ana  = real.(exp.((-im * 2.0 - 0.5 / 2) .* ts))
plot(ts, [num ana]; label = ["Re ⟨a⟩ (numeric)" "Re ⟨a⟩ (analytical)"],
     xlabel = "t", linewidth = 2)
