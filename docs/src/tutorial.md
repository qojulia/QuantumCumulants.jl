# Tutorial

The basic usage is probably best illustrated with a brief example. In the following, we solve a simple model for a single-atom laser.

We start by loading the package, defining some symbolic parameters and the photonic annihilation operator `a` as well as the atomic transition operator `σ`, which denotes a transition from level `j` to level `i` as `σ(i,j)`. This allows us to quickly write down the Hamiltonian and the collapse operators of the system with their corresponding decay rates.

```@example tutorial
using Latexify # hide
set_default(double_linebreak=true) # hide
using QuantumCumulants

# Define parameters
@variables Δ::Real g::Real γ::Real κ::Real ν::Real

# Define hilbert space
hf = FockSpace(:cavity)
ha = NLevelSpace(:atom, (:g, :e))
h = hf ⊗ ha

# Define the fundamental operators
@qnumbers a::Destroy(h)
σ(i, j) = Transition(h, :σ, i, j)

# Hamiltonian
H = Δ*a'*a + g*(a'*σ(:g, :e) + a*σ(:e, :g))

# Collapse operators
J = [a, σ(:g, :e), σ(:e, :g)]
rates = [κ, γ, ν]
nothing # hide
```

Now, we define a list of operators of which we want to compute the mean-field equations. We will only consider products of two operators. This is because later we will compute the dynamics of the system up to second order.

```@example tutorial
# Derive a set of equations
ops = [a'*a, σ(:e, :e), a'*σ(:g, :e)]
eqs = meanfield(ops, H, J; rates=rates)
```

To obtain a closed set of equations, we expand higher-order products to second order.

```@example tutorial
# Expand the above equations to second order
eqs_expanded = cumulant_expansion(eqs, 2)
```

The first-order contributions are always zero and can therefore be neglected. You can try adding `a` and `σ(:g,:e)` to the list of operators `ops` in order to see that yourself. Even more conveniently, [`complete`](@ref) automatically finds all missing averages and derives the corresponding equations, giving a closed system in one call.

```@example tutorial
# Close the system by deriving equations for any missing averages
eqs_full = complete(eqs_expanded)
```

Finally, we convert the [`MeanfieldEquations`](@ref) to a `System` from [ModelingToolkitBase](https://github.com/SciML/ModelingToolkitBase.jl), which can be solved numerically with [OrdinaryDiffEq](https://github.com/JuliaDiffEq/OrdinaryDiffEq.jl).

```@example tutorial
# Generate a ModelingToolkitBase System
using ModelingToolkitBase
sys = mtkcompile(System(eqs_full; name=:laser))

# Solve the system using the OrdinaryDiffEq package
using OrdinaryDiffEqTsit5
u0 = zeros(ComplexF64, length(eqs_full.states))
p0 = [Δ => 0.0, g => 1.5, γ => 0.25, κ => 1.0, ν => 4.0]
prob = ODEProblem(sys, merge(initial_values(eqs_full, u0), Dict(p0)), (0.0, 10.0))
sol = solve(prob, Tsit5())
nothing # hide
```

Numeric trajectories of any operator-average are recovered with [`get_solution`](@ref), which substitutes the symbolic average into the compiled solution and returns a callable in `t`:

```@example tutorial
using Plots
ts = range(0.0, 10.0; length=200)
n  = real.(get_solution(sol, a'*a, eqs_full).(ts))
pe = real.(get_solution(sol, σ(:e, :e), eqs_full).(ts))
plot(ts, n, label="Photon number", xlabel="t")
plot!(ts, pe, label="Excited state population")
savefig("tutorial.svg") # hide
nothing # hide
```

![Photon number and excited state population](tutorial.svg)
