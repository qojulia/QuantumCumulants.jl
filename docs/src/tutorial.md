# Tutorial

The basic usage is probably best illustrated with a brief example. In the following, we solve a simple model for a single-atom laser.

We start by loading the package, defining some symbolic parameters and the photonic annihilation operator `a` as well as the atomic transition operator `σ`, which denotes a transition from level `j` to level `i` as `σ(i,j)`. This allows us to quickly write down the Hamiltonian and the collapse operators of the system with their corresponding decay rates.

```@example tutorial
using Latexify # hide
set_default(double_linebreak=true) # hide
using Qumulants

# Define parameters
@cnumbers Δ g γ κ ν

# Define hilbert space
hf = FockSpace(:cavity)
ha = NLevelSpace(:atom,(:g,:e))
h = hf ⊗ ha

# Define the fundamental operators
@qnumbers a::Destroy(h) σ::Transition(h)

# Hamiltonian
H = Δ*a'*a + g*(a'*σ(:g,:e) + a*σ(:e,:g))

# Collapse operators
J = [a,σ(:g,:e),σ(:e,:g)]
rates = [κ,γ,ν]
nothing # hide
```

Now, we define a list of operators of which we want to compute the Heisenberg equations. We will only consider products of two operators. This is because later we will compute the dynamics of the system up to second order.

```@example tutorial
# Derive a set of Heisenberg equations
ops = [a'*a,σ(:e,:e),a'*σ(:g,:e)]
he = heisenberg(ops,H,J;rates=rates)
```

The equations derived above are differential equations for operators. In order to convert them to *c*-number equations, we need to average over them. To obtain a closed set of equations, we expand higher-order products to second order.

```@example tutorial
# Average the above equations and expand to second order
he_avg = average(he,2)
```

The first-order contributions are always zero and can therefore be neglected. You can try adding `a` and `σ(:g,:e)` to the list of operators `ops` in order to see that yourself. Or, even more conveniently, you can use `complete(he_avg)`, which will automatically find all missing averages and compute the corresponding equations.

Here, though, we will proceed by finding the missing averages, and neglecting them as zero using the `substitute` function.

```@example tutorial
# Find the missing averages
missed = find_missing(he_avg)

# Substitute them as zero
subs = Dict(missed .=> 0)
he_nophase = qsimplify(substitute(he_avg, subs))
```

Finally, we can convert the [`HeisenbergEquation`](@ref) to an `ODESystem` as defined in [ModelingToolkit](https://github.com/SciML/ModelingToolkit.jl) which can be solved numerically with [OrdinaryDiffEq](https://github.com/JuliaDiffEq/OrdinaryDiffEq.jl).

```@example tutorial
# Generate an ODESystem
using ModelingToolkit
sys = ODESystem(he_nophase)

# Solve the system using the OrdinaryDiffEq package
using OrdinaryDiffEq
u0 = zeros(ComplexF64,length(ops))
p = (Δ, g, γ, κ, ν)
p0 = p .=> (0, 1.5, 0.25, 1, 4)
prob = ODEProblem(sys,u0,(0.0,10.0),p0)
sol = solve(prob,RK4())
nothing # hide
```

The photon number of our laser and the excited state population of the atom are now stored in the first two fields of `sol.u`.

```@example tutorial
using Plots
n = real.(getindex.(sol.u, 1))
pe = real.(getindex.(sol.u, 2))
plot(sol.t, n, label="Photon number", xlabel="t")
plot!(sol.t, pe, label="Excited state population")
savefig("tutorial.svg") # hide
nothing # hide
```

![Photon number and excited state population](tutorial.svg)
