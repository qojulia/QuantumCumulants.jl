```@meta
EditURL = "../../../examples/optomechanical-cooling.jl"
```

# Optomechanical Cooling

In this example, we show how to implement a cooling scheme based on radiation pressure coupling of light to a mechanical oscillator, such as a membrane. The oscillator is placed inside an optical cavity. The cavity is driven by a laser and the resulting radiation pressure of the cavity field effectively couples the photons in the cavity mode to the vibrational phonons of the mechanical oscillator mode. This model is based on the one studied in [C. Genes, et. al., Phys. Rev. A 77, 033804 (2008)](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.77.033804), and the Hamiltonian reads

$H = -\hbar\Delta a^\dagger a + \hbar\omega_m b^\dagger b + \hbar Ga^\dagger a \left(b + b^\dagger\right) + \hbar E \left(a + a^\dagger\right),$

where $\Delta = \omega_\ell - \omega_c$ is the detuning between the driving laser ($\omega_\ell$) and the cavity ($\omega_c$). The amplitude of the laser is denoted by $E$, the resonance frequency of the mechanical oscillator by $\omega_m$, and the radiation pressure coupling is given by $G$. Additionally, photons leak out of the cavity at a rate $\kappa$.
We start by loading the needed packages and specifying the model.

````@example optomechanical-cooling
using QuantumCumulants
using OrdinaryDiffEqLowOrderRK, ModelingToolkitBase
using Plots

hc = FockSpace(:cavity) # Hilbertspace
hm = FockSpace(:motion)
h = hc ⊗ hm

@qnumbers a::Destroy(h, 1) b::Destroy(h, 2) # Operators


@variables Δ ωm E G κ # Parameters


H = -Δ * a' * a + ωm * b' * b + G * a' * a * (b + b') + E * (a + a') # Hamiltonian


J = [a] # Jump operators & rates
rates = [κ]
nothing # hide
````

We are specifically interested in the average number of photons $\langle a^\dagger a \rangle$ and phonons $\langle b^\dagger b \rangle$. Thus, we first derive the equations for these two averages. We restrict our description to a second order cumulant expansion.

````@example optomechanical-cooling
ops = [a' * a, b' * b] # Derive equations
eqs = meanfield(ops, H, J; rates = rates, order = 2)
````

To get a closed set of equations we automatically complete the system.

````@example optomechanical-cooling
eqs_completed = complete!(deepcopy(eqs)) # Complete equations
````

To calculate the dynamics we create a system of ordinary differential equations, which can be used by [DifferentialEquations.jl](https://diffeq.sciml.ai/stable/).

````@example optomechanical-cooling
sys = System(eqs_completed; name = :sys)
sys_c = mtkcompile(sys)
nothing # hide
````

Finally, we need to define the numerical parameters and the initial state of the system. We will consider the membrane at room temperature. Its vibrational mode is in a thermal state with an average number of phonons that can be estimated from $k_B T = n_\mathrm{vib}\hbar \omega_m$. If the resonator has a resonance frequency of $\omega_m = 10\mathrm{MHz}$, then the number of phonons at room temperature ($T\approx 300K$) is approximately $n_\mathrm{vib} \approx 4\times 10^6$.

````@example optomechanical-cooling
u0 = initial_values(eqs_completed; defaults = Dict(average(b' * b) => 4.0e6 + 0im)) # Initial state (4e6 phonons)

p0 = Dict{Num, ComplexF64}(Δ => -10.0 + 0im, ωm => 1.0 + 0im, E => 200.0 + 0im, G => 0.0125 + 0im, κ => 20.0 + 0im) # System parameters
prob = ODEProblem(sys_c, merge(u0, p0), (0.0, 60000.0))
sol = solve(prob, RK4())
nothing # hide
````

````@example optomechanical-cooling
t = real.(sol.t) # Plot results
phonons = real.(get_solution(sol, b'b, eqs_completed).(sol.t))
T = 7.5e-5 * phonons
photons = real.(get_solution(sol, a'a, eqs_completed).(sol.t))

p1 = plot(t, T, ylabel = "T in K", legend = false)
p2 = plot(t, photons, xlabel = "t⋅ωm", ylabel = "⟨a⁺a⟩", legend = false)
plot(p1, p2, layout = (2, 1), size = (650, 400))
````

## Package versions

These results were obtained using the following versions:

````@example optomechanical-cooling
using InteractiveUtils
versioninfo()

using Pkg
Pkg.status(
    ["QuantumCumulants", "OrdinaryDiffEqLowOrderRK", "ModelingToolkitBase", "Plots"],
    mode = PKGMODE_MANIFEST,
)
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

