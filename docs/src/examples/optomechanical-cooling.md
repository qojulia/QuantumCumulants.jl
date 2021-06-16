# Optomechanical Cooling

In this example, we show how to implement a cooling scheme based on radiation pressure coupling of light to a mechanical oscillator, such as a membrane. The oscillator is placed inside an optical cavity. The cavity is driven by a laser and the resulting radiation pressure of the cavity field effectively couples the photons in the cavity mode to the vibrational phonons of the mechanical oscillator mode. This model is based on the one studied in [C. Genes, et. al., Phys. Rev. A 77, 033804 (2008)](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.77.033804), and the Hamiltonian reads

$H = -\hbar\Delta a^\dagger a + \hbar\omega_m b^\dagger b + \hbar Ga^\dagger a \left(b + b^\dagger\right) + \hbar E \left(a + a^\dagger\right),$

where $\Delta = \omega_\ell - \omega_c$ is the detuning between the driving laser ($\omega_\ell$) and the cavity ($\omega_c$). The amplitude of the laser is denoted by $E$, the resonance frequency of the mechanical oscillator by $\omega_m$, and the radiation pressure coupling is given by $G$. Additionally, photons leak out of the cavity at a rate $\kappa$.
We start by loading the needed packages and specifying the model.


```@example optomechanics
using QuantumCumulants
using OrdinaryDiffEq, ModelingToolkit
using Plots

# Hilbertspace
hc = FockSpace(:cavity)
hm = FockSpace(:motion)
h = hc ⊗ hm

# Operators
@qnumbers a::Destroy(h,1) b::Destroy(h,2)

# Parameters
@cnumbers Δ ωm E G κ

# Hamiltonian
H = -Δ*a'*a + ωm*b'*b + G*a'*a*(b + b') + E*(a + a')

# Jump operators & rates
J = [a]
rates = [κ]
nothing # hide
```

We are specifically interested in the average number of photons $\langle a^\dagger a \rangle$ and phonons $\langle b^\dagger b \rangle$. Thus we first derive the equations for these two averages.


```@example optomechanics
# Derive equations
ops = [a'*a, b'*b]
eqs = meanfield(ops,H,J;rates=rates,order=2)
```


To get a closed set of equations we automatically complete the system, where we restrict our description to second order cumulant expansion.


```@example optomechanics
# Complete equations
eqs_completed = complete(eqs)
```


To calculate the dynamics we create a system of ordinary differential equations, which can be used by [DifferentialEquations.jl](https://diffeq.sciml.ai/stable/).


```@example optomechanics
sys = ODESystem(eqs_completed)
nothing # hide
```

Finally we need to define the numerical parameters and the initial value of the system. We will consider the membrane at room temperature. Its vibrational mode is in a thermal state with an average number of phonons that can be estimated from $k_B T = n_\mathrm{vib}\hbar \omega_m$. If the resonator has a resonance frequency of $\omega_m = 10\mathrm{MHz}$, then the number of phonons at room temperature ($T\approx 300K$) is approximately $n_\mathrm{vib} \approx 4\times 10^6$.


```@example optomechanics
# Initial state
u0 = zeros(ComplexF64, length(eqs_completed))
u0[2] = 4e6 # Initial number of phonons
# System parameters
p0 = (Δ=>-10, ωm=>1, E=>200, G=>0.0125, κ=>20)
prob = ODEProblem(sys,u0,(0.0,60000),p0)
sol = solve(prob,RK4())
nothing # hide
```


```@example optomechanics
# Plot results
t = real.(sol.t)
phonons = real.(sol[b'b])
T = 7.5e-5*phonons
photons = real.(sol[a'a])

p1 = plot(t, T, ylabel="T in K", legend=false)
p2 = plot(t, photons, xlabel="t⋅ωm", ylabel="⟨a⁺a⟩", legend=false)
plot(p1, p2, layout=(2,1), size=(500,400))
savefig("opto_cooling.svg") # hide
```


![svg](opto_cooling.svg)
