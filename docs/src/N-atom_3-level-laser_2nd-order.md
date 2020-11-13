# N-atom 3-level laser

This example describes a simple second order laser system consisting of $N$(=2) 3-level atoms coupled to a single mode cavity. An auxiliary state $|3\rangle$, which decays fast into the upper lasing state $|2\rangle$, is coherently pumped to achieve population inversion on the lasing transition $|1\rangle \leftrightarrow |2\rangle$. The Hamiltonian of this system is given by

$H = -\Delta_c a^\dagger a  - \sum_{i=1}^N \left[ \Delta_3^i \sigma^{33}_i  + g_i (a^\dagger \sigma^{12}_i + a\sigma^{21}_i) + \Omega_i (\sigma^{31}_i + \sigma^{13}_i) \right] $

Including dissipative processes as e.g. the atomic decay or photon losses through the cavity mirrors makes it an open quantum system. In the Schrödinger picture we calulate the dynamic of such open quantum systems with a density matrix $\rho$ following the master equation (see e.g. https://docs.qojulia.org/)

$\frac{d}{dt} \rho = - \frac{i}{\hbar} \left[ H, \rho \right] + \mathcal{L}[\rho]$,

with $\mathcal{L}[\rho] = \frac{\gamma}{2} (2 J \rho J^\dagger - J^\dagger J \rho - \rho J^\dagger J)$ the Liouvillian superoperator in standard Lindblad form for a dissipativ process with jump operator $J$ and rate $R$.

With Qumulants.jl we describe the system dynamics with averages in the Heisenberg picture. In the Heisenberg picture we replace the master equation by the quantum Langevin equation and the density matrix by operators. Since we are only interested in averages of operators we can immediately omit the stochastic terms of the quantum Langevin equation. Thus we get the following equation for the time evolution of a system operator average $\langle O \rangle$ (if $O$ is not explicitly time depending)

$\frac{d}{dt} \langle O \rangle = \frac{i}{\hbar} \left[ H, O \right] + \bar{\mathcal{L}}[O]$.

The Liovillian superoperator in the Heisenber picture $\bar{\mathcal{L}}[O]$ is almost the same as in the Schrödinger picture, except that $J$ and $J^\dagger$ are swapped in the first term, this means $\bar{\mathcal{L}}[O] = \frac{\gamma}{2} (2 J^\dagger O J - J^\dagger J O - O J^\dagger J)$, for a dissipativ process with jump operator $J$ and rate $R$.

For our system we have four different dissipative processes with the jump operators $a$, $\sigma^{12}_i$, $\sigma^{13}_i$ and $\sigma^{23}_i$ and corresponding rates $\kappa$, $\Gamma^i_{12}$, $\Gamma^i_{13}$ and $\Gamma^i_{23}$, respectively.

We start by loading the needed packages.


```julia
using Qumulants
using OrdinaryDiffEq
using Plots;
```

Then we define the symbolic parameters of the system, the Hilbertspace and the necessary operators. We define a atomic transition operator function $\sigma(i,j,k)$ for the transition from $|j \rangle$ to $|i \rangle$ of atom $k$. Since we only have one FockSpace we do not need to specify the Hilbertspace where the Destroy operator acts on. For the different atomic transitions we need to specify this, since there are more than one NLevelSpace. This is done by the "aon" field of the operator. In our case it is the $k+1$ in the Transition function.


```julia
# Parameters
N = 2 #number of atoms
κ, g, Γ23, Γ13, Γ12, Ω, Δc, Δ3 = parameters("κ g Γ_{23} Γ_{13} Γ_{12} Ω Δ_c Δ_3")

# Hilbertspace
hf = FockSpace(:cavity)
ha = ⊗([NLevelSpace(Symbol(:atom,i),3) for i=1:N]...)
h = hf ⊗ ha

# Operators
a = Destroy(h,:a)
σ(i,j,k) = Transition(h,Symbol("σ_{$k}"),i,j,k+1);
```

Now we create he Hamiltonian and the Jumps with the corresponding rates of our laser system. We assume here that all atoms are identical.


```julia
# Hamiltonian
H = -Δc*a'a + sum(g*(a'*σ(1,2,i) + a*σ(2,1,i)) for i=1:N) + sum(Ω*(σ(3,1,i) + σ(1,3,i)) for i=1:N) - sum(Δ3*σ(3,3,i) for i=1:N)

# Jumps
J = [a;[σ(1,2,i) for i=1:N];[σ(1,3,i) for i=1:N];[σ(2,3,i) for i=1:N]]

# Rates
rates = [κ;[Γ12 for i=1:N];[Γ13 for i=1:N];[Γ23 for i=1:N]];
```

Later we will complete the system automatically, which has the disadvantage that the equations are not ordered. Therefore we define a list of interseting operators, which we want to use later. Note that at least one operator(-product) is needed, which should have the order of the desired average-order. We derive the equations for these operators average them and automatically complete the system of equations.


```julia
# list of operators
ops = [a'a, σ(2,2,1), σ(3,3,1)]

he = heisenberg(ops,H,J; rates=rates)
he_avg_ = average(he,2) #second order average
he_avg = complete(he_avg_; multithread=true); #automatically complete the system
```


```julia
he_avg_
```




\begin{align}
\frac{d}{dt} \langle a^\dagger  a\rangle  =& -1.0 i g \left( \langle a^\dagger  \sigma_{1}^{12}\rangle  + \langle a^\dagger  \sigma_{2}^{12}\rangle  \right) + 1.0 i g \left( \langle a  \sigma_{1}^{21}\rangle  + \langle a  \sigma_{2}^{21}\rangle  \right) -1.0 \kappa \langle a^\dagger  a\rangle  \\
\frac{d}{dt} \langle \sigma_{1}^{22}\rangle  =& \langle \sigma_{1}^{33}\rangle  \Gamma_{{23}} -1.0 \langle \sigma_{1}^{22}\rangle  \Gamma_{{12}} + 1.0 i g \langle a^\dagger  \sigma_{1}^{12}\rangle  -1.0 i g \langle a  \sigma_{1}^{21}\rangle  \\
\frac{d}{dt} \langle \sigma_{1}^{33}\rangle  =& -1.0 i \langle \sigma_{1}^{31}\rangle  \Omega + 1.0 i \langle \sigma_{1}^{13}\rangle  \Omega -1.0 \langle \sigma_{1}^{33}\rangle  \left( \Gamma_{{13}} + \Gamma_{{23}} \right)
\end{align}




To calculate the time evolution we create a set of coupled ordinary differential equations which can be solved by DifferentialEquations.jl.


```julia
# list of symbolic parameters
ps = (g, Γ23, Γ13, Γ12, Ω, Δc, Δ3, κ)

meta_f = build_ode(he_avg, ps)

# function for DifferentialEquations.jl
f = Meta.eval(meta_f);
```

Finally we calculate the time evoultion after definig an initial state and numerical values for the parameters.


```julia
# initial state
u0 = zeros(ComplexF64, length(he_avg))

Γ12n = 1.0
Γ23n = 20Γ12n
Γ13n = 2Γ12n
Ωn = 5Γ13n
gn = 2Γ12n
Δcn = 0.0
Δ3n = 0.0
κn = 0.5Γ12n

p0 = (gn, Γ23n, Γ13n, Γ12n, Ωn, Δcn, Δ3n, κn)
tend = 10.0/κn

prob = ODEProblem(f,u0,(0.0,tend),p0)
sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8);
```

We plot the average photon number and the Population inversion.


```julia
n_t = real.(getindex.(sol.u, 1))
σ22m11_t = real.(2*getindex.(sol.u, 2) .+ getindex.(sol.u, 2) .-1 );#σ11 + σ22 + σ33 = 𝟙
```


```julia
plot(sol.t, n_t, xlabel="tΓ₁₂", ylabel="⟨a⁺a⟩", legend = false)
savefig("photon-number.svg");
```


```julia
plot(sol.t, σ22m11_t, xlabel="tΓ₁₂", ylabel="⟨σ22⟩ - ⟨σ11⟩", legend = false)
savefig("population-inversion.svg");
```


```julia

```
![](photon-number.svg)
![](population-inversion.svg)
