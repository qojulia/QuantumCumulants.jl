# Superradiant Laser

Using symmetry properties of a system can reduce the number of needed equations dramatically. A common approximation for laser systems to handle sufficiently big atom numbers is to assume that several atoms in the system behave completely identically. This means all the identical atoms have the same averages.

In this example we describe a so-called superradiant laser, where we assume all atoms to be identical. This model has been described in [D. Meiser et al., Phys. Rev. Lett. 102, 163601 (2009):](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.102.163601) The Hamiltonian of this system is
```math
\begin{equation}
H = - \hbar \Delta a^\dagger a +  \hbar \sum\limits_{j=1}^{N}  g_j (a^\dagger \sigma^{12}_{j} + a \sigma^{21}_{j}) ,
\end{equation}
```
where $\Delta = \omega_a - \omega_c$ is the detuning between the cavity ($\omega_c$) and the atomic ($\omega_a$) resonance frequency, the atom cavity coupling of the $j$-th atom is denoted by $g_j$. Additionally there are dissipative processes in the system: Atoms are incoherently pumped with the rate $R$, they decay individually with the rate $\Gamma$ and are affected by individual atomic dephasing with the rate $\nu$. Furthermore, photons leak out of the system with the rate $\kappa$.

We start by loading the packages.


```@example superradiant_laser_indexed
using QuantumCumulants
using OrdinaryDiffEq, SteadyStateDiffEq, ModelingToolkit
using Plots
```

Due to the implementation of symbolic indices and sums we only need to define the Hilbert space for one atom, even though we will simulate a system for several thousand.
Creating an operator with an $\texttt{Index}$ is done with the constructor $\texttt{IndexedOperator}$.


```@example superradiant_laser_indexed
# Hilbertspace
hc = FockSpace(:cavity)
ha = NLevelSpace(:atom,2)
h = hc ⊗ ha

# operators
@qnumbers a::Destroy(h)
σ(α,β,i) = IndexedOperator(Transition(h, :σ, α, β),i)
nothing # hide
```

Now we define the indices and the parameters of the system. An $\texttt{Index}$ needs the system Hilbert space, a symbol, an upper bound and the specific Hilbert space of the indexed operator. $\texttt{IndexedVariable}$ creates indexed variables. Actually we wouldn't need indexed variable in this example, this is just for demonstration purposes.


```@example superradiant_laser_indexed
@cnumbers N Δ κ Γ R ν
g(i) = IndexedVariable(:g, i)

i = Index(h,:i,N,ha)
j = Index(h,:j,N,ha)
```


We define the Hamiltonian using symbolic sums and define the individual dissipative processes. For an indexed jump operator the (symbolic) sum is build in the Liouvillian.


```@example superradiant_laser_indexed
# Hamiltonian
H = -Δ*a'a + Σ(g(i)*( a'*σ(1,2,i) + a*σ(2,1,i) ),i)

# Jump operators with corresponding rates
J = [a, σ(1,2,i), σ(2,1,i), σ(2,2,i)]
rates = [κ, Γ, R, ν]
nothing # hide
```

First we want to derive the equation for $\langle a^\dagger a \rangle$ and $\langle \sigma_j^{22} \rangle$. Note that you can only use indices on the LHS which haven't been used for the Hamiltonian and the jumps.


```@example superradiant_laser_indexed
# Derive equations
ops = [a'*a, σ(2,2,j)]
eqs = meanfield(ops,H,J;rates=rates,order=2)
nothing # hide
```



```math
\begin{align}
\frac{d}{dt} \langle a^\dagger  a\rangle  =& 1 i \underset{i}{\overset{N}{\sum}} {g}_{i}  \langle a  {\sigma}_{i}^{{21}}\rangle  -1 i \underset{i}{\overset{N}{\sum}} {g}_{i}  \langle a^\dagger  {\sigma}_{i}^{{12}}\rangle  -1.0 \kappa \langle a^\dagger  a\rangle  \\
\frac{d}{dt} \langle {\sigma}_{j}^{{22}}\rangle  =& R -1.0 R \langle {\sigma}_{j}^{{22}}\rangle  + 1 i {g}_{j} \langle a^\dagger  {\sigma}_{j}^{{12}}\rangle  -1 i {g}_{j} \langle a  {\sigma}_{j}^{{21}}\rangle  -1.0 \Gamma \langle {\sigma}_{j}^{{22}}\rangle
\end{align}
```



To get a closed set of equations we automatically complete the system. Since this system is phase invariant we know that all averages with a phase are zero, therefore we exclude these terms with a filter function. To be able to dispatch on all kind of sums containing averages we defined the Union $\texttt{AvgSums}$.


```@example superradiant_laser_indexed
# custom filter function
φ(x::Average) = φ(x.arguments[1])
φ(::Destroy) = -1
φ(::Create) =1
φ(x::QTerm) = sum(map(φ, x.args_nc))
φ(x::Transition) = x.i - x.j
φ(x::IndexedOperator) = x.op.i - x.op.j
φ(x::SingleSum) = φ(x.term)
φ(x::AvgSums) = φ(arguments(x))
phase_invariant(x) = iszero(φ(x))

# Complete equations
eqs_c = complete(eqs; filter_func=phase_invariant)
nothing # hide
```



```math
\begin{align}
\frac{d}{dt} \langle a^\dagger  a\rangle  =& 1 i \underset{i}{\overset{N}{\sum}} {g}_{i}  \langle a  {\sigma}_{i}^{{21}}\rangle  -1 i \underset{i}{\overset{N}{\sum}} {g}_{i}  \langle a^\dagger  {\sigma}_{i}^{{12}}\rangle  -1.0 \kappa \langle a^\dagger  a\rangle  \\
\frac{d}{dt} \langle {\sigma}_{j}^{{22}}\rangle  =& R -1.0 R \langle {\sigma}_{j}^{{22}}\rangle  + 1 i {g}_{j} \langle a^\dagger  {\sigma}_{j}^{{12}}\rangle  -1 i {g}_{j} \langle a  {\sigma}_{j}^{{21}}\rangle  -1.0 \Gamma \langle {\sigma}_{j}^{{22}}\rangle  \\
\frac{d}{dt} \langle a^\dagger  {\sigma}_{j}^{{12}}\rangle  =& 1 i \underset{i{\ne}j}{\overset{N}{\sum}} {g}_{i}  \langle {\sigma}_{i}^{{21}}  {\sigma}_{j}^{{12}}\rangle  + 1 i {g}_{j} \langle {\sigma}_{j}^{{22}}\rangle  -1 i {g}_{j} \langle a^\dagger  a\rangle  -0.5 R \langle a^\dagger  {\sigma}_{j}^{{12}}\rangle  -0.5 \Gamma \langle a^\dagger  {\sigma}_{j}^{{12}}\rangle  -0.5 \kappa \langle a^\dagger  {\sigma}_{j}^{{12}}\rangle  -0.5 \nu \langle a^\dagger  {\sigma}_{j}^{{12}}\rangle  -1 i \Delta \langle a^\dagger  {\sigma}_{j}^{{12}}\rangle  + 2 i {g}_{j} \langle {\sigma}_{j}^{{22}}\rangle  \langle a^\dagger  a\rangle  \\
\frac{d}{dt} \langle {\sigma}_{j}^{{12}}  {\sigma}_{k}^{{21}}\rangle  =& \left( -1.0 R -1.0 \Gamma \right) \langle {\sigma}_{j}^{{12}}  {\sigma}_{k}^{{21}}\rangle  -1 i {g}_{j} \langle a  {\sigma}_{k}^{{21}}\rangle  + 1 i {g}_{k} \langle a^\dagger  {\sigma}_{j}^{{12}}\rangle  -1.0 \nu \langle {\sigma}_{j}^{{12}}  {\sigma}_{k}^{{21}}\rangle  + 2 i {g}_{j} \langle {\sigma}_{j}^{{22}}\rangle  \langle a  {\sigma}_{k}^{{21}}\rangle  -2 i {g}_{k} \langle {\sigma}_{k}^{{22}}\rangle  \langle a^\dagger  {\sigma}_{j}^{{12}}\rangle
\end{align}
```



As mentioned before, we assume that all atoms behave identically. This means that e.g. the excited state population is equal for all atoms, hence we only need to calculate it for the first $\langle \sigma^{22}_1 \rangle = \langle \sigma^{22}_j \rangle$. Furthermore, it is clear that a sum over $N$ identical objects can be replaced by $N$ times the object. The function $\texttt{scale()}$ uses these rules to simplify the equations.


```@example superradiant_laser_indexed
eqs_sc = scale(eqs_c)
nothing # hide
```



```math
\begin{align}
\frac{d}{dt} \langle a^\dagger  a\rangle  =& -1.0 \kappa \langle a^\dagger  a\rangle  -1 i N g_{1} \langle a^\dagger  {\sigma}_{1}^{{12}}\rangle  + 1 i N g_{1} \langle a  {\sigma}_{1}^{{21}}\rangle  \\
\frac{d}{dt} \langle {\sigma}_{1}^{{22}}\rangle  =& R -1.0 R \langle {\sigma}_{1}^{{22}}\rangle  -1.0 \Gamma \langle {\sigma}_{1}^{{22}}\rangle  + 1 i g_{1} \langle a^\dagger  {\sigma}_{1}^{{12}}\rangle  -1 i g_{1} \langle a  {\sigma}_{1}^{{21}}\rangle  \\
\frac{d}{dt} \langle a^\dagger  {\sigma}_{1}^{{12}}\rangle  =& -0.5 R \langle a^\dagger  {\sigma}_{1}^{{12}}\rangle  -0.5 \Gamma \langle a^\dagger  {\sigma}_{1}^{{12}}\rangle  -0.5 \kappa \langle a^\dagger  {\sigma}_{1}^{{12}}\rangle  -0.5 \nu \langle a^\dagger  {\sigma}_{1}^{{12}}\rangle  + 1 i g_{1} \langle {\sigma}_{1}^{{22}}\rangle  -1 i g_{1} \langle a^\dagger  a\rangle  -1 i \Delta \langle a^\dagger  {\sigma}_{1}^{{12}}\rangle  + 1 i g_{1} \left( -1 + N \right) \langle {\sigma}_{1}^{{12}}  {\sigma}_{2}^{{21}}\rangle  + 2 i g_{1} \langle {\sigma}_{1}^{{22}}\rangle  \langle a^\dagger  a\rangle  \\
\frac{d}{dt} \langle {\sigma}_{1}^{{12}}  {\sigma}_{2}^{{21}}\rangle  =& \left( -1.0 R -1.0 \Gamma \right) \langle {\sigma}_{1}^{{12}}  {\sigma}_{2}^{{21}}\rangle  + 1 i g_{1} \langle a^\dagger  {\sigma}_{1}^{{12}}\rangle  -1 i g_{1} \langle a  {\sigma}_{1}^{{21}}\rangle  -1.0 \nu \langle {\sigma}_{1}^{{12}}  {\sigma}_{2}^{{21}}\rangle  -2 i g_{1} \langle {\sigma}_{1}^{{22}}\rangle  \langle a^\dagger  {\sigma}_{1}^{{12}}\rangle  + 2 i g_{1} \langle {\sigma}_{1}^{{22}}\rangle  \langle a  {\sigma}_{1}^{{21}}\rangle
\end{align}
```



To calculate the dynamics of the system we create a system of ordinary differential equations, which can be used by [DifferentialEquations.jl](https://diffeq.sciml.ai/stable/).


```@example superradiant_laser_indexed
@named sys = ODESystem(eqs_sc)
nothing # hide
```

Finally we need to define the numerical parameters and the initial value of the system. We will consider $2 \cdot 10^5$ Strontium atoms which are repumped with a rate of $R = 1\text{Hz}$ on the clock transition ($\Gamma = 1 \text{mHz}$). The atom-cavity coupling rate is $g = 1\text{Hz}$, the cavity has a linewidth of $\kappa = 5\text{kHz}$ and is detuned from the atomic resonance by $\Delta = 2.5\text{Hz}$.


```@example superradiant_laser_indexed
# Initial state
u0 = zeros(ComplexF64, length(eqs_sc))
# System parameters
N_ = 2e5
Γ_ = 1.0 #Γ=1mHz
Δ_ = 2500Γ_ #Δ=2.5Hz
g_ = 1000Γ_ #g=1Hz
κ_ = 5e6*Γ_ #κ=5kHz
R_ = 1000Γ_ #R=1Hz
ν_ = 1000Γ_ #ν=1Hz

ps = [N, Δ, g(1), κ, Γ, R, ν]
p0 = [N_, Δ_, g_, κ_, Γ_, R_, ν_]

prob = ODEProblem(sys,u0,(0.0, 1.0/50Γ_), ps.=>p0)
nothing # hide
```


```@example superradiant_laser_indexed
# Solve the numeric problem
sol = solve(prob,Tsit5(),maxiters=1e7)

# Plot time evolution
t = sol.t
n = real.(sol[a'a])
s22 = real.(sol[σ(2,2,1)])
# Plot
p1 = plot(t, n, xlabel="tΓ", ylabel="⟨a⁺a⟩", legend=false)
p2 = plot(t, s22, xlabel="tΓ", ylabel="⟨σ22⟩", legend=false)
plot(p1, p2, layout=(1,2), size=(700,300))
```


## Spectrum

We calculate the spectrum here with the Laplace transform of the two-time correlation function. This is implemented with the function $\texttt{Spectrum}$.


```@example superradiant_laser_indexed
corr = CorrelationFunction(a', a, eqs_c; steady_state=true, filter_func=phase_invariant)
corr_sc = scale(corr)
S = Spectrum(corr_sc, ps)
nothing # hide
```

The set of equations for the correlation function is given by


```@example superradiant_laser_indexed
corr_sc.de
nothing # hide
```



```math
\begin{align}
\frac{d}{d\tau} \langle a^\dagger  a_0\rangle  =& -0.5 \kappa \langle a^\dagger  a_0\rangle  -1 i \Delta \langle a^\dagger  a_0\rangle  + 1 i N g_{1} \langle {\sigma}_{1}^{{21}}  a_0\rangle  \\
\frac{d}{d\tau} \langle {\sigma}_{1}^{{21}}  a_0\rangle  =& -0.5 R \langle {\sigma}_{1}^{{21}}  a_0\rangle  + 1 i g_{1} \langle a^\dagger  a_0\rangle  -0.5 \Gamma \langle {\sigma}_{1}^{{21}}  a_0\rangle  -0.5 \nu \langle {\sigma}_{1}^{{21}}  a_0\rangle  -2 i g_{1} \langle {\sigma}_{1}^{{22}}\rangle  \langle a^\dagger  a_0\rangle
\end{align}
```



To ensure we are in the steady state we use a steady solver to calculate it. To this end we need to define the $\texttt{SteadyStateProblem}$ and specify the desired method. We also need to increase the $\texttt{maxiters}$ and the solver accuracy to handle this numerically involved problem.


```@example superradiant_laser_indexed
prob_ss = SteadyStateProblem(sys,sol.u[end],ps.=>p0)
sol_ss = solve(prob_ss, SSRootfind())
nothing # hide
```

The spectrum is then calculated with


```@example superradiant_laser_indexed
ω = [-10:0.01:10;]Γ_
spec = S(ω,sol_ss.u,p0)
spec_n = spec ./ maximum(spec)
δ = abs(ω[(findmax(spec)[2])])
nothing # hide
```



```@example superradiant_laser_indexed
plot(ω, spec_n, xlabel="ω/Γ", legend=false, size=(500,300))
```


Beside the narrow linewidth we can also see another key feature of the superradiant laser here, namely the very weak cavity pulling. At a detunig of $\Delta = 2500\Gamma$ there is only a shift of the laser light from the atomic resonance frequency of $\delta = 1\Gamma$.
