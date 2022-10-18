# Superradiant Laser using Symbolic Summations

We can use the implemented indexing and summation features to calculate the Superradiant Laser example. Here we take advantage of these functionalities to simplify the definition of our Hamiltonian and the completion of the resulting equations of motion. The Hamiltonian of the system is once again given as:

$H = - \hbar \Delta a^\dagger a +  \hbar \sum\limits_{j=1}^{N}  g_j (a^\dagger \sigma^{12}_{j} + a \sigma^{21}_{j}),$

where $\Delta = \omega_a - \omega_c$ is the detuning between the cavity ($\omega_c$) and the atomic ($\omega_a$) resonance frequency, the atom cavity coupling of the atom $j$ is denoted by $g_j$. Additionally there are dissipative processes in the system, namely: Atoms are incoherently pumped with the rate $R$, they decay individually with the rate $\Gamma$ and are affected by individual atomic dephasing with the rate $\nu$. Photons also leak out of the system with the rate $\kappa$.

First of all we need to load the needed packages:


```@example superradiant_laser_indexed
using QuantumCumulants
using OrdinaryDiffEq, SteadyStateDiffEq, ModelingToolkit, DifferentialEquations
using Plots
nothing #hide
```

We continue by defining some of our parameters needed as $@cnumbers$, as well as the order of the cumulant expansion we want to have. Further more we define the Hilberspaces for our system, in this case our system consists of a **FockSpace** (the cavity) and a **NLevelSpace** (the atoms). It is important to note here, that we do not need to construct multiple **NLevelSpaces** for different atoms, since we can use the features of indices, which allow us to distinguish between different atoms in the calculation.


```@example superradiant_laser_indexed
order = 2 #order of the cumulant expansion
@cnumbers N Δ g κ Γ R ν

# Hilbertspace
hc = FockSpace(:cavity)
ha = NLevelSpace(:atom,2)

h = hc ⊗ ha
nothing #hide
```

    ℋ(cavity) ⊗ ℋ(atom)

In the next step we define those indices by using the constructor **Index**. Here we create several indices at once, since we want to use different indices for different atoms. However, each of the defined **Index** objects has the same range **N** and are defined on the same sub-Hilbertspace **ha**. We can further use these indices to define **IndexedOperator** objects, which are operators associated with an **Index**. Here we also create the Hamiltonian of the System using the symbolic Summation function **Σ**, which creates summation objects.


```@example superradiant_laser_indexed
# Indices and Operators
i = Index(h,:i,N,ha)
k = Index(h,:k,N,ha)
l = Index(h,:l,N,ha)

@qnumbers a::Destroy(h)
σ(i,j,k) = IndexedOperator(Transition(h,:σ,i,j),k)

# Define the Hamiltonian
H = -Δ*a'a + g*(Σ(a'*σ(1,2,i),i) + Σ(a*σ(2,1,i),i))
nothing #hide
```

```math
\underset{i}{\overset{N}{\sum}} g  a^\dagger  {\sigma}_{i}^{{12}} + \underset{i}{\overset{N}{\sum}} g  a  {\sigma}_{i}^{{21}} -1 \Delta a^\dagger a
```


In the next step we proceed by defining the dissipation operators with their corresponding dissipation rates and the operators, for which we want to calculate the equations of motion. In the final line of this code-block we calculate the meanfield eqautions using the **indexedMeanfield** function, which uses the features of symbolic summations.


```@example superradiant_laser_indexed
# Define Jump-Operators with corresponding rates
J = [a,σ(1,2,l),σ(2,1,l),σ(2,2,l)]
rates = [κ, Γ, R, ν]

# Define Operators, for which the meanfield shall be calculated
ops = [a'*a,σ(2,2,k)]

# It is best-practice to use every Index-Entity in only one context

# create Meanfield-Equations with given order for the given operators
eqs = indexed_meanfield(ops,H,J;rates=rates,order=order)
nothing #hide
```

```math
\begin{align}
\frac{d}{dt} \langle a^\dagger  a\rangle  =& 1 i \underset{i}{\overset{N}{\sum}} g  \langle a  {\sigma}_{i}^{{21}}\rangle  -1 i \underset{i}{\overset{N}{\sum}} g  \langle a^\dagger  {\sigma}_{i}^{{12}}\rangle  -1.0 \kappa \langle a^\dagger  a\rangle  \\
\frac{d}{dt} \langle {\sigma}_{k}^{{22}}\rangle  =& R -1.0 R \langle {\sigma}_{k}^{{22}}\rangle  + 1 i g \langle a^\dagger  {\sigma}_{k}^{{12}}\rangle  -1.0 \Gamma \langle {\sigma}_{k}^{{22}}\rangle  -1 i g \langle a  {\sigma}_{k}^{{21}}\rangle 
\end{align}
```



To get a closed set of equations we automatically complete the system. Since this system is phase invariant we know that all averages with a phase are zero, therefore we exclude these terms with a filter function.


```@example superradiant_laser_indexed
# custom filter function
# Using the filter function defined below, one can reduce the size of the complete system to only contain phase-invariant terms
φ(x::Average) = φ(x.arguments[1])
φ(::Destroy) = -1
φ(::Create) =1
φ(x::QTerm) = sum(map(φ, x.args_nc))
φ(x::Transition) = x.i - x.j
φ(x::IndexedOperator) = x.op.i - x.op.j
φ(x::IndexedSingleSum) = φ(x.term)
φ(x::AvgSums) = φ(arguments(x))
phase_invariant(x) = iszero(φ(x))


# We use the extraIndices keyword to provide names for indices, that are needed for intermediate calculation
eqs_c = complete(eqs;filter_func=phase_invariant,scaling=false,extra_indices=[:q])
nothing # hide
```

```math
\begin{align}
\frac{d}{dt} \langle a^\dagger  a\rangle  =& 1 i \underset{i}{\overset{N}{\sum}} g  \langle a  {\sigma}_{i}^{{21}}\rangle  -1 i \underset{i}{\overset{N}{\sum}} g  \langle a^\dagger  {\sigma}_{i}^{{12}}\rangle  -1.0 \kappa \langle a^\dagger  a\rangle  \\
\frac{d}{dt} \langle {\sigma}_{k}^{{22}}\rangle  =& R -1.0 R \langle {\sigma}_{k}^{{22}}\rangle  + 1 i g \langle a^\dagger  {\sigma}_{k}^{{12}}\rangle  -1.0 \Gamma \langle {\sigma}_{k}^{{22}}\rangle  -1 i g \langle a  {\sigma}_{k}^{{21}}\rangle  \\
\frac{d}{dt} \langle a^\dagger  {\sigma}_{k}^{{12}}\rangle  =& 1 i \underset{i{\ne}k}{\overset{N}{\sum}} g  \langle {\sigma}_{i}^{{21}}  {\sigma}_{k}^{{12}}\rangle  + 1 i g \langle {\sigma}_{k}^{{22}}\rangle  -1 i g \langle a^\dagger  a\rangle  -0.5 R \langle a^\dagger  {\sigma}_{k}^{{12}}\rangle  -0.5 \Gamma \langle a^\dagger  {\sigma}_{k}^{{12}}\rangle  -0.5 \kappa \langle a^\dagger  {\sigma}_{k}^{{12}}\rangle  -0.5 \nu \langle a^\dagger  {\sigma}_{k}^{{12}}\rangle  -1 i \Delta \langle a^\dagger  {\sigma}_{k}^{{12}}\rangle  + 2 i g \langle {\sigma}_{k}^{{22}}\rangle  \langle a^\dagger  a\rangle  \\
\frac{d}{dt} \langle {\sigma}_{k}^{{12}}  {\sigma}_{q}^{{21}}\rangle  =& -1.0 R \langle {\sigma}_{k}^{{12}}  {\sigma}_{q}^{{21}}\rangle  -1.0 \Gamma \langle {\sigma}_{k}^{{12}}  {\sigma}_{q}^{{21}}\rangle  -1.0 \nu \langle {\sigma}_{k}^{{12}}  {\sigma}_{q}^{{21}}\rangle  + 1.0 i g \langle a^\dagger  {\sigma}_{k}^{{12}}\rangle  -1.0 i g \langle a  {\sigma}_{q}^{{21}}\rangle  -2.0 i g \langle {\sigma}_{q}^{{22}}\rangle  \langle a^\dagger  {\sigma}_{k}^{{12}}\rangle  + 2.0 i g \langle {\sigma}_{k}^{{22}}\rangle  \langle a  {\sigma}_{q}^{{21}}\rangle 
\end{align}
```



The equations, that we got in the previous step can now be easily adjusted for different purposes. We assume here all atoms in the system behave identically. We can then simply reduce the above equations to equations for specific atoms (atom 1,2,3...) by using the **scale** function.


```@example superradiant_laser_indexed
# Now one can easily scale the Equations above
eqs_sc = scale(eqs_c)
nothing #hide
```

```math
\begin{align}
\frac{d}{dt} \langle a^\dagger  a\rangle  =& -1.0 \kappa \langle a^\dagger  a\rangle  -1 i N g \langle a^\dagger  {\sigma}_{1}^{{12}}\rangle  + 1 i N g \langle a^\dagger  {\sigma}_{1}^{{12}}\rangle ^* \\
\frac{d}{dt} \langle {\sigma}_{1}^{{22}}\rangle  =& R -1.0 R \langle {\sigma}_{1}^{{22}}\rangle  + 1 i g \langle a^\dagger  {\sigma}_{1}^{{12}}\rangle  -1 i g \langle a^\dagger  {\sigma}_{1}^{{12}}\rangle ^* -1.0 \Gamma \langle {\sigma}_{1}^{{22}}\rangle  \\
\frac{d}{dt} \langle a^\dagger  {\sigma}_{1}^{{12}}\rangle  =& 1 i g \langle {\sigma}_{1}^{{22}}\rangle  -1 i g \langle a^\dagger  a\rangle  -0.5 R \langle a^\dagger  {\sigma}_{1}^{{12}}\rangle  -0.5 \Gamma \langle a^\dagger  {\sigma}_{1}^{{12}}\rangle  -0.5 \kappa \langle a^\dagger  {\sigma}_{1}^{{12}}\rangle  -1 i \Delta \langle a^\dagger  {\sigma}_{1}^{{12}}\rangle  -0.5 \nu \langle a^\dagger  {\sigma}_{1}^{{12}}\rangle  + 2 i g \langle {\sigma}_{1}^{{22}}\rangle  \langle a^\dagger  a\rangle  + 1 i g \left( -1 + N \right) \langle {\sigma}_{1}^{{12}}  {\sigma}_{2}^{{21}}\rangle  \\
\frac{d}{dt} \langle {\sigma}_{1}^{{12}}  {\sigma}_{2}^{{21}}\rangle  =& -1.0 R \langle {\sigma}_{1}^{{12}}  {\sigma}_{2}^{{21}}\rangle  -1.0 \Gamma \langle {\sigma}_{1}^{{12}}  {\sigma}_{2}^{{21}}\rangle  -1.0 \nu \langle {\sigma}_{1}^{{12}}  {\sigma}_{2}^{{21}}\rangle  + 1.0 i g \langle a^\dagger  {\sigma}_{1}^{{12}}\rangle  -1.0 i g \langle a^\dagger  {\sigma}_{1}^{{12}}\rangle ^* -2.0 i g \langle {\sigma}_{1}^{{22}}\rangle  \langle a^\dagger  {\sigma}_{1}^{{12}}\rangle  + 2.0 i g \langle {\sigma}_{1}^{{22}}\rangle  \langle a^\dagger  {\sigma}_{1}^{{12}}\rangle ^*
\end{align}
```



To calculate the dynamics of the system we create a system of ordinary differential equations, which can be used by [DifferentialEquations.jl](https://diffeq.sciml.ai/stable/). Furthermore we define the numerical parameters as well as the initial state of the system and solve the system in the next step.


```@example superradiant_laser_indexed
# define the ODE System and Problem
@named sys = ODESystem(eqs_sc)

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

ps = [N, Δ, g, κ, Γ, R, ν]
p0 = [N_, Δ_, g_, κ_, Γ_, R_, ν_]

prob = ODEProblem(sys,u0,(0.0, 1.0/50Γ_), ps.=>p0)
nothing #hide
```


```@example superradiant_laser_indexed
# Solve the Problem
sol = solve(prob,maxiters=1e7)

# Plot time evolution
t = sol.t
n = real.(sol[a'a])
s22 = real.(sol[σ(2,2,1)])
# Plot
p1 = plot(t, n, xlabel="tΓ", ylabel="⟨a⁺a⟩", legend=false)
p2 = plot(t, s22, xlabel="tΓ", ylabel="⟨σ22⟩", legend=false)
plot(p1, p2, layout=(1,2), size=(700,300))
savefig("superradiant_laser_indexed.svg") # hide
```

![svg](superradiant_laser_indexed.svg)    

# Spectrum

We can now calculate the spectrum using the Laplace transform of the two-ime correlation function. This is done here with the **Spectrum** function.


```@example superradiant_laser_indexed
# For the Spectrum 
# setting the scaling keyword attribute to true gives us again a output in similar form to scale(eqs_c)
corr = CorrelationFunction(a', a, eqs_c; steady_state=true, filter_func=phase_invariant,scaling=true);
S = Spectrum(corr, ps)
nothing #hide
```

The set of equations for the correlation function is then given by:


```@example superradiant_laser_indexed
corr.de
nothing #hide
```

```math
\begin{align}
\frac{d}{d\tau} \langle a^\dagger  a_0\rangle  =& -0.5 \kappa \langle a^\dagger  a_0\rangle  -1 i \Delta \langle a^\dagger  a_0\rangle  + 1 i N g \langle a_0  {\sigma}_{1}^{{21}}\rangle  \\
\frac{d}{d\tau} \langle a_0  {\sigma}_{1}^{{21}}\rangle  =& -0.5 R \langle a_0  {\sigma}_{1}^{{21}}\rangle  -0.5 \Gamma \langle a_0  {\sigma}_{1}^{{21}}\rangle  -0.5 \nu \langle a_0  {\sigma}_{1}^{{21}}\rangle  + 1 i g \langle a^\dagger  a_0\rangle  -2 i g \langle {\sigma}_{1}^{{22}}\rangle  \langle a^\dagger  a_0\rangle 
\end{align}
```

To ensure we are in the steady state we use a steady solver to calculate it. To this end we need to define the SteadyStateProblem and specify the desired method. We also need to increase the maxiters and the solver accuracy to handle this numerically involved problem.


```@example superradiant_laser_indexed
prob_ss = SteadyStateProblem(prob)
sol_ss = solve(prob_ss, DynamicSS(Tsit5(); abstol=1e-8, reltol=1e-8),
    reltol=1e-14, abstol=1e-14, maxiters=5e7)
nothing #hide
```


```@example superradiant_laser_indexed
ω = [-10:0.01:10;]Γ_
spec = S(ω,sol_ss.u,p0)
spec_n = spec ./ maximum(spec)
δ = abs(ω[(findmax(spec)[2])])
nothing #hide
```


```@example superradiant_laser_indexed
plot(ω, spec_n, xlabel="ω/Γ", legend=false, size=(500,300))
savefig("superradiant_laser_indexed_spectrum.svg") # hide
```


    
![svg](superradiant_laser_indexed_spectrum.svg)
