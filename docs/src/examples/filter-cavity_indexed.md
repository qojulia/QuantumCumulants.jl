# Laser with Filter Cavities

An intuitive and straightforward approach to calculate the spectrum of a laser is to filter the emitted light. We can do this by coupling filter cavities with different detunings to the main cavity and observe the photon number in these 'filters', see for example [K. Debnath et al., Phys Rev A 98, 063837 (2018)](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.98.063837).

The main goal of this example is to combine two indexed Hilbert spaces, where one will be scaled and the other evaluated. The model is basically the same as for the [superradiant laser](https://qojulia.github.io/QuantumCumulants.jl/stable/examples/superradiant_laser_indexed/) example, but with the additional terms due to the filter cavities. The Hamiltonian of this system is

```math
\begin{equation}
H = - \Delta a^\dagger a + g \sum\limits_{j=1}^{N} (a^\dagger \sigma^{12}_{j} + a \sigma^{21}_{j}) - \sum\limits_{i=1}^{M} \delta_i b_i^\dagger b_i +  g_f \sum\limits_{i=1}^{M} (a^\dagger b_i + a b_i^\dagger),
\end{equation}
```

where $\delta_i$ is the detuning of the $i$-th filter cavity and $g_f$ the coupling to the normal cavity. Furthermore, their decay rate is $\kappa_f$.

We start by loading the packages.


```@example filter_cavity_indexed
using QuantumCumulants
using OrdinaryDiffEq, ModelingToolkit
using Plots
```

We create the parameters of the system including the $\texttt{IndexedVariable}$ $\delta_i$. For the atoms and filter cavities we only need one Hilbert space each. We define the indices for each Hilbert space and use them to create $\texttt{IndexedOperators}$.


```@example filter_cavity_indexed
# Parameters
@cnumbers κ g gf κf R Γ Δ ν N M
δ(i) = IndexedVariable(:δ, i)

# Hilbertspace
hc = FockSpace(:cavity)
hf = FockSpace(:filter)
ha = NLevelSpace(:atom, 2)
h = hc ⊗ hf ⊗ ha

# Indices and Operators
i = Index(h,:i,M,hf)
j = Index(h,:j,N,ha)

@qnumbers a::Destroy(h,1)
b(k) = IndexedOperator(Destroy(h,:b,2), k)
σ(α,β,k) = IndexedOperator(Transition(h,:σ,α,β,3), k)
nothing # hide
```

We define the Hamiltonian using symbolic sums and define the individual dissipative processes. For an indexed jump operator the (symbolic) sum is build in the Liouvillian, in this case corresponding to individual decay processes.


```@example filter_cavity_indexed
# Hamiltonian
H = Δ*Σ(σ(2,2,j),j) + Σ(δ(i)*b(i)'b(i),i) +
    gf*(Σ(a'*b(i) + a*b(i)',i)) + g*(Σ(a'*σ(1,2,j) + a*σ(2,1,j),j))

# Jumps & rates
J = [a, b(i), σ(1,2,j), σ(2,1,j), σ(2,2,j)]
rates = [κ, κf, Γ, R, ν]
nothing # hide
```

We derive the equation for $\langle a^\dagger a \rangle$ and complete the system automatically in second order.


```@example filter_cavity_indexed
eqs = meanfield(a'a,H,J;rates=rates,order=2)
nothing # hide
```




```math
\begin{align}
\frac{d}{dt} \langle a^\dagger  a\rangle  =& 1 i \underset{i}{\overset{M}{\sum}} gf  \langle a  {b}_{i}^\dagger\rangle  -1 i \underset{i}{\overset{M}{\sum}} gf  \langle a^\dagger  {b}_{i}\rangle  + 1 i \underset{j}{\overset{N}{\sum}} g  \langle a  {\sigma}_{j}^{{21}}\rangle  -1 i \underset{j}{\overset{N}{\sum}} g  \langle a^\dagger  {\sigma}_{j}^{{12}}\rangle  -1.0 \kappa \langle a^\dagger  a\rangle
\end{align}
```





```@example filter_cavity_indexed
eqs_c = complete(eqs);
nothing # hide
```

Now we assume that all atoms behave identically, but we want to obtain the equations for 20 different filter cavities. To this end we $\texttt{scale}$ the Hilbert space of the atoms and $\texttt{evaluate}$ the filter cavities. Specifying the Hilbert space is done with the kwarg $\texttt{h}$, which can either be the specific Hilbert space or its acts-on number. Evaluating the filter cavities requires a numeric upper bound for the used $\texttt{Index}$, we provide this with a dictionary on the kwarg $\texttt{limits}$.


```@example filter_cavity_indexed
M_ = 20
eqs_sc = scale(eqs_c;h=[ha]) #h=[3]
eqs_eval = evaluate(eqs_sc; limits=Dict(M=>M_)) #h=[hf]
println("Number of eqs.: $(length(eqs_eval))")
```

To calculate the dynamics of the system we create a system of ordinary differential equations, which can be used by [DifferentialEquations.jl](https://diffeq.sciml.ai/stable/). Finally we need to define the numerical parameters and the initial state of the system.


```@example filter_cavity_indexed
@named sys = ODESystem(eqs_eval)
nothing # hide
```


```@example filter_cavity_indexed
# Initial state
u0 = zeros(ComplexF64, length(eqs_eval))

# System parameters
N_ = 200
Γ_ = 1.0
Δ_ = 0Γ_
g_ = 1Γ_
κ_ = 100Γ_
R_ = 10Γ_
ν_ = 1Γ_

gf_ = 0.1Γ_
κf_ = 0.1Γ_
δ_ls = [0:1/M_:1-1/M_;]*10Γ_

ps = [Γ, κ, g, κf, gf, R, [δ(i) for i=1:M_]..., Δ, ν, N]
p0 = [Γ_, κ_, g_, κf_, gf_, R_, δ_ls..., Δ_, ν_, N_]

prob = ODEProblem(sys,u0,(0.0, 10.0/κf_), ps.=>p0)
nothing # hide
```


```@example filter_cavity_indexed
# Solve the numeric problem
sol = solve(prob, Tsit5(); abstol=1e-10, reltol=1e-10, maxiters=1e7)

t = sol.t
n = abs.(sol[a'a])
n_b(i) =  abs.(sol[b(i)'b(i)])
n_f = [abs(sol[b(i)'b(i)][end]) for i=1:M_] ./ (abs(sol[b(1)'b(1)][end]))
nothing # hide
```


```@example filter_cavity_indexed
# Plot results
p1 = plot(t, n_b(1), alpha=0.5, ylabel="⟨bᵢ⁺bᵢ⟩", legend=false)
for i=2:M_
    plot!(t, n_b(i), alpha=0.5, legend=false)
end
#p1 = plot!(twinx(), t, n, xlabel="tΓ", ylabel="⟨a⁺a⟩", legend=false)

p2 = plot([-reverse(δ_ls);δ_ls], [reverse(n_f);n_f], xlabel="δ/Γ", ylabel="intensity", legend=false)
plot(p1, p2, layout=(1,2), size=(700,300))
```
