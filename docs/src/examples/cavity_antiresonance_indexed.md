# Cavity Antiresonance using Symbolic Summations

In this example we investigate a system of $N$ closely spaced quantum emitters inside a coherently driven single mode cavity. The model is descriped in [D. Plankensteiner, et. al., Phys. Rev. Lett. 119, 093601 (2017)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.119.093601).
The Hamiltonian of this system is composed of three parts $H = H_c + H_a + H_{\mathrm{int}}$, the driven cavity $H_c$, the dipole-dipole interacting atoms $H_a$ and the atom-cavity interaction $H_\mathrm{int}$:

```math
\begin{align*}
H_\mathrm{c} = \Delta_c a^\dagger a + \eta (a^\dagger + a) \\
\\
H_a = \Delta_a \sum\limits_{i} \sigma_i^{22} + \sum\limits_{i \neq j} \Omega_{ij} \sigma_i^{21} \sigma_j^{12} \\
\\
H_\mathrm{int} = \sum\limits_{i} g_i (a^\dagger \sigma_i^{12} + a \sigma_i^{21})
\end{align*}
```

Additionally the system features two decay channels, the lossy cavity with photon decay rate $\kappa$ and collective atomic emission described by the decay-rate matrix $\Gamma_{ij}$. In our example we will consider the case of $N=2$ atoms.

We start by loading the packages.



```@example cavity_antiresonance_indexed
using QuantumCumulants
using OrdinaryDiffEq, SteadyStateDiffEq, ModelingToolkit
using Plots
nothing #hide
```

The Hilbert space for this system is given by one cavity mode and $N$ two-level atoms and the parameters $g_j, \, \Gamma_{ij}$ and $\Omega_{ij}$ are defined as **IndexedVariable** of atom $i$ and $j$. We will describe the system in first order mean-field.


```@example cavity_antiresonance_indexed
order = 1
@cnumbers Δc η Δa κ N

hc = FockSpace(:cavity)
ha = NLevelSpace(Symbol(:atom),2)
h = hc ⊗ ha

#define indices
i = Index(h,:i,N,ha)
j = Index(h,:j,N,ha)
k = Index(h,:k,N,ha)

N_n = 2 #number of atoms

#define indexed variables
gi = IndexedVariable(:g,i)
Γ_ij = IndexedVariable(:Γ,i,j)  
Ω_ij = IndexedVariable(:Ω,i,j;identical=false) #false indicates that the 2 indices can never be the same
nothing #hide
```

Now we create the operators on the composite Hilbert space using a **IndexedOperator** constructor to assign each **Transition** operator an index $k$.


```@example cavity_antiresonance_indexed
@qnumbers a::Destroy(h)
σ(x,y,k) = IndexedOperator(Transition(h,:σ,x,y),k)
nothing #hide
```


We define the Hamiltonian and Liouvillian. For the collective atomic decay we can write the corresponding jump process with a rate-matrix $R$ and a vector $J$ of jump operators, such that an operator $\mathcal{O}$ follows the equation

$\dot{\mathcal{O}} = \sum_{ij} R_{ij} (J_i^\dagger \mathcal{O} J_j - J_i^\dagger J_j \mathcal{O}/2 -  \mathcal{O} J_i^\dagger J_j/2) + \mathrm{noise}.$


```@example cavity_antiresonance_indexed
# Hamiltonian
Hc = Δc*a'a + η*(a' + a)

Ha = Δa*Σ(σ(2,2,i),i) + Σ(Ω_ij*σ(2,1,i)*σ(1,2,j),j,i)

Hi = Σ(gi*(a'*σ(1,2,i) + a*σ(2,1,i)),i)
H = Hc + Ha + Hi

# Jump operators & and rates
J = [a, [σ(1,2,i),σ(1,2,j)] ] 
rates = [κ,Γ_ij]
nothing #hide
```

Using the Hamiltonian and Liouvillian we derive the system of equations in first order mean-field.



```@example cavity_antiresonance_indexed
ops = [a, σ(2,2,k), σ(1,2,k)]
eqs = indexed_meanfield(ops,H,J;rates=rates,order=order)
nothing #hide
```

```math
\begin{align}
\frac{d}{dt} \langle a\rangle  =& -1 i \eta -1 i \underset{i}{\overset{N}{\sum}} {g}_{i}  \langle {\sigma}_{i}^{{12}}\rangle  -1 i {\Delta}c \langle a\rangle  -0.5 \kappa \langle a\rangle  \\
\frac{d}{dt} \langle {\sigma}_{k}^{{22}}\rangle  =& -0.5 \underset{i{\ne}j,k}{\overset{N}{\sum}} {\Gamma}_{ik}  \langle {\sigma}_{i}^{{21}}\rangle   \langle {\sigma}_{k}^{{12}}\rangle  -0.5 \underset{j{\ne}k}{\overset{N}{\sum}} {\Gamma}_{kj}  \langle {\sigma}_{k}^{{21}}\rangle   \langle {\sigma}_{j}^{{12}}\rangle  + 1 i \underset{i{\ne}j,k}{\overset{N}{\sum}} {\Omega}_{ik}  \langle {\sigma}_{i}^{{21}}\rangle   \langle {\sigma}_{k}^{{12}}\rangle  -1 i \underset{j{\ne}i,k}{\overset{N}{\sum}} {\Omega}_{kj}  \langle {\sigma}_{k}^{{21}}\rangle   \langle {\sigma}_{j}^{{12}}\rangle  -1.0 {\Gamma}_{kk} \langle {\sigma}_{k}^{{22}}\rangle  -1 i {g}_{k} \langle a\rangle  \langle {\sigma}_{k}^{{21}}\rangle  + 1 i {g}_{k} \langle a^\dagger\rangle  \langle {\sigma}_{k}^{{12}}\rangle  \\
\frac{d}{dt} \langle {\sigma}_{k}^{{12}}\rangle  =& \underset{j{\ne}k}{\overset{N}{\sum}} {\Gamma}_{kj}  \langle {\sigma}_{j}^{{12}}\rangle   \langle {\sigma}_{k}^{{22}}\rangle  -1 i \underset{j{\ne}i,k}{\overset{N}{\sum}} {\Omega}_{kj}  \langle {\sigma}_{j}^{{12}}\rangle  -0.5 \underset{j}{\overset{N}{\sum}} {\Gamma}_{kj}  \langle {\sigma}_{j}^{{12}}\rangle  + 2 i \underset{j{\ne}i,k}{\overset{N}{\sum}} {\Omega}_{kj}  \langle {\sigma}_{j}^{{12}}\rangle   \langle {\sigma}_{k}^{{22}}\rangle  -1 i {g}_{k} \langle a\rangle  -1 i {\Delta}a \langle {\sigma}_{k}^{{12}}\rangle  + 2 i {g}_{k} \langle a\rangle  \langle {\sigma}_{k}^{{22}}\rangle 
\end{align}
```

We complete the set of equations automatically.

```@example cavity_antiresonance_indexed
eqs_comp = complete(eqs;extra_indices=[:q])
nothing #hide
```

We now evaluate the symbolic indices inside these equations to their numerical value by calling **evaluate** and create an ordinary differential equation system in order to solve it numerically.


```@example cavity_antiresonance_indexed
eqs_ = evaluate(eqs_comp;mapping=(N=>N_n)) #use the numerical value of N in mapping keyword
@named sys = ODESystem(eqs_)
nothing #hide
```

Finally we need to define the initial state of the system and the numerical parameters. In the end we want to obtain the transmission rate $T$ of our system. For this purpose we calculate the steady state photon number in the cavity $|\langle a \rangle|^2$ for different laser frequencies.


```@example cavity_antiresonance_indexed
u0 = zeros(ComplexF64, length(eqs_))
# parameter
Γ_ = 1.0
d = 2π*0.08 #0.08λ
θ = π/2

Ωij_(i,j) = Γ_*(-3/4)*( (1-(cos(θ))^2)*cos(d)/d-(1-3*(cos(θ))^2)*(sin(d)/(d^2)+(cos(d)/(d^3))) )
function Γij_(i,j)
    i==j ? Γ_ : Γ_*(3/2)*( (1-(cos(θ))^2)*sin(d)/d+(1-3*(cos(θ))^2)*((cos(d)/(d^2))-(sin(d)/(d^3))))
end

ΓMatrix = [Γij_(i,j) for i = 1:N_n, j=1:N_n]
ΩMatrix = [Ωij_(i,j) for i = 1:N_n, j=1:N_n]

g_ = 2Γ_
κ_ = 20Γ_
Δa_ = 0Γ_
Δc_ = 0Γ_
η_ = κ_/100

g_v = [g_*(-1)^j for j=1:N_n]
ps = [Δc, η, Δa, κ, gi, Γ_ij, Ω_ij]
nothing #hide
```

The transmission rate $T$ with respect to the pump laser detuning is given by the relative steady state intra-cavity photon number $n(\Delta)/n_\mathrm{max}$. We qualitatively reproduce the antiresonance from [D. Plankensteiner, et. al., Phys. Rev. Lett. 119, 093601 (2017)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.119.093601) for two atoms.


```@example cavity_antiresonance_indexed
Δ_ls = [-10:0.05:10;]Γ_
n_ls = zeros(length(Δ_ls))

for i=1:length(Δ_ls)
    Δc_i = Δ_ls[i]
    Δa_i = Δc_i + Ωij_(1,2) #cavity on resonace with the shifted collective emitter
    p0_i = [Δc_i, η_, Δa_i, κ_, g_v, ΓMatrix, ΩMatrix]
    ps_ = value_map(ps,p0_i;mapping=(N=>N_n)) #Combine all the parameters + values to one list for solving
    prob = ODEProblem(sys,u0,(0.0, 20Γ_), ps_);
    prob_ss = SteadyStateProblem(prob);
    sol_ss = solve(prob_ss, DynamicSS(Tsit5(); abstol=1e-8, reltol=1e-8),
        reltol=1e-14, abstol=1e-14, maxiters=5e7)
    n_ls[i] = abs2(sol_ss[a])
end
nothing #hide
```


```@example cavity_antiresonance_indexed
T = n_ls ./ maximum(n_ls)
p = plot(Δ_ls, T, xlabel="Δ/Γ", ylabel="T", legend=false)
savefig("cavity_antiresonance_indexed.svg") # hide
```

![svg](cavity_antiresonance_indexed.svg)
