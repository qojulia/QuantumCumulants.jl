```@meta
EditURL = "../../../examples/waveguide.jl"
```

# Waveguide Energy Transfer

In this example, we investigate the collective behaviour of atomic ensembles coupled to a waveguide. The waveguide-mediated dipole-dipole Hamiltonian is
```math
    H = \sum_{i \neq j}^{N} \Omega_{ij}^{1\mathrm{D}} \sigma^+_i \sigma^-_j,
```
with the coupling rate $\Omega^{1\mathrm{D}}_{ij}$.
The incoherent part describing the dissipative processes is accounted for by the Lindblad term
```math
    \mathcal{L}[\rho]=\frac{1}{2}\sum_{i,j}^N \Gamma^{1\mathrm{D}}_{ij} (2 \sigma^-_i \rho \sigma^+_j-\sigma^+_i \sigma^-_j \rho - \rho \sigma^+_i \sigma^-_j)
```
with the collective fiber-mediated decay rates $\Gamma^{1\mathrm{D}}_{ij}$. For $M$ atomic ensembles, where each of the $N_i$ atoms within an ensemble behaves identically, we can use $M$ spin-$N_i/2$ systems to describe them. This means
```math
    H = \sum_{i \neq j}^{M} {\Omega}_{ij}^{1\mathrm{D}} S^+_i S^-_j \hspace{0.75cm} \mathrm{and} \hspace{0.75cm} \mathcal{L}[\rho]=\frac{1}{2}\sum_{i,j}^M {\Gamma}^{1\mathrm{D}}_{ij} (2 S^-_i \rho S^+_j-S^+_i S^-_j \rho - \rho S^+_i S^-_j),
```
with $S^{\pm}_i = S^x \pm i S^y = \sum_{k=1}^{N_i} \sigma^{\pm}_k$, where $S^{x,y,z}_i = \frac{1}{2} \sum_{k=1}^{N_i} \sigma^{x,y,z}_k $.

After loading the needed packages we define the Hilbert space and spin operators of the $M$ atomic ensembles coupled to the waveguide. We also define the numerical Parameters and utilize that $\Omega_{ij} = \Omega_{ji}$ and $\Gamma_{ij} = \Gamma_{ji}$. Due to this the derived equations are simplified further.

````@example waveguide
using QuantumCumulants
using ModelingToolkit
using OrdinaryDiffEq
using QuantumOptics
using Plots

M_p = 2 # number of pumped spins
M_np = 2 # number of non-pumped spins
M = M_p + M_np

h_spin(i) = SpinSpace("spin_$(i)") # Hilbert space
h = tensor([h_spin(i) for i = 1:M]...)

Sx(i) = Spin(h, "S$(i)", 1, i) # Operators

Sy(i) = Spin(h, "S$(i)", 2, i)
Sz(i) = Spin(h, "S$(i)", 3, i)

Sm(i) = (Sx(i) - 1im*Sy(i))
Sp(i) = (Sx(i) + 1im*Sy(i))

Ω(i, j) = i>j ? cnumber("Ω_$(j)_$(i)") : cnumber("Ω_$(i)_$(j)") # Parameter
Γ(i, j) = i>j ? cnumber("Γ_$(j)_$(i)") : cnumber("Γ_$(i)_$(j)")
nothing # hide
````

With the symbolic operators and Parameters we define the Hamiltonian and the collective jump operators with the corresponding rates. For rates written in matrix form the program automatically assumes collective dissipation according to
```math
    \mathcal{L}[\rho]=\frac{1}{2}\sum_{i,j}^M R_{ij} (2 J_i \rho J^+_j-J^+_i \rho J_j - \rho J^+_i J_j),
```
with the jump operator $J_i$ and the corresponding decay rate matrix $R_{ij}$. This implementation is similar to the one in [QuantumOptics.jl](https://docs.qojulia.org/timeevolution/master/) for collective dissipation.

````@example waveguide
H = sum((i≠j)*Ω(i, j)*Sp(i)*Sm(j) for i = 1:M for j = 1:M) # Hamiltonian

J = [Sm(c1) for c1 = 1:M] # Jump operators and rate matrix
rates = [Γ(c1, c2) for c1 = 1:M, c2 = 1:M]
nothing # hide
````

Now we create a complete list of operators, which we use to derive the second-order equations. This is considerable faster than to automatically complete the equations. We only show one of the 90 rather lengthy equations.

````@example waveguide
S(i) = [Sx(i), Sy(i), Sz(i)] # create list of operators
SiSi(i) = [Sx(i)Sx(i), Sx(i)Sy(i), Sx(i)Sz(i), Sy(i)Sy(i), Sy(i)Sz(i), Sz(i)Sz(i)]
ops = [];
[push!(ops, S(i)...) for i = 1:M]
[push!(ops, SiSi(i)...) for i = 1:M]
for i = 1:M, j = i:M
    if i≠j
        for α = 1:3, β = 1:3
            push!(ops, S(i)[α]*S(j)[β])
        end
    end
end
println("Number of equations = $(length(ops))")

eqs = meanfield(ops, H, J; rates = rates, order = 2) # derive equations
meanfield(Sz(1), H, J; rates = rates, order = 2)
nothing # hide
````

```math
\begin{align}
\frac{d}{dt} \langle {S1}^{{z}}\rangle  =& -1.0 \langle {S1}^{{y}}  {S4}^{{y}}\rangle  \Gamma_{1\_4} -1.0 \langle {S1}^{{x}}  {S3}^{{x}}\rangle  \Gamma_{1\_3} -1.0 \langle {S1}^{{y}}  {S1}^{{y}}\rangle  \Gamma_{1\_1} -1.0 \langle {S1}^{{y}}  {S2}^{{y}}\rangle  \Gamma_{1\_2} -1.0 \langle {S1}^{{z}}\rangle  \Gamma_{1\_1} -1.0 \langle {S1}^{{x}}  {S1}^{{x}}\rangle  \Gamma_{1\_1} + 2 \langle {S1}^{{y}}  {S2}^{{x}}\rangle  \Omega_{1\_2} + 2 \langle {S1}^{{y}}  {S4}^{{x}}\rangle  \Omega_{1\_4} -1.0 \langle {S1}^{{x}}  {S2}^{{x}}\rangle  \Gamma_{1\_2} -1.0 \langle {S1}^{{x}}  {S4}^{{x}}\rangle  \Gamma_{1\_4} -2 \langle {S1}^{{x}}  {S3}^{{y}}\rangle  \Omega_{1\_3} -2 \langle {S1}^{{x}}  {S4}^{{y}}\rangle  \Omega_{1\_4} -1.0 \langle {S1}^{{y}}  {S3}^{{y}}\rangle  \Gamma_{1\_3} + 2 \langle {S1}^{{y}}  {S3}^{{x}}\rangle  \Omega_{1\_3} -2 \langle {S1}^{{x}}  {S2}^{{y}}\rangle  \Omega_{1\_2}
\end{align}
```

We define the numerical system Parameters and the initial state. For spin operators the initial state is unfortunately not fully trivial with all zeros. Therefore, we use the [numeric conversion](https://qojulia.github.io/QuantumCumulants.jl/stable/implementation/#Numeric-averages-and-conversion) which allows us to define the quantum state in [QuantumOptics.jl](https://docs.qojulia.org/quantumobjects/states/) and convert it to the correct initial average values.

````@example waveguide
Γ_ = 1.0  # Γ1D # System Parameters
Ω_ = Γ_/2 # Ω1D

x_p = [i/M_p for i = 0:(M_p-1)]*2π # positions
x_np = [i/M_np for i = 0:(M_np-1)]*2π
x_ls = [x_p; x_np]

Γij(i, j) = cos(abs(x_ls[i] - x_ls[j]))
Γ_ls = [Γij(i, j) for i = 1:M for j = 1:M]*Γ_
Ωij(i, j) = sin(abs(x_ls[i] - x_ls[j]))
Ω_ls = [Ωij(i, j) for i = 1:M for j = 1:M]*Ω_

ps = [[Γ(i, j) for i = 1:M for j = 1:M]; [Ω(i, j) for i = 1:M for j = 1:M];]
p0 = [Γ_ls; Ω_ls;]
nothing # hide
````

Note that the function [coherentspinstate](https://docs.qojulia.org/api/#QuantumOpticsBase.coherentspinstate) can only reliably produce states for spins with a length of about $10^4$, for [spindown](https://docs.qojulia.org/api/#QuantumOpticsBase.spindown) or [spinup](https://docs.qojulia.org/api/#QuantumOpticsBase.spinup) there is no restriction. Also, we created the function `LazyKet`, which allows you to define product states without calculating the tensor product directly. This makes it possible to define initial (product) states for very large quantum systems. The function [`initial_values`](@ref)(eqs, $\psi$) calculates the initial values for a set of equations `eqs` with the initial quantum state $\psi$.

````@example waveguide
N_p = 8000 # number of pumped atoms
N_np = 2000 # number of non-pumped atoms
bs1 = SpinBasis(N_p/2/M_p) # spin basis for pumped atoms
bs2 = SpinBasis(N_np/2/M_np) # spin basis for non-pumped atoms
b = tensor([bs1 for i = 1:M_p]..., [bs2 for i = 1:M_np]...)

Θ = π/3 # corresponds to an excitation of 75%
ψ_p = coherentspinstate(bs1, Θ, 0.0) # initially pumped atoms
ψ_np = spindown(bs2) # initially non-pumped atoms
ψ0 = LazyKet(b, ([ψ_p for i = 1:M_p]..., [ψ_np for i = 1:M_np]...))
u0 = initial_values(eqs, ψ0)  # initial state
nothing # hide
````

Finally, we create the ODE-problem, calculate the dynamics and plot the results. We see a fast transfer of the population from the initially excited to the other atomic ensemble.

````@example waveguide
@named sys = System(eqs)
dict = merge(Dict(unknowns(sys) .=> u0), Dict(ps .=> p0))
prob = ODEProblem(sys, dict, (0.0, 8e-3))
sol = solve(prob, Tsit5(), abstol = 1e-6, reltol = 1e-6)
nothing # hide
````

````@example waveguide
t = sol.t * 1e3
Sz_t(i) = real.(sol[Sz(i)])
N_ls = [[N_p for i = 1:M_p] ./ M_p; [N_np for i = 1:M_np] ./ M_np]
σ22_t(i) = Sz_t(i)/N_ls[i] .+ 1/2
σ22_p_t = sum(σ22_t(i) for i = 1:M_p) ./ M_p
σ22_np_t = sum(σ22_t(i) for i = (M_p+1):(M_p+M_np)) ./ M_np

p = plot(xlabel = "Γ1D t ⋅ 10³", ylabel = "⟨σ²²⟩")
plot!(p, t, σ22_p_t, label = "pumped")
plot!(p, t, σ22_np_t, label = "non-pumped")
plot(p, size = (500, 300))
````

## Package versions

These results were obtained using the following versions:

````@example waveguide
using InteractiveUtils
versioninfo()

using Pkg
Pkg.status(["SummationByPartsOperators", "OrdinaryDiffEq"], mode = PKGMODE_MANIFEST)
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

