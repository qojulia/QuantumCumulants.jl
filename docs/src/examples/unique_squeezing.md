# Unique Steady-State Squeezing

In this example we show the unique squeezing observed in a driven Dicke model described by $N$ two-level systems coupled to a quantized harmonic oscillator. First we present the full dynamics with a second order cumulant expansion. The Hamiltonian describing the system is

```math
\begin{align}
H = \omega a^\dagger a + \frac{\Omega}{2} \sum_j  \sigma^j_z + \frac{g}{2} \sum_j  (a^\dagger + a) \sigma^j_x + \eta ( a \, e^{i \omega_\mathrm{d} t} + a^\dagger e^{-i \omega_\mathrm{d} t}),
\end{align}
```

for $N = 1$ it describes the driven quantum Rabi model. Additionally the system features two decay channels, losses of the harmonic oscillator with rate $\kappa$ and relaxation of the two-level system with rate $\gamma$.

We start by loading the packages.


```@example unique_squeezing
using QuantumCumulants
using OrdinaryDiffEq, ModelingToolkit
using Plots
nothing # hide
```

We define the Hilbert space and the symbolic parameters of the system.


```@example unique_squeezing
# Define hilbert space
hf = FockSpace(:harmonic)
ha = NLevelSpace(Symbol(:spin),2)
h = hf ⊗ ha

# Paramter
@cnumbers ω Ω ωd η κ g γ N
@syms t::Real # time
nothing # hide
```

On the Hilbert space we create the destroy operator $a$ of the harmonic oscillator and the (indexed) transition operator $\sigma_i^{xy}$ for the $i$-th two-level system.


```@example unique_squeezing
@qnumbers a::Destroy(h)
σ(x,y,i) = IndexedOperator(Transition(h,:σ,x,y),i)
nothing # hide
```

With the symbolic parameters, operators and indices we define the Hamiltonian and Liouvillian of the system.


```@example unique_squeezing
# Indices
i = Index(h,:i,N,ha)
j = Index(h,:j,N,ha)

# Hamiltonian
Hf =  ω*a'*a + η*(a'*exp(-1im*ωd*t) + a*exp(1im*ωd*t) )
Ha =  Ω*Σ(σ(2,2,i)-σ(1,1,i),i)/2
Hi =  g*Σ((σ(1,2,i)+σ(2,1,i))*(a + a'),i)/2
H = Hf + Ha + Hi

# Jump operators & and rates
J = [a, σ(1,2,i)]
rates = [κ, γ]
nothing # hide
```

First we derive the mean-field equations in second order for $\langle a \rangle$, $\langle a^\dagger a \rangle$ and $\langle \sigma^{22}_j \rangle$, then we complete the system to obtain a closed set of equations.


```@example unique_squeezing
eqs = meanfield([a, a'a, σ(2,2,j)],H,J;rates=rates,order=2)
nothing # hide
```

```math
\begin{align}
\frac{d}{dt} \langle a\rangle  =& -0.5 i \left( \underset{i}{\overset{N}{\sum}} g  \langle {\sigma}_{i}^{{12}}\rangle  + \underset{i}{\overset{N}{\sum}} g  \langle {\sigma}_{i}^{{21}}\rangle  \right) -1 i \eta e^{-1 i t {\omega}d} -0.5 \kappa \langle a\rangle  -1 i \omega \langle a\rangle  \\
\frac{d}{dt} \langle a^\dagger  a\rangle  =& 0.5 i \left( \underset{i}{\overset{N}{\sum}} g  \langle a  {\sigma}_{i}^{{12}}\rangle  + \underset{i}{\overset{N}{\sum}} g  \langle a  {\sigma}_{i}^{{21}}\rangle  \right) -0.5 i \left( \underset{i}{\overset{N}{\sum}} g  \langle a^\dagger  {\sigma}_{i}^{{12}}\rangle  + \underset{i}{\overset{N}{\sum}} g  \langle a^\dagger  {\sigma}_{i}^{{21}}\rangle  \right) -1.0 \kappa \langle a^\dagger  a\rangle  -1 i \eta \langle a^\dagger\rangle  e^{-1 i t {\omega}d} + 1 i \eta \langle a\rangle  e^{1 i t {\omega}d} \\
\frac{d}{dt} \langle {\sigma}_{j}^{{22}}\rangle  =& -1.0 \gamma \langle {\sigma}_{j}^{{22}}\rangle  -0.5 i g \left( \langle a^\dagger  {\sigma}_{j}^{{21}}\rangle  + \langle a  {\sigma}_{j}^{{21}}\rangle  \right) + 0.5 i g \left( \langle a^\dagger  {\sigma}_{j}^{{12}}\rangle  + \langle a  {\sigma}_{j}^{{12}}\rangle  \right)
\end{align}
```


```@example unique_squeezing
eqs_c = complete(eqs)
length(eqs_c)
```

All two-level systems behave identically, due to this permutation symmetry of the system we can scale-up the equations.


```@example unique_squeezing
eqs_sc = scale(eqs_c)
scale(eqs) # Example scaling on the first three equations
nothing # hide
```

```math
\begin{align}
\frac{d}{dt} \langle a\rangle  =& -0.5 i \left( N g \langle {\sigma}_{1}^{{12}}\rangle  + N g \langle {\sigma}_{1}^{{21}}\rangle  \right) -1 i \eta e^{-1 i t {\omega}d} -0.5 \kappa \langle a\rangle  -1 i \omega \langle a\rangle  \\
\frac{d}{dt} \langle a^\dagger  a\rangle  =& -0.5 i \left( N g \langle a^\dagger  {\sigma}_{1}^{{12}}\rangle  + N g \langle a^\dagger  {\sigma}_{1}^{{21}}\rangle  \right) + 0.5 i \left( N g \langle a  {\sigma}_{1}^{{12}}\rangle  + N g \langle a  {\sigma}_{1}^{{21}}\rangle  \right) -1.0 \kappa \langle a^\dagger  a\rangle  -1 i \eta \langle a^\dagger\rangle  e^{-1 i t {\omega}d} + 1 i \eta \langle a\rangle  e^{1 i t {\omega}d} \\
\frac{d}{dt} \langle {\sigma}_{1}^{{22}}\rangle  =& 0.5 i g \left( \langle a^\dagger  {\sigma}_{1}^{{12}}\rangle  + \langle a  {\sigma}_{1}^{{12}}\rangle  \right) -0.5 i g \left( \langle a^\dagger  {\sigma}_{1}^{{21}}\rangle  + \langle a  {\sigma}_{1}^{{21}}\rangle  \right) -1.0 \gamma \langle {\sigma}_{1}^{{22}}\rangle
\end{align}
```

To calculate the dynamics of the system we create a system of ordinary differential equations with its initial state and numerical parameters.


```@example unique_squeezing
# symbolic ordinary differential equation system
@named sys = ODESystem(eqs_sc)

# initial state
u0 = zeros(ComplexF64, length(eqs_sc));

# Parameters
ω_ = 1.0
Ω_ = 2e3ω_
N_ = 1
gc_ = sqrt(Ω_*ω_/N) # renormalization of coupling to keep the system intensive
g_ = 0.9gc_
η_ = 4ω_
κ_ = ω_
γ_ = ω_
ωd_ = sqrt(1-g_^2/gc_^2)*ω_

# symbolic and numeric parameter list
ps = [ω , Ω , ωd , g , η , κ , γ , N ]
p0 = [ω_, Ω_, ωd_, g_, η_, κ_, γ_, N_]
nothing # hide
```

We solve the dynamics for four different numbers of two-level systems $N = [1, 10, 20, 100]$.


```@example unique_squeezing
sol_ls = []
N_ls = [1,2,10,100]
for N_ in N_ls
    p0_ = [ω_, Ω_, ωd_, g_, η_, κ_, γ_, N_]
    prob = ODEProblem(sys,u0,(0.0, 4π/ωd_), ps.=>p0_)
    sol = solve(prob,Tsit5(),reltol=1e-10,abstol=1e-10)
    push!(sol_ls,sol)
end
```


```@example unique_squeezing
# plot results
c_ls=[:black, :red, :blue, :cyan]
p1 = plot(xlabel="ω t", ylabel="Δ² O")
p2 = plot(xlabel="ω t", ylabel="⟨σz⟩")
for i=1:length(N_ls)
    sol = sol_ls[i]
    t_ = sol.t

    sqx = sol[a'*a'] + sol[a*a] + 2*sol[a'*a] .+ 1 - (sol[a'] + sol[a]).^2
    sqy = sol[a'*a'] + sol[a*a] - 2*sol[a'*a] .- 1 - (sol[a'] - sol[a]).^2
    plot!(p1,t_,real.(sqx),label="N = $(N_ls[i])",color=c_ls[i])
    plot!(p1,t_,-real.(sqy),ls=:dash,label=nothing,color=c_ls[i])

    s22 = sol[σ(2,2,1)]
    plot!(p2,t_,real.(2s22 .- 1),color=c_ls[i],label=nothing)
end
plot(p1, p2, layout=(1,2), size=(700,250),bottom_margin=5*Plots.mm, left_margin=5*Plots.mm)
```

## Effective model

For a suffeciently low excitation we can adiabatically elminate the dynamics of the two-level system(s). This leads to an effective Hamiltonian

```math
\begin{align}
H_\mathrm{a} = \omega a^\dagger a - \frac{g^2}{4 \Omega}(a + a^\dagger)^2 + \eta ( a \, e^{i \omega_\mathrm{d} t} + a^\dagger e^{-i \omega_\mathrm{d} t}).
\end{align}
```

We calculate now the dynamics for this effective model and compare it with the full system. Note that this Hamiltonian is quadratic, which means that a second order description is exact.


```@example unique_squeezing
# effective Hamiltonian
@cnumbers gΩ # g^2/4Ω
H_a = Hf - gΩ*(a + a')^2

eqs_a = meanfield([a, a'a, a*a],H_a,[a];rates=[κ],order=2)
nothing # hide
```

```math
\begin{align}
\frac{d}{dt} \langle a\rangle  =& 2 i g\Omega \left( \langle a^\dagger\rangle  + \langle a\rangle  \right) -1 i \eta e^{-1 i t {\omega}d} -0.5 \kappa \langle a\rangle  -1 i \omega \langle a\rangle  \\
\frac{d}{dt} \langle a^\dagger  a\rangle  =& -1.0 \kappa \langle a^\dagger  a\rangle  + 2 i g\Omega \langle a^\dagger  a^\dagger\rangle  -2 i g\Omega \langle a  a\rangle  -1 i \eta \langle a^\dagger\rangle  e^{-1 i t {\omega}d} + 1 i \eta \langle a\rangle  e^{1 i t {\omega}d} \\
\frac{d}{dt} \langle a  a\rangle  =& 2 i g\Omega + 4 i g\Omega \left( \langle a^\dagger  a\rangle  + \langle a  a\rangle  \right) -1.0 \kappa \langle a  a\rangle  -2 i \omega \langle a  a\rangle  -2 i \eta \langle a\rangle  e^{-1 i t {\omega}d}
\end{align}
```


```@example unique_squeezing
# symbolic ordinary differential equation system
@named sys_a = ODESystem(eqs_a)

# initial state
u0_a = zeros(ComplexF64, length(eqs_a))

# Additional parameter
gΩ_ = g_^2/(4Ω_)

# symbolic and numeric parameter list
ps_a = [ω , ωd , η , κ , N , gΩ ]
p0_a = [ω_, ωd_, η_, κ_, N_, gΩ_]

# define and solve numeric ordinary differential equation problem
prob_a = ODEProblem(sys_a,u0_a,(0.0, 4π/ωd_), ps_a.=>p0_a)
sol_a = solve(prob_a,Tsit5(),reltol=1e-8,abstol=1e-8)
nothing # hide
```


```@example unique_squeezing
# plot results
sol = sol_ls[4]
t_ = sol.t
sqx = sol[a'*a'] + sol[a*a] + 2*sol[a'*a] .+ 1 - (sol[a'] + sol[a]).^2
sqy = sol[a'*a'] + sol[a*a] - 2*sol[a'*a] .- 1 - (sol[a'] - sol[a]).^2

t_a = sol_a.t
sqx_a = sol_a[a'*a'] + sol_a[a*a] + 2*sol_a[a'*a] .+ 1 - (sol_a[a'] + sol_a[a]).^2
sqy_a = sol_a[a'*a'] + sol_a[a*a] - 2*sol_a[a'*a] .- 1 - (sol_a[a'] - sol_a[a]).^2

p = plot(xlabel="ω t", ylabel="Δ² O")
plot!(p,t_,real.(sqx),label="X - Full model")
plot!(p,t_,-real.(sqy),label="P - Full model",ls=:dash)
plot!(p,t_a,real.(sqx_a),label="X - Effective model")
plot!(p,t_a,-real.(sqy_a),label="P - Effective model",ls=:dash)
plot(p, size=(500,200))
```
