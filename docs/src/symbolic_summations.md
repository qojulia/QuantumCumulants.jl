# Symbolic Indexing

## Implementation

### Index

The main tool one needs to define first, before using any type of symbolic summations is the **Index** object. This object consists of four different fields, that all need to be specified upon construction. These fields consist of a [`HilbertSpace`](@ref), a **name**, which is just a generic Symbol, a **Range**, which can either consist of again a Symbol or a concrete number and a specific [`HilbertSpace`](@ref), which defines the space on which operators, that inherit the **Index** entity act on. This in specific means that an **Index** entity, that is acting in an [`NLevelSpace`](@ref) can only be inherited by [`Transition`](@ref) operators. In the example below, two indices are defined equivalently, aswell as a third one being defined acting on the [`FockSpace`](@ref) of the defined [`ProductSpace`](@ref) **h**.


```@example symbolic_summations
using QuantumCumulants

@cnumbers N

ha = NLevelSpace(:atoms,2)
hc = FockSpace(:cavity)
h = hc ⊗ ha

i = Index(h,:i,N,ha)
i2 = Index(h,:i,N,2) #equivalent definition

n = Index(h,:n,5,hc)
nothing #hide
```

### IndexedOperators

Operators, such as the quantum harmonic destruction operator [`Destroy`](@ref) or the [`Transition`](@ref) operator can be associated with an [`Index`](@ref) of the corresponding Hilbertspace by creating a so called [`IndexedOperator`](@ref). This new object consists of two fields, namely the operator itself and an Index. Below, there are two **IndexedOperator** entities created on the two different hilbertspaces, defined previously.


```@example symbolic_summations
σ(x,y,z) = IndexedOperator(Transition(h,:σ,x,y),z)
a(z) = IndexedOperator(Destroy(h,:a),z)
nothing #hide
```




In the above example, we defined both IndexedOperators **σ** and **a** as callable Instances with the attribute-variable **z**. These can now be used to easily create operators, that act spcifically with their associated Index.


```@example symbolic_summations
a(n)*σ(2,2,i)
nothing #hide
```

```math
{a}_{n} {\sigma}_{i}^{{22}}
```


Similar to Operators, one can also create so called [`IndexedVariable`](@ref) objects, which consist simply of a name and an [`Index`](@ref). These can be used to store variables that are associated with an Index.


```@example symbolic_summations
gi = IndexedVariable(:g,i)
nothing #hide
```

### Summations

As for now, we only created single Instances of indexed operators. These operators and variables can now be used to define symbolic summations, which can then again be used in defining a Hamiltonian and deriving equations of motion for specific opearator averages.

Such a summation needs two arguments to be constructed, the **term**, over which the summation shall sum over, and an [`Index`](@ref), over which the sum runs. As an example, we define below a simple sum over a single IndexedOperator.


```@example symbolic_summations
∑(σ(2,2,i),i)
nothing #hide
```
```math
\underset{i}{\overset{N}{\sum}}{σ}_{i}^{{22}}
```

As can be seen above, a sum with a single running-index can be created using the **∑** (\sum) command. Other equivalent functions are **Σ** (\Sigma) and the **SingleSum()** constructor. Similar to this one can also create summations over up to two different running-indices:


```@example symbolic_summations
∑(a(n)*σ(2,1,i),i,n)
nothing #hide
```

```math
\underset{i}{\overset{N}{\sum}} \underset{n}{\overset{5}{\sum}} {a}_{n}  {σ}_{i}^{{21}}
```


These two running-indices do not need to act on different Hilbertspaces, however, when one chooses indices, acting on the same Hilbertspace, one can observe, that an immediate simplification occurs, as shown below.


```@example symbolic_summations
k = Index(h,:k,N,ha) # some additional Indices
l = Index(h,:l,N,ha) 

∑(σ(2,1,k)*σ(1,2,l),k,l)
nothing #hide
```
```math
\underset{k{\ne}l}{\overset{N}{\sum}} \underset{l{\ne}k}{\overset{N}{\sum}} {\sigma}_{l}^{{12}}  {\sigma}_{k}^{{21}} + \underset{k}{\overset{N}{\sum}} {\sigma}_{k}^{{22}}
```


What happened above, is a simple simplification, that occurs when two indices, acting on the same Hilbertspace meet each other inside of a summation. The special case, where the numeric values of both indices can be the same, i.e `l`=`k`, is calculated immediately. This can also be observed, when a symbolic sum is multiplied with an [`IndexedOperator`](@ref), that is acting on the same Hilbertspace as the summation-index.


```@example symbolic_summations
∑(σ(2,2,k),k) * σ(2,1,l) 
nothing #hide
```
```math
\underset{k{\ne}l}{\overset{N}{\sum}} {\sigma}_{k}^{{22}}  {\sigma}_{l}^{{21}} + {\sigma}_{l}^{{21}}
```

At the current state of devolopment, the possibility to create summations with 3 or more running indices is not yet implemented.

## Short Example

As a short example, we will briefly go over the entire process of defining a Hamiltonian over the derivation of equations up to solving these equations numerically.

For this example we will consider **N** 2-level atoms in a single mode Cavity.

We start by defining all Indices and Operators we need.


```@example symbolic_summations
using QuantumCumulants

ha = NLevelSpace(:atoms,2)
hc = FockSpace(:cavity)
h = hc ⊗ ha

@cnumbers N Δ κ γ ν

i = Index(h,:i,N,ha)
j = Index(h,:j,N,ha)

@qnumbers b::Destroy(h)
σ(x,y,z) = IndexedOperator(Transition(h,:σ,x,y),z)
gi = IndexedVariable(:g,i)

H = Δ*b'*b + ∑(gi*(b*σ(2,1,i) + b'*σ(1,2,i)),i)
nothing #hide
```

We continue by defining the starting operators **ops**, for which the first set of equations is calculated, aswell as the Jump operators **J** with their corresponding rates. It is important to note here, that for IndexedOperators, for which these equations are calculated need to have an [`Index`](@ref), which is not yet used in a summation within the Hamiltonian **H**. We can then create the first set by simply calling the **meanfield** function.


```@example symbolic_summations
ops = [b'b, σ(2,2,j)]
J = [b, σ(1,2,i),σ(2,1,i)] 
rates = [κ,γ,ν]

eqs = meanfield(ops,H,J;rates=rates,order=2)
nothing #hide
```

```math
\begin{align}
\frac{d}{dt} \langle b^\dagger  b\rangle  =& 1 i \underset{i}{\overset{N}{\sum}} {g}_{i}  \langle b  {\sigma}_{i}^{{21}}\rangle  -1 i \underset{i}{\overset{N}{\sum}} {g}_{i}  \langle b^\dagger  {\sigma}_{i}^{{12}}\rangle  -1.0 \kappa \langle b^\dagger  b\rangle  \\
\frac{d}{dt} \langle {\sigma}_{j}^{{22}}\rangle  =& \nu -1.0 \gamma \langle {\sigma}_{j}^{{22}}\rangle  -1.0 \nu \langle {\sigma}_{j}^{{22}}\rangle  + 1 i {g}_{j} \langle b^\dagger  {\sigma}_{j}^{{12}}\rangle  -1 i {g}_{j} \langle b  {\sigma}_{j}^{{21}}\rangle 
\end{align}
```

We can then complete the set of equations by simply calling the **complete** function. When completing equations, that contain indexed operators or variables, additional indices, that are not yet in use, are introduces automatically for calculation. One can also individually set these indices via the keyword **extra_indices**, which takes either a set of Symbols, or specifically defined indices. 


```@example symbolic_summations
eqs_comp = complete(eqs)
nothing #hide
```

```math
\begin{align}
\frac{d}{dt} \langle b^\dagger  b\rangle  =& 1 i \underset{i}{\overset{N}{\sum}} {g}_{i}  \langle b  {\sigma}_{i}^{{21}}\rangle  -1 i \underset{i}{\overset{N}{\sum}} {g}_{i}  \langle b^\dagger  {\sigma}_{i}^{{12}}\rangle  -1.0 \kappa \langle b^\dagger  b\rangle  \\
\vdots
\end{align}
```

When we take a look now at one of the equations in the completed system, we can see, that the equation still contains symbolic summations and indices. Also up until now, we have not specified a numerical value for the upper boundaries we use for summations, meaning the used **N** still has no numerical assignment. To create now numerically solvable equations we can simply use the [`evaluate`](ref) function. Within this function we also specify the numerical value for **N**, in this case we chose 3 atoms.


```@example symbolic_summations
evaled = evaluate(eqs_comp;limits=(N=>3))
nothing #hide
```

```math
\begin{align}
\frac{d}{dt} \langle {\sigma}_{1}^{{22}}\rangle  =& \nu + 1 i g_{1} \langle b^\dagger  {\sigma}_{1}^{{12}}\rangle  -1.0 \gamma \langle {\sigma}_{1}^{{22}}\rangle  -1.0 \nu \langle {\sigma}_{1}^{{22}}\rangle  -1 i g_{1} {\langle b^\dagger  {\sigma}_{1}^{{12}}\rangle ^{*}} \\
\vdots
\end{align}
```

The last set of equations are now in a numerical solvable form, that we can convert to an `ODESystem` as defined in [ModelingToolkit](https://github.com/SciML/ModelingToolkit.jl) which can be solved numerically with [OrdinaryDiffEq](https://github.com/JuliaDiffEq/OrdinaryDiffEq.jl). We are further more utilizing the fact, that we declared `gi` as a [`IndexedVariable`](@ref) to give each atom a different coupling strength `g`. This can now be done by using the [`value_map`](@ref) function to create a parameter mapping for the `ODEProblem`. In this case we give as coupling strength a vector with different values.


```@example symbolic_summations
# Generate an ODESystem
using ModelingToolkit
@named sys = ODESystem(evaled)

# Solve the system using the OrdinaryDiffEq package
using OrdinaryDiffEq
u0 = zeros(ComplexF64,length(evaled))
p = [Δ, gi, γ, κ, ν]
p0 = [0,[0.75,1.2,1.5],0.25,1,1.5]
p_ = value_map(p,p0;limits=(N=>3))
prob = ODEProblem(sys,u0,(0.0,10.0),p_)
sol = solve(prob,RK4())
nothing #hide
```

Just as with variables in [ModelingToolkit](https://github.com/SciML/ModelingToolkit.jl), the solution of the respective averages can be accessed with a `getindex` method. In the following we extract and plot the photon number and the atomic excited state population of the 1st atom by indexing the solution:


```@example symbolic_summations
using Plots
n = real.(sol[b'*b])
pe = [real.(sol[σ(2,2,i)]) for i = 1:3]
pl = plot(sol.t, n, label="Photon number", xlabel="t")
for i = 1:3
    plot!(sol.t, pe[i], label="Excited state population of atom $(i)")
end
pl
savefig("symbolic_summation.svg") # hide
nothing # hide
```
    
![Photon number and excited state populations](symbolic_summation.svg)
    


