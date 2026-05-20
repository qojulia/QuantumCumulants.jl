# # Cavity Antiresonance

# In this example we investigate a system of $N$ closely spaced quantum emitters inside a coherently driven single mode cavity. The model is described in [D. Plankensteiner, et. al., Phys. Rev. Lett. 119, 093601 (2017)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.119.093601).
# The Hamiltonian of this system is composed of three parts $H = H_c + H_a + H_{\mathrm{int}}$, the driven cavity $H_c$, the dipole-dipole interacting atoms $H_a$ and the atom-cavity interaction $H_\mathrm{int}$:

# ```math
# \begin{align}
# H_\mathrm{c} &= \hbar \Delta_c a^\dagger a + \hbar \eta (a^\dagger + a) \\
# &\\
# H_a &= \hbar \Delta_a \sum\limits_{j} \sigma_j^{22} + \hbar \sum\limits_{i \neq j} \Omega_{ij} \sigma_i^{21} \sigma_j^{12}
# &\\
# H_\mathrm{int} &= \hbar \sum\limits_{j} g_j (a^\dagger \sigma_j^{12} + a \sigma_j^{21})
# \end{align}
# ```

# Additionally the system features two decay channels, the lossy cavity with photon decay rate $\kappa$ and collective atomic emission described by the decay-rate matrix $\Gamma_{ij}$.

# We start by loading the packages.

using QuantumCumulants
using OrdinaryDiffEq, ModelingToolkitBase
using Plots

# The Hilbert space for this system is given by one cavity mode and $N$ two-level atoms. Here we use symbolic indices, sums and double sums to define the system.
# The parameters $g_j, \, \Gamma_{ij}$ and $\Omega_{ij}$ are defined as indexed variables of atom $i$ and $j$. We will describe the system in first order mean-field.

hc = FockSpace(:cavity) # Hilbert space
ha = NLevelSpace(Symbol(:atom), 2)
h = hc ⊗ ha

@variables N Δc η Δa κ # Parameters
g(i) = IndexedVariable(:g, i)
Γ(i, j) = DoubleIndexedVariable(:Γ, i, j)
Ω(i, j) = DoubleIndexedVariable(:Ω, i, j; identical = false)


i = Index(h, :i, N, ha) # Indices
j = Index(h, :j, N, ha)


# The kwarg ’identical=false’ for the double indexed variable specifies that $\Omega_{ij} = 0$ for $i = j$.
# Now we create the operators on the composite Hilbert space using the $\texttt{IndexedOperator}$ constructor, which assigns each $\texttt{Transition}$ operator an $\texttt{Index}$.

@qnumbers a::Destroy(h)
σ(x, y, k) = IndexedOperator(Transition(h, :σ, x, y), k)
nothing # hide

# We define the Hamiltonian and Liouvillian. For the collective atomic decay we write the corresponding dissipative processes with a double indexed variable $R_{ij}$ and an indexed jump operator $J_j$, such that an operator average $\langle \mathcal{O} \rangle$ follows the equation

# ```math
# \begin{equation}
# \langle \dot{\mathcal{O}} \rangle = \sum_{ij} R_{ij} \left( \langle J_i^\dagger \mathcal{O} J_j \rangle - \frac{1}{2} \langle J_i^\dagger J_j \mathcal{O} \rangle - \frac{1}{2} \langle \mathcal{O} J_i^\dagger J_j \rangle \right).
# \end{equation}
# ```

# The inner dipole-dipole sum excludes the diagonal `i == j` by passing the
# `non_equal` vector `[i]` to the inner `Σ` (SQA v0.5 replaced the old
# `non_equal=true` keyword with this explicit form).
Hc = Δc * a'a + η * (a' + a) # Hamiltonian
Ha = Δa * Σ(σ(2, 2, i), i) + Σ(Σ(Ω(i, j) * σ(2, 1, i) * σ(1, 2, j), j, [i]), i)
Hi = Σ(g(i) * (a' * σ(1, 2, i) + a * σ(2, 1, i)), i)
H = Hc + Ha + Hi

J = [a, σ(1, 2, i)] # Jump operators
rates = [κ, Γ(i, j)]
nothing # hide

# We derive the system of equations in first order mean-field.

eqs = meanfield(a, H, J; rates = rates, order = 1)
complete!(eqs)
nothing # hide

# To create the equations for a specific number of atoms we use the function $\texttt{evaluate}$.
# In v1.0 the `evaluate(eqs; limits=(N=>n,))` pass that unrolls the indexed
# sums for a fixed $N$ is still being ported (see TODO.md), so this port
# stops at the symbolic equations. The numeric transmission sweep below will
# be re-enabled once `evaluate` lands.

# N_ = 2
# eqs_ = evaluate(eqs; limits = (N=>N_))
# @named sys = System(eqs_)
