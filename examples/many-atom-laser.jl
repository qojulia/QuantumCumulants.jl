# # Many-atom laser

# This example describes a second order laser system consisting of $N$ three-level atoms coupled to a single mode cavity. An auxiliary state $|3\rangle$, which quickly decays into the upper lasing state $|2\rangle$, is coherently pumped to achieve population inversion on the lasing transition $|1\rangle \leftrightarrow |2\rangle$. The Hamiltonian of this system is given by

# ```math
# H = -\Delta_{c} a^{\dagger} a  - \sum_{i=1}^N \left[ \Delta_3^i \sigma_i^{33}  + g_i (a^{\dagger} \sigma_i^{12} + a\sigma_i^{21}) + \Omega_i (\sigma_i^{31} + \sigma_i^{13}) \right].
# ```

# Including dissipative processes as, e.g. the atomic decay or photon losses through the cavity mirrors, makes it an open quantum system. In the Schrödinger picture we would compute the dynamics of such open quantum systems with a density matrix $\rho$ according to a master equation (see e.g. [https://docs.qojulia.org/](https://docs.qojulia.org/)),

# $\frac{d}{dt} \rho = - \frac{i}{\hbar} \left[ H, \rho \right] + \mathcal{L}[\rho],$

# with $\mathcal{L}[\rho] = \frac{\gamma}{2} (2 J \rho J^\dagger - J^\dagger J \rho - \rho J^\dagger J)$ the Liouvillian superoperator in standard Lindblad form for a dissipative process with jump operator $J$ and rate $R$.

# With **QuantumCumulants.jl** we describe the system dynamics with averages, which are deduced from the operator equations of motion in the Heisenberg picture. In the Heisenberg picture open systems are described by the quantum Langevin equation. Assuming white noise, we can omit the stochastic terms of the quantum Langevin equation when computing averages. Thus we get the following equation for the time evolution of a system operator average $\langle O \rangle$ (if $O$ is not explicitly time dependent):

# $\frac{d}{dt} \langle O \rangle = \frac{i}{\hbar} \left[ H, O \right] + \bar{\mathcal{L}}[O].$

# The superoperator $\bar{\mathcal{L}}[O]$ is similar to the Lindblad term in the Schrödinger picture, except that $J$ and $J^\dagger$ are swapped in the first term, i.e. $\bar{\mathcal{L}}[O] = \frac{\gamma}{2} (2 J^\dagger O J - J^\dagger J O - O J^\dagger J)$, for a dissipative process with jump operator $J$ and rate $R$.

# For our system we have four different dissipative processes with the jump operators $a$, $\sigma^{12}_i$, $\sigma^{13}_i$ and $\sigma^{23}_i$, and the corresponding decay rates $\kappa$, $\Gamma^i_{12}$, $\Gamma^i_{13}$ and $\Gamma^i_{23}$, respectively.

# We start by loading the needed packages.

using QuantumCumulants
using ModelingToolkitBase, OrdinaryDiffEqTsit5
using Plots

# Then we define the symbolic parameters of the system, the Hilbertspace and the necessary operators. We define an atomic transition operator function $\sigma(i,j,k)$ for the transition from $|j \rangle$ to $|i \rangle$ of atom $k$. Since we only have one [`FockSpace`](@ref) we do not need to specify the Hilbertspace on which the [`Destroy`](@ref) operator acts. For the different atomic transitions, however, we need to specify this, since there is more than one [`NLevelSpace`](@ref). This information is stored in the `.aon` field of each operator.

N = 2 # number of atoms
@variables κ g Γ23 Γ13 Γ12 Ω Δc Δ3

hf = FockSpace(:cavity) # Hilbertspace
ha = ⊗([NLevelSpace(Symbol(:atom, i), 3) for i in 1:N]...)
h = hf ⊗ ha

a = Destroy(h, :a) # Operators
σ(i, j, k) = Transition(h, Symbol("σ_{$k}"), i, j, k + 1)
nothing # hide

# Now we create the Hamiltonian and the jumps with the corresponding rates of our laser system. We assume here that all atoms are identical.

H =
    -Δc * a'a +
    sum(g * (a' * σ(1, 2, i) + a * σ(2, 1, i)) for i in 1:N) +
    sum(Ω * (σ(3, 1, i) + σ(1, 3, i)) for i in 1:N) - sum(Δ3 * σ(3, 3, i) for i in 1:N) # Hamiltonian

J = [a; [σ(1, 2, i) for i in 1:N]; [σ(1, 3, i) for i in 1:N]; [σ(2, 3, i) for i in 1:N]] # Jumps

rates = [κ; [Γ12 for i in 1:N]; [Γ13 for i in 1:N]; [Γ23 for i in 1:N]] # Rates
nothing # hide

# Later we will complete the system automatically, which has the disadvantage that the equations are not ordered. Therefore we define a list of interesting operators, which we want to use later. Note that at least one operator(-product) is needed. We derive the equations for these operators, average them, and automatically complete the system of equations.

ops = [a'a, σ(2, 2, 1), σ(3, 3, 1)] # list of operators

eqs = meanfield(ops, H, J; rates = rates, order = 2) #second order average

complete!(eqs) # automatically complete the system
nothing # hide

# To calculate the time evolution we create a Julia function which can be used by DifferentialEquations.jl to solve the set of ordinary differential equations.

# Build a System out of the MeanfieldEquations
sys = System(eqs; name = :laser)
sys_c = mtkcompile(sys)
nothing # hide

# Finally, we compute the time evolution after defining an initial state and numerical values for the parameters.

u0 = initial_values(eqs) # initial state

Γ12n = 1.0
Γ23n = 20Γ12n
Γ13n = 2Γ12n
Ωn = 5Γ13n
gn = 2Γ12n
Δcn = 0.0
Δ3n = 0.0
κn = 0.5Γ12n

ps = (g, Γ23, Γ13, Γ12, Ω, Δc, Δ3, κ) # list of parameters
p0 = Dict(ps .=> (gn, Γ23n, Γ13n, Γ12n, Ωn, Δcn, Δ3n, κn))
tend = 10.0 / κn

prob = ODEProblem(sys_c, merge(u0, p0), (0.0, tend))
sol = solve(prob, Tsit5(), reltol = 1.0e-8, abstol = 1.0e-8)
nothing # hide

# We plot the average photon number and the population inversion of the lasing transition.

n_t = real.(get_solution(sol, a' * a, eqs).(sol.t))
σ22_t = real.(get_solution(sol, σ(2, 2, 1), eqs).(sol.t))
σ33_t = real.(get_solution(sol, σ(3, 3, 1), eqs).(sol.t))
σ22m11_t = 2 .* σ22_t .+ σ33_t .- 1 #σ11 + σ22 + σ33 = 𝟙

p1 = plot(sol.t, n_t, xlabel = "tΓ₁₂", ylabel = "⟨a⁺a⟩", legend = false) # Plot
p2 = plot(sol.t, σ22m11_t, xlabel = "tΓ₁₂", ylabel = "⟨σ22⟩ - ⟨σ11⟩", legend = false)
plot(p1, p2, layout = (1, 2), size = (800, 320), left_margin = 5Plots.mm, bottom_margin = 5Plots.mm, dpi = 300)

# ## Package versions

# These results were obtained using the following versions:

using InteractiveUtils
versioninfo()

using Pkg
Pkg.status(
    ["QuantumCumulants", "OrdinaryDiffEqTsit5", "ModelingToolkitBase", "Plots"],
    mode = PKGMODE_MANIFEST,
)
