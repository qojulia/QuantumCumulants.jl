# # Unique Steady-State Squeezing

# In this example we show the unique squeezing observed in a driven Dicke model described by $N$ two-level systems coupled to a quantized harmonic oscillator [[K. Gietka. et. al., Phys. Rev. Lett. 131, 223604 (2023)]](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.131.223604).  First we present the full dynamics with a second order cumulant expansion. The Hamiltonian describing the system is

# ```math
# \begin{align}
# H = \omega a^\dagger a + \frac{\Omega}{2} \sum_j  \sigma^j_z + \frac{g}{2} \sum_j  (a^\dagger + a) \sigma^j_x + \eta ( a \, e^{i \omega_{d}\, t} + a^\dagger e^{-i \omega_{d} \,t}),
# \end{align}
# ```

# for $N = 1$ it describes the driven quantum Rabi model. Additionally the system features a decay channel, losses of the harmonic oscillator with rate $\kappa$.

# We start by loading the packages.


using QuantumCumulants
using OrdinaryDiffEq, ModelingToolkitBase
using Plots
nothing # hide

# We define the Hilbert space and the symbolic parameters of the system.

hf = FockSpace(:harmonic) # Define hilbert space
ha = NLevelSpace(Symbol(:spin), 2)
h = hf ⊗ ha

@variables ω Ω ωd η κ g γ N ξ # Parameter
t = first(ModelingToolkitBase.@independent_variables t) # time iv used by meanfield internally
nothing # hide

# On the Hilbert space we create the destroy operator $a$ of the harmonic oscillator and the (indexed) transition operator $\sigma_i^{xy}$ for the $i$-th two-level system.

@qnumbers a::Destroy(h)
σ(x, y, i) = IndexedOperator(Transition(h, :σ, x, y), i)
nothing # hide

# With the symbolic parameters, operators and indices we define the Hamiltonian and Liouvillian of the system. Note, however, that in the strong coupling regime $g \sim g_c\equiv \sqrt{\omega \Omega}$ the driving term and jump operators have to be redefined. For a strongly interacting system, the ground state is very different from the ground state of a non-interacting system. Therefore, using jump operators of a non-interacting system would lead to extraction of energy from the ground state of a strongly interacting system. The correct operators are the ones that diagonalise the Hamiltonian with adiabatiacally eliminated spins.

b = a * cosh(ξ) + a' * sinh(ξ) # Operators diagonalizing the Hamiltonian,i.e. approximate new eigenmodes of the system

i = Index(h, :i, N, ha) # Indices
j = Index(h, :j, N, ha)

Hf = ω * a' * a + η * (b' * exp(-1im * ωd * t) + b * exp(1im * ωd * t))
Ha = Ω * Σ(σ(2, 2, i) - σ(1, 1, i), i) / 2
Hi = g * Σ((σ(1, 2, i) + σ(2, 1, i)) * (a + a'), i) / 2
H = Hf + Ha + Hi # Hamiltonian

J = [b, σ(1, 2, i)] # Jump operators & and rates
rates = [κ, γ]

ps = [ω, Ω, ωd, g, η, κ, γ, N, ξ] # symbolic and numeric parameter list

# First we derive the mean-field equations in second order for $\langle a \rangle$, $\langle a^\dagger a \rangle$ and $\langle \sigma^{22}_j \rangle$, then we complete the system to obtain a closed set of equations.

eqs = meanfield([a, a'a, σ(2, 2, j)], H, J; rates = rates, order = 2, iv = t)
nothing # hide

# ```math
# \begin{align}
# \frac{d}{dt} \langle a\rangle  =& -0.5 i \left( \underset{i}{\overset{N}{\sum}} g  \langle {\sigma}_{i}^{{12}}\rangle  + \underset{i}{\overset{N}{\sum}} g  \langle {\sigma}_{i}^{{21}}\rangle  \right) -1 i \omega \langle a\rangle  -1 i \eta \sinh\left( \xi \right) e^{1 i t {\omega}d} -1 i \eta e^{-1 i t {\omega}d} {cosh(\xi)^{*}} -0.5 \kappa \cosh\left( \xi \right) \langle a\rangle  {cosh(\xi)^{*}} + 0.5 \kappa \sinh\left( \xi \right) \langle a\rangle  {sinh(\xi)^{*}} \\
# \frac{d}{dt} \langle a^\dagger  a\rangle  =& 0.5 i \left( \underset{i}{\overset{N}{\sum}} g  \langle a  {\sigma}_{i}^{{12}}\rangle  + \underset{i}{\overset{N}{\sum}} g  \langle a  {\sigma}_{i}^{{21}}\rangle  \right) -0.5 i \left( \underset{i}{\overset{N}{\sum}} g  \langle a^\dagger  {\sigma}_{i}^{{12}}\rangle  + \underset{i}{\overset{N}{\sum}} g  \langle a^\dagger  {\sigma}_{i}^{{21}}\rangle  \right) + \kappa \sinh\left( \xi \right) {sinh(\xi)^{*}} + \kappa \sinh\left( \xi \right) {sinh(\xi)^{*}} \langle a^\dagger  a\rangle  -1 i \eta \sinh\left( \xi \right) \langle a^\dagger\rangle  e^{1 i t {\omega}d} + 1 i \eta \cosh\left( \xi \right) \langle a\rangle  e^{1 i t {\omega}d} -1 i \eta \langle a^\dagger\rangle  e^{-1 i t {\omega}d} {cosh(\xi)^{*}} + 1 i \eta \langle a\rangle  e^{-1 i t {\omega}d} {sinh(\xi)^{*}} -1.0 \kappa \cosh\left( \xi \right) {cosh(\xi)^{*}} \langle a^\dagger  a\rangle  \\
# \frac{d}{dt} \langle {\sigma}_{j}^{{22}}\rangle  =& -1.0 \gamma \langle {\sigma}_{j}^{{22}}\rangle  -0.5 i g \left( \langle a^\dagger  {\sigma}_{j}^{{21}}\rangle  + \langle a  {\sigma}_{j}^{{21}}\rangle  \right) + 0.5 i g \left( \langle a^\dagger  {\sigma}_{j}^{{12}}\rangle  + \langle a  {\sigma}_{j}^{{12}}\rangle  \right)
# \end{align}
# ```

eqs_c = complete(eqs)
length(eqs_c)

# All two-level systems behave identically, due to this permutation symmetry of the system we can scale-up the equations.
eqs_sc = scale(eqs_c)
scale(eqs) # Example scaling on the first three equations
nothing # hide

# ```math
# \begin{align}
# \frac{d}{dt} \langle a\rangle  =& -0.5 i \left( N g \langle {\sigma}_{1}^{{21}}\rangle  + N g \langle {\sigma}_{1}^{{12}}\rangle  \right) -1 i \omega \langle a\rangle  -1 i \eta \sinh\left( \xi \right) e^{1 i t {\omega}d} -1 i \eta e^{-1 i t {\omega}d} {cosh(\xi)^{*}} -0.5 \kappa \cosh\left( \xi \right) \langle a\rangle  {cosh(\xi)^{*}} + 0.5 \kappa \sinh\left( \xi \right) \langle a\rangle  {sinh(\xi)^{*}} \\
# \frac{d}{dt} \langle a^\dagger  a\rangle  =& -0.5 i \left( N g \langle a^\dagger  {\sigma}_{1}^{{21}}\rangle  + N g \langle a^\dagger  {\sigma}_{1}^{{12}}\rangle  \right) + 0.5 i \left( N g \langle a  {\sigma}_{1}^{{21}}\rangle  + N g \langle a  {\sigma}_{1}^{{12}}\rangle  \right) + \kappa \sinh\left( \xi \right) {sinh(\xi)^{*}} + \kappa \sinh\left( \xi \right) {sinh(\xi)^{*}} \langle a^\dagger  a\rangle  -1 i \eta \sinh\left( \xi \right) \langle a^\dagger\rangle  e^{1 i t {\omega}d} + 1 i \eta \cosh\left( \xi \right) \langle a\rangle  e^{1 i t {\omega}d} -1 i \eta \langle a^\dagger\rangle  e^{-1 i t {\omega}d} {cosh(\xi)^{*}} + 1 i \eta \langle a\rangle  e^{-1 i t {\omega}d} {sinh(\xi)^{*}} -1.0 \kappa \cosh\left( \xi \right) {cosh(\xi)^{*}} \langle a^\dagger  a\rangle  \\
# \frac{d}{dt} \langle {\sigma}_{1}^{{22}}\rangle  =& -0.5 i g \left( \langle a^\dagger  {\sigma}_{1}^{{21}}\rangle  + \langle a  {\sigma}_{1}^{{21}}\rangle  \right) + 0.5 i g \left( \langle a^\dagger  {\sigma}_{1}^{{12}}\rangle  + \langle a  {\sigma}_{1}^{{12}}\rangle  \right) -1.0 \gamma \langle {\sigma}_{1}^{{22}}\rangle
# \end{align}
# ```

# To calculate the dynamics of the system we create a system of ordinary differential equations with its initial state and numerical parameters.


sys = System(eqs_sc; name = :sys) # symbolic ordinary differential equation system
sys = mtkcompile(sys)
u0 = zeros(ComplexF64, length(eqs_sc)); # initial state

ω_ = 1.0 # Parameters
Ω_ = 2.0e3ω_
η_ = 4ω_
κ_ = ω_
γ_ = ω_

# Derived quantities (`gc`, `g`, `ωd`, `ξ`) depend on `N`. We pick the
# largest `N` (the well-converged case) as the canonical global value used
# by the post-solve effective-model analysis below; the loop redefines them
# locally per `N_`.

N_global = 100
gc_ = sqrt(Ω_ * ω_ / N_global)
g_ = 0.9gc_
ωd_ = sqrt(1 - g_^2 / gc_^2) * ω_
ξ_ = 1 / 4 * log(1 - N_global * g_^2 / (ω_ * Ω_))

# We solve the dynamics for four different numbers of two-level systems $N = [1, 2, 10, 100]$.

sol_ls = []
N_ls = [1, 2, 10, 100]
for N_ in N_ls
    local gc_l = sqrt(Ω_ * ω_ / N_)
    local g_l = 0.9gc_l
    local ωd_l = sqrt(1 - g_l^2 / gc_l^2) * ω_
    local ξ_l = 1 / 4 * log(1 - N_ * g_l^2 / (ω_ * Ω_))
    local p0 = [ω_, Ω_, ωd_l, g_l, η_, κ_, γ_, N_, ξ_l]
    local u0_dict = Dict{Any, Any}(unknowns(sys) .=> u0)
    local dict = merge(u0_dict, Dict{Any, Any}(ps .=> p0))
    local prob = ODEProblem(sys, dict, (0.0, 4π / ωd_l))
    local sol = solve(prob, Tsit5(); saveat = π / 30ωd_l, reltol = 1.0e-10, abstol = 1.0e-10)
    push!(sol_ls, sol)
end

#

c_ls = [:black, :red, :blue, :cyan] # plot results
p1 = plot(xlabel = "ω t", ylabel = "Δ² O")
p2 = plot(xlabel = "ω t", ylabel = "⟨σz⟩")
for k in 1:length(N_ls)
    sol = sol_ls[k]
    t_ = sol.t

    adag_a = get_solution(sol, a' * a, eqs_sc).(t_)
    aa = get_solution(sol, a * a, eqs_sc).(t_)
    adag_adag = get_solution(sol, a' * a', eqs_sc).(t_)
    a_ = get_solution(sol, a, eqs_sc).(t_)
    adag = get_solution(sol, a', eqs_sc).(t_)

    sqx = adag_adag + aa + 2 * adag_a .+ 1 - (adag + a_) .^ 2
    sqy = adag_adag + aa - 2 * adag_a .- 1 - (adag - a_) .^ 2
    plot!(p1, t_, real.(sqx), label = "N = $(N_ls[k])", color = c_ls[k])
    plot!(p1, t_, -real.(sqy), ls = :dash, label = nothing, color = c_ls[k])

    s22 = get_solution(sol, σ(2, 2, i), eqs_sc).(t_)
    plot!(p2, t_, real.(2s22 .- 1), color = c_ls[k], label = nothing)
end
plot(
    p1,
    p2,
    layout = (1, 2),
    size = (700, 250),
    bottom_margin = 5 * Plots.mm,
    left_margin = 5 * Plots.mm,
)

# ## Effective model

# For a sufficiently low excitation we can adiabatically eliminate the dynamics of the two-level system(s). This leads to an effective Hamiltonian

# ```math
# \begin{align}
# H_\mathrm{a} = \omega a^\dagger a - \frac{g^2}{4 \Omega}(a + a^\dagger)^2 + \eta ( a \, e^{i \omega_{d} \, t} + a^\dagger e^{-i \omega_{d} \, t}).
# \end{align}
# ```

# We calculate now the dynamics for this effective model and compare it with the full system. Note that this Hamiltonian is quadratic, which means that a second order description is exact.

@variables gΩ # g^2/4Ω
H_a = Hf - N * gΩ * (a + a')^2 # effective Hamiltonian, N is added for the sake of intensitivity

eqs_a = meanfield([a, a'a, a * a], H_a, [b]; rates = [κ], order = 2)
nothing # hide

# ```math
# \begin{align}
# \frac{d}{dt} \langle a\rangle  =& -1 i \omega \langle a\rangle  + 2 i N g\Omega \left( \langle a^\dagger\rangle  + \langle a\rangle  \right) -1 i \eta \sinh\left( \xi \right) e^{1 i t {\omega}d} -1 i \eta e^{-1 i t {\omega}d} {cosh(\xi)^{*}} -0.5 \kappa \cosh\left( \xi \right) \langle a\rangle  {cosh(\xi)^{*}} + 0.5 \kappa \sinh\left( \xi \right) \langle a\rangle  {sinh(\xi)^{*}} \\
# \frac{d}{dt} \langle a^\dagger  a\rangle  =& \kappa \sinh\left( \xi \right) {sinh(\xi)^{*}} + 2 i N g\Omega \langle a^\dagger  a^\dagger\rangle  -2 i N g\Omega \langle a  a\rangle  + \kappa \sinh\left( \xi \right) {sinh(\xi)^{*}} \langle a^\dagger  a\rangle  -1 i \eta \sinh\left( \xi \right) \langle a^\dagger\rangle  e^{1 i t {\omega}d} + 1 i \eta \cosh\left( \xi \right) \langle a\rangle  e^{1 i t {\omega}d} -1 i \eta \langle a^\dagger\rangle  e^{-1 i t {\omega}d} {cosh(\xi)^{*}} + 1 i \eta \langle a\rangle  e^{-1 i t {\omega}d} {sinh(\xi)^{*}} -1.0 \kappa \cosh\left( \xi \right) {cosh(\xi)^{*}} \langle a^\dagger  a\rangle  \\
# \frac{d}{dt} \langle a  a\rangle  =& 2 i N g\Omega -2 i \omega \langle a  a\rangle  + 4 i N g\Omega \left( \langle a^\dagger  a\rangle  + \langle a  a\rangle  \right) -1.0 \kappa \sinh\left( \xi \right) {cosh(\xi)^{*}} + \kappa \sinh\left( \xi \right) {sinh(\xi)^{*}} \langle a  a\rangle  -1.0 \kappa \cosh\left( \xi \right) {cosh(\xi)^{*}} \langle a  a\rangle  -2 i \eta \sinh\left( \xi \right) \langle a\rangle  e^{1 i t {\omega}d} -2 i \eta \langle a\rangle  e^{-1 i t {\omega}d} {cosh(\xi)^{*}}
# \end{align}
# ```


sys_a = System(eqs_a; name = :sys_a) # symbolic ordinary differential equation system
sys_a = mtkcompile(sys_a)

u0_a = zeros(ComplexF64, length(eqs_a)) # initial state

gΩ_ = g_^2 / (4Ω_) # Additional parameter
N_ = 69 # the final result does not depend on N

ps_a = [ω, ωd, η, κ, N, gΩ, ξ] # symbolic parameter list (matches sys_a)
p0_a = [ω_, ωd_, η_, κ_, N_, gΩ_, ξ_]

u0_a_dict = initial_values(eqs_a, u0_a)
dict_a = merge(u0_a_dict, Dict(ps_a .=> p0_a))
prob_a = ODEProblem(sys_a, dict_a, (0.0, 4π / ωd_)) # define and solve numeric ordinary differential equation problem
sol_a = solve(prob_a, Tsit5(); saveat = π / 30ωd_, reltol = 1.0e-8, abstol = 1.0e-8)
nothing # hide

sol = sol_ls[4] # plot results
t_ = sol.t
adag_a = get_solution(sol, a' * a, eqs_sc).(t_)
aa = get_solution(sol, a * a, eqs_sc).(t_)
adag_adag = get_solution(sol, a' * a', eqs_sc).(t_)
a_ = get_solution(sol, a, eqs_sc).(t_)
adag = get_solution(sol, a', eqs_sc).(t_)
sqx = adag_adag + aa + 2 * adag_a .+ 1 - (adag + a_) .^ 2
sqy = adag_adag + aa - 2 * adag_a .- 1 - (adag - a_) .^ 2

t_a = sol_a.t
adag_a_a = get_solution(sol_a, a' * a, eqs_a).(t_a)
aa_a = get_solution(sol_a, a * a, eqs_a).(t_a)
adag_adag_a = get_solution(sol_a, a' * a', eqs_a).(t_a)
a_a = get_solution(sol_a, a, eqs_a).(t_a)
adag_a2 = get_solution(sol_a, a', eqs_a).(t_a)
sqx_a = adag_adag_a + aa_a + 2 * adag_a_a .+ 1 - (adag_a2 + a_a) .^ 2
sqy_a = adag_adag_a + aa_a - 2 * adag_a_a .- 1 - (adag_a2 - a_a) .^ 2

p = plot(xlabel = "ω t", ylabel = "Δ² O")
plot!(p, t_, real.(sqx), label = "X - Full model")
plot!(p, t_, -real.(sqy), label = "P - Full model", ls = :dash)
plot!(p, t_a, real.(sqx_a), label = "X - Effective model")
plot!(p, t_a, -real.(sqy_a), label = "P - Effective model", ls = :dash)
plot(p, size = (500, 200))

# ## Package versions

# These results were obtained using the following versions:

using InteractiveUtils
versioninfo()

using Pkg
Pkg.status(["QuantumCumulants", "OrdinaryDiffEq"], mode = PKGMODE_MANIFEST)
