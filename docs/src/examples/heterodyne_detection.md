```@meta
EditURL = "../../../examples/heterodyne_detection.jl"
```

# Heterodyne detection of emission from atomic ensemble
This example implements the heterodyne detection of the light emitted by an ensemble of atoms, described in [H. Yu, et. al., Phys. Rev. Lett 133, 073601 (2024)](https://doi.org/10.1103/PhysRevLett.133.073601). We describe the stochastic master equation time evolution that governs the measurement backaction on the system.

````@example heterodyne_detection
using QuantumCumulants
using ModelingToolkitBase
using OrdinaryDiffEqLowOrderRK
using StochasticDiffEq
using StochasticDiffEq.SciMLBase: ReturnCode
using Plots
using Random # hide
````

The system we discuss is an ensemble of $N$ equivalent two level systems with transition operators $\sigma^{\alpha \beta}_j$ and frequency $\omega_a$) coupled to a cavity mode $\hat a$, all with the same coupling $g$ and free space decay rate $\gamma$. The cavity has frequency $\omega_c$ and decay rate $\kappa$. We are in the frame rotating with the frequency $\omega_a=\omega_c$.

The system Hamiltonian is

$H=\omega_c\hat a^\dagger\hat a+\omega_a\sum_j\sigma^{22}_j+\hat a^\dagger\sigma^{12}_j+g\hat a\sigma^{21}_j.$

````@example heterodyne_detection
@variables N ωa γ η χ ωc κ g ξ ωl
@register_symbolic pulse(t)

hc = FockSpace(:resonator)
ha = NLevelSpace(:atom, 2)
h = hc ⊗ ha

j = Index(h, :j, N, ha)
k = Index(h, :k, N, ha)

@qnumbers a::Destroy(h, 1)
σ(α, β, k) = IndexedOperator(Transition(h, :σ, α, β, 2), k)

eqs_seed = meanfield(a, ωc * a' * a, [a]; rates = [κ])
t = eqs_seed.iv

H = ωc * a' * a + ωa * Σ(σ(2, 2, j), j) + g * a' * Σ(σ(1, 2, j), j) + g * a * Σ(σ(2, 1, j), j);
nothing #hide
````

We include decay of the cavity mode, where we choose the decay operator to be $\hat a\exp(\rm{i}\omega_lt)$ for convenience (we will see later why this is reasonable, and you can convince yourself that it gives the same decay term as $\hat a$). The individual decay of the atoms with rate $\gamma$ is given by the decay operator $\sigma^{12}_j$.

We also include an incoherent pump $\sigma^{21}_j$ with amplitude $\eta$, which is a pulse that is on between $t_0$ and $t_0+t_1$. The last Lindblad term is a dephasing term with strength $\chi$ corresponding to the operator $\sigma^{22}_j$.

````@example heterodyne_detection
J = [a * exp(1.0im * ωl * t), σ(1, 2, j), σ(2, 1, j), σ(2, 2, j)]
rates = [κ, γ, η * pulse(t), 2 * χ]

pulse(t) = (t > t0 && t < t0 + t1) * 1.0
nothing # hide
````

The measurement terms are included by defining measurement efficiencies for all decay channels. A channel with efficiency zero is ignored, i.e. not measured. The operator $\hat a\exp(\rm{i}\omega_lt)$ corresponding to heterodyne detection with local oscillator frequency $\omega_l$ is measured with efficiency $\xi$.

Such measurement terms are then included in the equation of motion for system as

```math
\begin{align}
d\langle\hat A\rangle=\sqrt{\xi\kappa/2}\langle \hat a^\dagger \rm{e}^{-\rm{i}\omega_lt}\hat A+\hat A\hat a\rm{e}^{\rm{i}\omega_lt} -\hat a\langle\hat a \rm{e}^{\rm{i}\omega_lt}+\hat a^\dagger \rm{e}^{-\rm{i}\omega_lt}\rangle\rangle
\end{align}
```

for the operator $\hat a \rm{e}^{\rm{i}\omega_lt}$ corresponding to heterodyne detection.

````@example heterodyne_detection
efficiencies = [ξ, 0, 0, 0]
ops = [a', a' * a, σ(2, 2, k), σ(1, 2, k), a * a]
eqs = meanfield(ops, H, J; rates = rates, efficiencies = efficiencies, direction = Forward(), order = 2, iv = t)
nothing # hide
````

```math
\begin{align}
\frac{d}{dt} \langle a^\dagger\rangle &= -0.5 \kappa \langle a^\dagger\rangle + dW/dt \left( -1 \sqrt{\xi \kappa} e^{-1.0 i {\omega}l t} \langle a^\dagger\rangle ^{2} + \sqrt{\xi \kappa} e^{1.0 i {\omega}l t} \langle a^\dagger a\rangle -1 \langle a\rangle \sqrt{\xi \kappa} e^{1.0 i {\omega}l t} \langle a^\dagger\rangle + \sqrt{\xi \kappa} e^{-1.0 i {\omega}l t} \langle a^\dagger a^\dagger\rangle \right) + 1 i \underset{j}{\overset{N}{\sum}} g \langle {\sigma}{j}^{{21}}\rangle + 1 i {\omega}c \langle a^\dagger\rangle \ \frac{d}{dt} \langle a^\dagger a\rangle &= -1 i \underset{j}{\overset{N}{\sum}} g \langle a^\dagger {\sigma}{j}^{{12}}\rangle + 1 i \underset{j}{\overset{N}{\sum}} g \langle a {\sigma}{j}^{{21}}\rangle + \left( -1 \langle a\rangle \sqrt{\xi \kappa} e^{1.0 i {\omega}l t} \langle a^\dagger a\rangle + \sqrt{\xi \kappa} e^{-1.0 i {\omega}l t} \left( -2 \langle a\rangle \langle a^\dagger\rangle ^{2} + 2 \langle a^\dagger\rangle \langle a^\dagger a\rangle + \langle a\rangle \langle a^\dagger a^\dagger\rangle \right) + \left( \langle a a\rangle \langle a^\dagger\rangle + 2 \langle a\rangle \langle a^\dagger a\rangle -2 \langle a\rangle ^{2} \langle a^\dagger\rangle \right) \sqrt{\xi \kappa} e^{1.0 i {\omega}l t} -1 \sqrt{\xi \kappa} e^{-1.0 i {\omega}l t} \langle a^\dagger\rangle \langle a^\dagger a\rangle \right) dW/dt -1.0 \kappa \langle a^\dagger a\rangle \ \frac{d}{dt} \langle {\sigma}{k}^{{22}}\rangle &= -1.0 \mathrm{pulse}\left( t \right) \eta \langle {\sigma}{k}^{{22}}\rangle + 1 i \langle a^\dagger {\sigma}{k}^{{12}}\rangle g + dW/dt \left( \langle a {\sigma}{k}^{{22}}\rangle \sqrt{\xi \kappa} e^{1.0 i {\omega}l t} -1 \langle a\rangle \sqrt{\xi \kappa} e^{1.0 i {\omega}l t} \langle {\sigma}{k}^{{22}}\rangle -1 \sqrt{\xi \kappa} e^{-1.0 i {\omega}l t} \langle {\sigma}{k}^{{22}}\rangle \langle a^\dagger\rangle + \langle a^\dagger {\sigma}{k}^{{22}}\rangle \sqrt{\xi \kappa} e^{-1.0 i {\omega}l t} \right) -1.0 \gamma \langle {\sigma}{k}^{{22}}\rangle -1 i g \langle a {\sigma}{k}^{{21}}\rangle + \mathrm{pulse}\left( t \right) \eta \ \frac{d}{dt} \langle {\sigma}{k}^{{12}}\rangle &= -0.5 \gamma \langle {\sigma}{k}^{{12}}\rangle -1 i {\omega}a \langle {\sigma}{k}^{{12}}\rangle -1 \chi \langle {\sigma}{k}^{{12}}\rangle + dW/dt \left( \sqrt{\xi \kappa} \langle a^\dagger {\sigma}{k}^{{12}}\rangle e^{-1.0 i {\omega}l t} -1 \sqrt{\xi \kappa} e^{-1.0 i {\omega}l t} \langle {\sigma}{k}^{{12}}\rangle \langle a^\dagger\rangle -1 \langle a\rangle \sqrt{\xi \kappa} \langle {\sigma}{k}^{{12}}\rangle e^{1.0 i {\omega}l t} + \langle a {\sigma}{k}^{{12}}\rangle \sqrt{\xi \kappa} e^{1.0 i {\omega}l t} \right) -1 i \langle a\rangle g -0.5 \mathrm{pulse}\left( t \right) \eta \langle {\sigma}{k}^{{12}}\rangle + 2 i \langle a {\sigma}{k}^{{22}}\rangle g \ \frac{d}{dt} \langle a a\rangle &= -1.0 \langle a a\rangle \kappa + \left( \left( 3 \langle a a\rangle \langle a\rangle -2 \langle a\rangle ^{3} \right) \sqrt{\xi \kappa} e^{1.0 i {\omega}l t} + \left( \langle a a\rangle \langle a^\dagger\rangle + 2 \langle a\rangle \langle a^\dagger a\rangle -2 \langle a\rangle ^{2} \langle a^\dagger\rangle \right) \sqrt{\xi \kappa} e^{-1.0 i {\omega}l t} -1 \langle a a\rangle \sqrt{\xi \kappa} e^{-1.0 i {\omega}l t} \langle a^\dagger\rangle -1 \langle a a\rangle \langle a\rangle \sqrt{\xi \kappa} e^{1.0 i {\omega}l t} \right) dW/dt -2 i \underset{j}{\overset{N}{\sum}} g \langle a {\sigma}_{j}^{{12}}\rangle -2 i \langle a a\rangle {\omega}c
\end{align}
```

The completion and scaling as previously discussed in other examples work exactly the same way for equations including noise terms.

````@example heterodyne_detection
eqs_c = complete(eqs; get_adjoints = false)
scaled_eqs = scale(eqs_c)
nothing # hide
````

Here we define the actual values for the system parameters. We then show that the deterministic time evolution without the noise terms is still accessible by using the constructor `System` for the stochastic system of equations and the syntax for the simulation of the time evolution is as usual.

````@example heterodyne_detection
# frequencies are in kHz
ωc_ = 0.0
κ_ = 2π * 1.13e3
ξ_ = 0.12
N_ = 5.0e4
ωa_ = 0.0
γ_ = 2π * 0.375
η_ = 2π * 20
χ_ = 0.016
g_ = 2π * 0.73
ωl_ = 2π * 1.0e3
t0 = 0.0
t1 = 20.0e-3
p = [N, ωa, γ, η, χ, ωc, κ, g, ξ, ωl]
p0 = [N_, ωa_, γ_, η_, χ_, ωc_, κ_, g_, ξ_, ωl_]
T_end = 0.1 # 0.1ms
````

To compare the stochastic evolution to the no-measurement case, build the deterministic system by passing `MeanfieldEquations(scaled_eqs)` to `System`. The constructor strips the noise drift, which also drops `ξ` from the compiled parameter list (it appeared only inside the `sqrt(ξ κ)` factor of the noise term). Because the physical parameter dict `p .=> p0` still carries `ξ`, we filter it through `parameter_map(sys, …)`; the helper keeps only the entries whose key is a live parameter or unknown of the compiled system and silently drops the rest.

````@example heterodyne_detection
sys = mtkcompile(System(MeanfieldEquations(scaled_eqs); name = :sys))
u0 = zeros(ComplexF64, length(scaled_eqs))
dict = parameter_map(sys, merge(Dict(unknowns(sys) .=> u0), Dict(p .=> p0)))
prob = ODEProblem(sys, dict, (0.0, T_end))
sol_det = solve(prob, RK4(), dt = T_end / 2.0e5)

graph = plot(xlabel = "Time (ms)", ylabel = "Cavity photon number")
plot!(graph, sol_det.t, real(get_solution(sol_det, a'a, scaled_eqs)(sol_det.t)), legend = false)
````

The stochastic time evolution is accessible via the constructor `SDESystem`, whose syntax is exactly the same as for the System, but with keyword args as defined in the [SDE tutorial](https://docs.sciml.ai/DiffEqDocs/stable/tutorials/sde_example/). We need to provide a noise process for the measurement. If the noise is white the appropriate noise process is a Wiener process. The `SDEProblem` is constructed just as the `ODEProblem`, but with an additional noise argument.

We can then make use of the `EnsembleProblem`, which automatically runs multiple instances of the stochastic equations of motion. The number of trajectories can be set in the solve call. See the tutorial linked above for more details of the function calls here.

````@example heterodyne_detection
sys_st = mtkcompile(System(scaled_eqs; name = :sys_st))
dict_st = parameter_map(sys_st, merge(Dict(unknowns(sys_st) .=> u0), Dict(p .=> p0)))

Random.seed!(2) # hide
noise = StochasticDiffEq.RealWienerProcess(0.0, 0.0)
prob_st = SDEProblem(sys_st, dict_st, (0.0, T_end); noise = noise)
sol_test = solve(prob_st, EM(); dt = T_end / 2.0e5);

plot(
    sol_test.t,
    real(get_solution(sol_test, a'a, scaled_eqs)(sol_test.t)),
    xlabel = "Time (ms)",
    ylabel = "Cavity photon number",
    legend = false,
)

Random.seed!(1) # hide
eprob = EnsembleProblem(prob_st)
traj = 200
tspan = range(0.0, T_end, length = 201)
sol = solve(
    eprob,
    StochasticDiffEq.EM(),
    dt = T_end / 2.0e5,
    save_noise = true,
    trajectories = traj,
    saveat = tspan,
)
nothing # hide
````

We plot the average of the cavity photon number for the stochastic and deterministic equation of motion with the trajectories in gray in the background. We can see that the dynamics of the system is indeed modified by the measurement backaction.

````@example heterodyne_detection
n_avg_ = zeros(length(tspan)) # average photon number
a_real_avg_ = zeros(length(tspan)) # average field (real part)
traj_succ = zeros(traj) # trajectories can be numerically unstable
for i in 1:traj
    sol_ = sol.u[i]
    if sol_.retcode == ReturnCode.Success
        n_avg_ .+= real(get_solution(sol_, a'a, scaled_eqs)(tspan))
        a_real_avg_ .+= real(get_solution(sol_, a, scaled_eqs)(tspan))
        traj_succ[i] = 1
    end
end

@show sum(traj_succ)
n_avg = n_avg_ / sum(traj_succ)
a_real_avg = a_real_avg_ / sum(traj_succ)
graph1 = plot(xlabel = "Time (ms)", ylabel = "Cavity photon number")
for i in 1:traj
    plot!(
        graph1,
        sol.u[i].t,
        real(get_solution(sol.u[i], a'a, scaled_eqs)(sol.u[i].t)),
        color = :grey,
        alpha = 0.75,
        label = nothing,
    )
end
plot!(graph1, tspan, n_avg, color = :red, label = nothing)
````

For the following cavity field amplitude we see that the ensemble average becomes zero, but the trajectories have finite cavity field amplitude.

````@example heterodyne_detection
graph2 = plot(xlabel = "Time (ms)", ylabel = "Real part cavity field")
for i in 1:traj
    plot!(graph2, tspan, real(get_solution(sol.u[i], a', scaled_eqs)(tspan)), color = :grey, alpha = 0.75, label = nothing)
end
plot!(graph2, tspan, a_real_avg, color = :red, label = nothing)
````

## Package versions

These results were obtained using the following versions:

````@example heterodyne_detection
using InteractiveUtils
versioninfo()

using Pkg
Pkg.status(
    ["QuantumCumulants", "ModelingToolkitBase", "OrdinaryDiffEq", "StochasticDiffEq", "Plots"],
    mode = PKGMODE_MANIFEST,
)
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

