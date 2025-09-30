```@meta
EditURL = "../../../examples/heterodyne_detection.jl"
```

# Heterodyne detection of emission from atomic ensemble
This example implements the heterodyne detection of the light emitted by an ensemble of atoms (see https://arxiv.org/abs/2211.13068). We thus describe the stochastic master equation time evolution that governs the measurement backaction on the system.

````@example heterodyne_detection
using QuantumCumulants
using SymbolicUtils
using Symbolics
using ModelingToolkit
using StochasticDiffEq
using DifferentialEquations
using PyPlot
using Statistics
````

## System definition
The system we discuss is an ensemble of $N$ equivalent atoms (two level systems with transition operators $\sigma_{1,2,j}$ and frequency $\omega_a$) coupled to a cavity mode $\hat a$, all with the same coupling $g$ and free space decay rate $\kappa$. The cavity has frequency $\omega_c$ and decay rate $\kappa$. We are in the frame rotating with the frequency $\omega_a=\omega_c$.

The system Hamiltonian is then
\begin{equation}
H=\omega_c\hat a^\dagger\hat a+\omega_a\sum_j\sigma_{2,2,j}+\hat a^\dagger\sigma_{1,2,j}+g\hat a\sigma_{2,1,j}.
\end{equation}

````@example heterodyne_detection
@cnumbers N ωa γ η χ ωc κ g ξ ωl
@syms t::Real
@register pulse(t)

hc = FockSpace(:resonator)
ha = NLevelSpace(:atom,2)
h = hc ⊗ ha

j = Index(h,:j,N,ha)
k = Index(h,:k,N,ha)

@qnumbers a::Destroy(h,1)
σ(α,β,k) = IndexedOperator(Transition(h,:σ,α,β,2), k)
H = ωc * a' * a + ωa * Σ(σ(2,2,j),j)+g*a'*Σ(σ(1,2,j),j)+g*a*Σ(σ(2,1,j),j);
nothing #hide
````

## Decay terms

We then include four terms in the Lindblad term of the Master equation. We include spontaneous decay of the cavity mode, where we choose the decay operator to be $\hat a\exp(\rm{i}\omega_lt)$ for convenience (we will see later why this is reasonable, and you can convince yourself that it gives the same decay term as $\hat a$). The free decay of the atoms with $\gamma$ is given by the decay operator $\sigma_{1,2,j}$.

We also include an incoherent pump $\sigma_{2,1,j}$ with amplitude $\eta$, which is a pulse that is on between $t_0$ and $t_0+t_1$. The last term is a dephasing term with strength $\chi$ corresponding to the operator $\sigma_{2,2,j}$.

````@example heterodyne_detection
J = [a*exp(1.0im*ωl*t),σ(1,2,j),σ(2,1,j),σ(2,2,j)]
rates = [κ,γ,η*pulse(t),2*χ]

function pulse(t)
    if t>t0 && t<t0+t1
        return 1
    else
        return 0
    end
end
````

## Stochastic measurement terms
We include the measurement terms by defining measurement efficiencies for all decay channels. A channel with efficiency zero is ignored, i.e. not measured. The operator $\hat a\exp(\rm{i}\omega_lt)$ corresponding to heterodyne detection with local oscillator frequency $\omega_l$ is then measured with efficiency $\xi$.

Such measurement terms are then included in the equation of motion for system operators is then included as

\begin{equation}
d\langle\hat A\rangle=\sqrt{\xi\kappa/2}\langle \hat a^\dagger \rm{e}^{-\rm{i}\omega_lt}\hat A+\hat A\hat a\rm{e}^{\rm{i}\omega_lt} -\hat a\langle\hat a \rm{e}^{\rm{i}\omega_lt}+\hat a^\dagger \rm{e}^{-\rm{i}\omega_lt}\rangle\rangle
\end{equation}

for the operator $\hat a \rm{e}^{\rm{i}\omega_lt}$  corresponding to heterodyne detection.

````@example heterodyne_detection
efficiencies = [ξ,0,0,0]
ops = [a',a'*a,σ(2,2,k),σ(1,2,k), a*a];
eqs = meanfield(ops,H,J; rates = rates, efficiencies = efficiencies, order = 2)
````

## Finding the equations of motion
The completion and scaling as previously discussed in other examples work exactly the same way for equations including noise terms.

````@example heterodyne_detection
eqs_c = indexed_complete(eqs)
scaled_eqs = scale(eqs_c)
````

## Deterministic time evolution
Here we define the actual values for the system parameters. We then show that the deterministic time evolution without the noise terms is still accessible by using the constructor System for the stochastic system of equations and the syntax for the simulation of the time evolution is as usual.

````@example heterodyne_detection
ωc_ = 0.0; κ_ = 2.0 * π * 2.26e6; ξ_ = 0.12; N_=5e4;
ωa_ = 0.0; γ_ = 2.0 * π * 375; η_ = 2.0 * π * 20e3; χ_ = 0.1;
g_ = 6.531*10^3; ωl_ = 2.0 * π * 10^3; t0=0.0;t1=20e-6
p = [N,ωa,γ,η,χ,ωc,κ,g,ξ,ωl]
p0 = [N_,ωa_,γ_,η_,χ_,ωc_,κ_,g_,ξ_,ωl_]

@named sys = System(scaled_eqs)
u0 = zeros(ComplexF64,length(scaled_eqs.equations))
dict = merge(Dict(unknowns(sys) .=> u0), Dict(p .=> p0))
prob = ODEProblem(sys, dict,(0.0,1e-3))
sol_det = solve(prob,RK4(),dt=1e-9)

plot(sol_det.t .* 1000, map(x -> real(x[2]), sol_det.u))
xlim(left = 0.0, right = 0.1)
ylabel("Cavity phonon number")
xlabel("Time (ms)")
````

## Stochastic time evolution
The stochastic time evolution is accessible via the constructor SDESystem, whose syntax is exactly the same as for the System, but with keyword args as defined in https://docs.sciml.ai/DiffEqDocs/stable/tutorials/sde_example/. We then need to provide a noise process for the measurement. If the noise is white the appropriate noise process is a Wiener process. The SDEProblem is then constructed just as the ODEProblem, but with an additional noise argument.

We can then make use of the EnsembleProblem which automatically runs multiple instances of the stochastic equations of motion. The number of trajectories can then be set in the solve call. See the documentation cited above for more details of the function calls here.

````@example heterodyne_detection
@named sys = SDESystem(scaled_eqs)
u0 = zeros(ComplexF64,length(scaled_eqs.equations))
noise = StochasticDiffEq.RealWienerProcess(0.0, 0.0)
prob = SDEProblem(sys, u0,(0.0,0.1*1e-3),p.=>p0,noise=noise)

eprob = EnsembleProblem(prob)
sol = solve(eprob,StochasticDiffEq.EM(), dt = 1e-9, save_noise=true,trajectories=256,saveat=0:1e-6:0.1*1e-3)
````

## Plot the average of the cavity number
Here we plot the average of the cavity number for the stochastic and deterministic equation of motion with the trajectories in gray in the background. We can see that the dynamics of the system is indeed modified by the measurement backaction.

````@example heterodyne_detection
us = fill(0.0, length(sol.u), length(sol.u[1].t))
for (i, el) in enumerate(sol.u)
    u = map(x->real(x[2]), el)
    if length(u) == size(us, 2) us[i,:] .+= u else continue end
    plot(el.t .* 1000, u, color = "grey", alpha = 0.1)
end

plot(sol.u[1].t .* 1000, mean(us, dims = 1)[1,:], color = "red", label = "stochastic")
plot(sol_det.t .* 1000, map(x -> real(x[2]), sol_det.u), color = "blue", label = "deterministic")
xlim(left = 0.0, right = 0.1)
ylim(bottom = 0, top = 500)
ylabel("Cavity phonon number")
xlabel("Time (ms)")
legend()
````

## Plot the average of the cavity field amplitude
For the cavity field amplitude we see that the ensemble average becomes zero, but the trajectories have finite cavity field amplitude

````@example heterodyne_detection
us = fill(0.0, length(sol.u), length(sol.u[1].t))
for (i, el) in enumerate(sol.u)
    u = map(x->real(x[1]), el)
    if length(u) == size(us, 2) us[i,:] .+= u else continue end
    plot(el.t .* 1000, u, color = "grey", alpha = 0.05)
end

plot(sol_det.t .* 1000, map(x -> real(x[1]), sol_det.u), color = "blue", label = "deterministic")
plot(sol.u[1].t .* 1000, mean(us, dims = 1)[1,:], color = "red", label = "stochastic")
xlim(left = 0.0, right = 0.1)
ylabel("Cavity field")
xlabel("Time (ms)")
legend()
````

## Package versions

These results were obtained using the following versions:

````@example heterodyne_detection
using InteractiveUtils
versioninfo()

using Pkg
Pkg.status(["SummationByPartsOperators", "OrdinaryDiffEq"],
           mode=PKGMODE_MANIFEST)
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

