# # Heterodyne detection of emission from atomic ensemble
# This example implements the heterodyne detection of the light emitted by an ensemble of atoms, described in [H. Yu, et. al., Phys. Rev. Lett 133, 073601 (2024)](https://doi.org/10.1103/PhysRevLett.133.073601). We describe the stochastic master equation time evolution that governs the measurement backaction on the system.

using QuantumCumulants
using SymbolicUtils
using Symbolics
using ModelingToolkit
using StochasticDiffEq
using OrdinaryDiffEq
using Plots
using Statistics
using Random # hide

# The system we discuss is an ensemble of $N$ equivalent two level systems with transition operators $\sigma^{\alpha \beta}_j$ and frequency $\omega_a$) coupled to a cavity mode $\hat a$, all with the same coupling $g$ and free space decay rate $\gamma$. The cavity has frequency $\omega_c$ and decay rate $\kappa$. We are in the frame rotating with the frequency $\omega_a=\omega_c$.

# The system Hamiltonian is 

# $H=\omega_c\hat a^\dagger\hat a+\omega_a\sum_j\sigma^{22}_j+\hat a^\dagger\sigma^{12}_j+g\hat a\sigma^{21}_j$.

@rnumbers N ωa γ η χ ωc κ g ξ ωl
@syms t::Real
@register_symbolic pulse(t)

hc = FockSpace(:resonator)
ha = NLevelSpace(:atom, 2)
h = hc ⊗ ha

j = Index(h, :j, N, ha)
k = Index(h, :k, N, ha)

@qnumbers a::Destroy(h, 1)
σ(α, β, k) = IndexedOperator(Transition(h, :σ, α, β, 2), k)
H = ωc * a' * a + ωa * Σ(σ(2, 2, j), j) + g*a'*Σ(σ(1, 2, j), j) + g*a*Σ(σ(2, 1, j), j);

# We include decay of the cavity mode, where we choose the decay operator to be $\hat a\exp(\rm{i}\omega_lt)$ for convenience (we will see later why this is reasonable, and you can convince yourself that it gives the same decay term as $\hat a$). The individual decay of the atoms with rate $\gamma$ is given by the decay operator $\sigma^{12}_j$.

# We also include an incoherent pump $\sigma^{21}_j$ with amplitude $\eta$, which is a pulse that is on between $t_0$ and $t_0+t_1$. The last Lindblad term is a dephasing term with strength $\chi$ corresponding to the operator $\sigma^{22}_j$.

J = [a*exp(1.0im*ωl*t), σ(1, 2, j), σ(2, 1, j), σ(2, 2, j)]
rates = [κ, γ, η*pulse(t), 2*χ]

function pulse(t)
    if t>t0 && t<t0+t1
        return 1
    else
        return 0
    end
end
nothing # hide

# The measurement terms are included by defining measurement efficiencies for all decay channels. A channel with efficiency zero is ignored, i.e. not measured. The operator $\hat a\exp(\rm{i}\omega_lt)$ corresponding to heterodyne detection with local oscillator frequency $\omega_l$ is measured with efficiency $\xi$.

# Such measurement terms are then included in the equation of motion for system as

# ```math
# \begin{align}
# d\langle\hat A\rangle=\sqrt{\xi\kappa/2}\langle \hat a^\dagger \rm{e}^{-\rm{i}\omega_lt}\hat A+\hat A\hat a\rm{e}^{\rm{i}\omega_lt} -\hat a\langle\hat a \rm{e}^{\rm{i}\omega_lt}+\hat a^\dagger \rm{e}^{-\rm{i}\omega_lt}\rangle\rangle
# \end{align}
# ```

# for the operator $\hat a \rm{e}^{\rm{i}\omega_lt}$ corresponding to heterodyne detection.

efficiencies = [ξ, 0, 0, 0]
ops = [a', a'*a, σ(2, 2, k), σ(1, 2, k), a*a]
eqs = meanfield(ops, H, J; rates = rates, efficiencies = efficiencies, order = 2)
nothing # hide

# ```math
# \begin{align} 
# \frac{d}{dt} \langle a^\dagger\rangle &= -0.5 \kappa \langle a^\dagger\rangle + dW/dt \left( -1 \sqrt{\xi \kappa} e^{-1.0 i {\omega}l t} \langle a^\dagger\rangle ^{2} + \sqrt{\xi \kappa} e^{1.0 i {\omega}l t} \langle a^\dagger a\rangle -1 \langle a\rangle \sqrt{\xi \kappa} e^{1.0 i {\omega}l t} \langle a^\dagger\rangle + \sqrt{\xi \kappa} e^{-1.0 i {\omega}l t} \langle a^\dagger a^\dagger\rangle \right) + 1 i \underset{j}{\overset{N}{\sum}} g \langle {\sigma}{j}^{{21}}\rangle + 1 i {\omega}c \langle a^\dagger\rangle \ \frac{d}{dt} \langle a^\dagger a\rangle &= -1 i \underset{j}{\overset{N}{\sum}} g \langle a^\dagger {\sigma}{j}^{{12}}\rangle + 1 i \underset{j}{\overset{N}{\sum}} g \langle a {\sigma}{j}^{{21}}\rangle + \left( -1 \langle a\rangle \sqrt{\xi \kappa} e^{1.0 i {\omega}l t} \langle a^\dagger a\rangle + \sqrt{\xi \kappa} e^{-1.0 i {\omega}l t} \left( -2 \langle a\rangle \langle a^\dagger\rangle ^{2} + 2 \langle a^\dagger\rangle \langle a^\dagger a\rangle + \langle a\rangle \langle a^\dagger a^\dagger\rangle \right) + \left( \langle a a\rangle \langle a^\dagger\rangle + 2 \langle a\rangle \langle a^\dagger a\rangle -2 \langle a\rangle ^{2} \langle a^\dagger\rangle \right) \sqrt{\xi \kappa} e^{1.0 i {\omega}l t} -1 \sqrt{\xi \kappa} e^{-1.0 i {\omega}l t} \langle a^\dagger\rangle \langle a^\dagger a\rangle \right) dW/dt -1.0 \kappa \langle a^\dagger a\rangle \ \frac{d}{dt} \langle {\sigma}{k}^{{22}}\rangle &= -1.0 \mathrm{pulse}\left( t \right) \eta \langle {\sigma}{k}^{{22}}\rangle + 1 i \langle a^\dagger {\sigma}{k}^{{12}}\rangle g + dW/dt \left( \langle a {\sigma}{k}^{{22}}\rangle \sqrt{\xi \kappa} e^{1.0 i {\omega}l t} -1 \langle a\rangle \sqrt{\xi \kappa} e^{1.0 i {\omega}l t} \langle {\sigma}{k}^{{22}}\rangle -1 \sqrt{\xi \kappa} e^{-1.0 i {\omega}l t} \langle {\sigma}{k}^{{22}}\rangle \langle a^\dagger\rangle + \langle a^\dagger {\sigma}{k}^{{22}}\rangle \sqrt{\xi \kappa} e^{-1.0 i {\omega}l t} \right) -1.0 \gamma \langle {\sigma}{k}^{{22}}\rangle -1 i g \langle a {\sigma}{k}^{{21}}\rangle + \mathrm{pulse}\left( t \right) \eta \ \frac{d}{dt} \langle {\sigma}{k}^{{12}}\rangle &= -0.5 \gamma \langle {\sigma}{k}^{{12}}\rangle -1 i {\omega}a \langle {\sigma}{k}^{{12}}\rangle -1 \chi \langle {\sigma}{k}^{{12}}\rangle + dW/dt \left( \sqrt{\xi \kappa} \langle a^\dagger {\sigma}{k}^{{12}}\rangle e^{-1.0 i {\omega}l t} -1 \sqrt{\xi \kappa} e^{-1.0 i {\omega}l t} \langle {\sigma}{k}^{{12}}\rangle \langle a^\dagger\rangle -1 \langle a\rangle \sqrt{\xi \kappa} \langle {\sigma}{k}^{{12}}\rangle e^{1.0 i {\omega}l t} + \langle a {\sigma}{k}^{{12}}\rangle \sqrt{\xi \kappa} e^{1.0 i {\omega}l t} \right) -1 i \langle a\rangle g -0.5 \mathrm{pulse}\left( t \right) \eta \langle {\sigma}{k}^{{12}}\rangle + 2 i \langle a {\sigma}{k}^{{22}}\rangle g \ \frac{d}{dt} \langle a a\rangle &= -1.0 \langle a a\rangle \kappa + \left( \left( 3 \langle a a\rangle \langle a\rangle -2 \langle a\rangle ^{3} \right) \sqrt{\xi \kappa} e^{1.0 i {\omega}l t} + \left( \langle a a\rangle \langle a^\dagger\rangle + 2 \langle a\rangle \langle a^\dagger a\rangle -2 \langle a\rangle ^{2} \langle a^\dagger\rangle \right) \sqrt{\xi \kappa} e^{-1.0 i {\omega}l t} -1 \langle a a\rangle \sqrt{\xi \kappa} e^{-1.0 i {\omega}l t} \langle a^\dagger\rangle -1 \langle a a\rangle \langle a\rangle \sqrt{\xi \kappa} e^{1.0 i {\omega}l t} \right) dW/dt -2 i \underset{j}{\overset{N}{\sum}} g \langle a {\sigma}_{j}^{{12}}\rangle -2 i \langle a a\rangle {\omega}c 
# \end{align}
# ```

# The completion and scaling as previously discussed in other examples work exactly the same way for equations including noise terms.

eqs_c = indexed_complete(eqs)
scaled_eqs = scale(eqs_c)
nothing # hide

# ```math
# \begin{align} 
# \frac{d}{dt} \langle a^\dagger\rangle &= -0.5 \kappa \langle a^\dagger\rangle + dW/dt \left( -1 \sqrt{\xi \kappa} e^{-1.0 i {\omega}l t} \langle a^\dagger\rangle ^{2} + \sqrt{\xi \kappa} e^{1.0 i {\omega}l t} \langle a^\dagger a\rangle -1 \langle a\rangle \sqrt{\xi \kappa} e^{1.0 i {\omega}l t} \langle a^\dagger\rangle + \sqrt{\xi \kappa} e^{-1.0 i {\omega}l t} \langle a^\dagger a^\dagger\rangle \right) + 1 i \langle {\sigma}{1}^{{21}}\rangle g N + 1 i {\omega}c \langle a^\dagger\rangle \ \frac{d}{dt} \langle a^\dagger a\rangle &= -1 i \langle a^\dagger {\sigma}{1}^{{12}}\rangle g N + \left( -1 \langle a\rangle \sqrt{\xi \kappa} e^{1.0 i {\omega}l t} \langle a^\dagger a\rangle + \sqrt{\xi \kappa} e^{-1.0 i {\omega}l t} \left( -2 \langle a\rangle \langle a^\dagger\rangle ^{2} + 2 \langle a^\dagger\rangle \langle a^\dagger a\rangle + \langle a\rangle \langle a^\dagger a^\dagger\rangle \right) + \left( \langle a a\rangle \langle a^\dagger\rangle + 2 \langle a\rangle \langle a^\dagger a\rangle -2 \langle a\rangle ^{2} \langle a^\dagger\rangle \right) \sqrt{\xi \kappa} e^{1.0 i {\omega}l t} -1 \sqrt{\xi \kappa} e^{-1.0 i {\omega}l t} \langle a^\dagger\rangle \langle a^\dagger a\rangle \right) dW/dt + 1 i g N \langle a {\sigma}{1}^{{21}}\rangle -1.0 \kappa \langle a^\dagger a\rangle \ \frac{d}{dt} \langle {\sigma}{1}^{{22}}\rangle &= -1.0 \mathrm{pulse}\left( t \right) \langle {\sigma}{1}^{{22}}\rangle \eta -1 i g \langle a {\sigma}{1}^{{21}}\rangle + dW/dt \left( \langle a {\sigma}{1}^{{22}}\rangle \sqrt{\xi \kappa} e^{1.0 i {\omega}l t} -1 \langle {\sigma}{1}^{{22}}\rangle \sqrt{\xi \kappa} e^{-1.0 i {\omega}l t} \langle a^\dagger\rangle + \langle a^\dagger {\sigma}{1}^{{22}}\rangle \sqrt{\xi \kappa} e^{-1.0 i {\omega}l t} -1 \langle {\sigma}{1}^{{22}}\rangle \langle a\rangle \sqrt{\xi \kappa} e^{1.0 i {\omega}l t} \right) -1.0 \langle {\sigma}{1}^{{22}}\rangle \gamma + 1 i \langle a^\dagger {\sigma}{1}^{{12}}\rangle g + \mathrm{pulse}\left( t \right) \eta \ \frac{d}{dt} \langle {\sigma}{1}^{{12}}\rangle &= 2 i \langle a {\sigma}{1}^{{22}}\rangle g + dW/dt \left( \sqrt{\xi \kappa} \langle a {\sigma}{1}^{{12}}\rangle e^{1.0 i {\omega}l t} -1 \langle a\rangle \sqrt{\xi \kappa} \langle {\sigma}{1}^{{12}}\rangle e^{1.0 i {\omega}l t} -1 \sqrt{\xi \kappa} \langle {\sigma}{1}^{{12}}\rangle e^{-1.0 i {\omega}l t} \langle a^\dagger\rangle + \langle a^\dagger {\sigma}{1}^{{12}}\rangle \sqrt{\xi \kappa} e^{-1.0 i {\omega}l t} \right) -1 i {\omega}a \langle {\sigma}{1}^{{12}}\rangle -0.5 \mathrm{pulse}\left( t \right) \langle {\sigma}{1}^{{12}}\rangle \eta -1 i \langle a\rangle g -0.5 \gamma \langle {\sigma}{1}^{{12}}\rangle -1 \chi \langle {\sigma}{1}^{{12}}\rangle \ \frac{d}{dt} \langle a a\rangle &= -1.0 \langle a a\rangle \kappa + \left( \left( 3 \langle a a\rangle \langle a\rangle -2 \langle a\rangle ^{3} \right) \sqrt{\xi \kappa} e^{1.0 i {\omega}l t} + \left( \langle a a\rangle \langle a^\dagger\rangle + 2 \langle a\rangle \langle a^\dagger a\rangle -2 \langle a\rangle ^{2} \langle a^\dagger\rangle \right) \sqrt{\xi \kappa} e^{-1.0 i {\omega}l t} -1 \langle a a\rangle \sqrt{\xi \kappa} e^{-1.0 i {\omega}l t} \langle a^\dagger\rangle -1 \langle a a\rangle \langle a\rangle \sqrt{\xi \kappa} e^{1.0 i {\omega}l t} \right) dW/dt -2 i g \langle a {\sigma}{1}^{{12}}\rangle N -2 i \langle a a\rangle {\omega}c \ \frac{d}{dt} \langle a^\dagger {\sigma}{1}^{{12}}\rangle &= -1 i g \langle a^\dagger a\rangle -0.5 \langle a^\dagger {\sigma}{1}^{{12}}\rangle \mathrm{pulse}\left( t \right) \eta + 2 i g \left( \langle {\sigma}{1}^{{22}}\rangle \langle a^\dagger a\rangle + \langle a {\sigma}{1}^{{22}}\rangle \langle a^\dagger\rangle -2 \langle {\sigma}{1}^{{22}}\rangle \langle a\rangle \langle a^\dagger\rangle + \langle a^\dagger {\sigma}{1}^{{22}}\rangle \langle a\rangle \right) + 1 i \left( -1 + N \right) g \langle {\sigma}{1}^{{21}} {\sigma}{2}^{{12}}\rangle -1 i {\omega}a \langle a^\dagger {\sigma}{1}^{{12}}\rangle + 1 i \langle {\sigma}{1}^{{22}}\rangle g -0.5 \langle a^\dagger {\sigma}{1}^{{12}}\rangle \left( \gamma + \kappa \right) + 1 i \langle a^\dagger {\sigma}{1}^{{12}}\rangle {\omega}c -1 \langle a^\dagger {\sigma}{1}^{{12}}\rangle \chi + dW/dt \left( -1 \langle a^\dagger {\sigma}{1}^{{12}}\rangle \langle a\rangle \sqrt{\xi \kappa} e^{1.0 i {\omega}l t} + \sqrt{\xi \kappa} e^{-1.0 i {\omega}l t} \left( -2 \langle {\sigma}{1}^{{12}}\rangle \langle a^\dagger\rangle ^{2} + \langle {\sigma}{1}^{{12}}\rangle \langle a^\dagger a^\dagger\rangle + 2 \langle a^\dagger {\sigma}{1}^{{12}}\rangle \langle a^\dagger\rangle \right) + \sqrt{\xi \kappa} e^{1.0 i {\omega}l t} \left( -2 \langle a\rangle \langle {\sigma}{1}^{{12}}\rangle \langle a^\dagger\rangle + \langle {\sigma}{1}^{{12}}\rangle \langle a^\dagger a\rangle + \langle a^\dagger {\sigma}{1}^{{12}}\rangle \langle a\rangle + \langle a {\sigma}{1}^{{12}}\rangle \langle a^\dagger\rangle \right) -1 \langle a^\dagger {\sigma}{1}^{{12}}\rangle \sqrt{\xi \kappa} e^{-1.0 i {\omega}l t} \langle a^\dagger\rangle \right) \ \frac{d}{dt} \langle a {\sigma}{1}^{{22}}\rangle &= -1 i \left( -1 + N \right) g \langle {\sigma}{1}^{{12}} {\sigma}{2}^{{22}}\rangle -1 i \langle a {\sigma}{1}^{{22}}\rangle {\omega}c + dW/dt \left( \sqrt{\xi \kappa} e^{-1.0 i {\omega}l t} \left( \langle {\sigma}{1}^{{22}}\rangle \langle a^\dagger a\rangle + \langle a {\sigma}{1}^{{22}}\rangle \langle a^\dagger\rangle -2 \langle {\sigma}{1}^{{22}}\rangle \langle a\rangle \langle a^\dagger\rangle + \langle a^\dagger {\sigma}{1}^{{22}}\rangle \langle a\rangle \right) -1 \langle a {\sigma}{1}^{{22}}\rangle \sqrt{\xi \kappa} e^{-1.0 i {\omega}l t} \langle a^\dagger\rangle -1 \langle a {\sigma}{1}^{{22}}\rangle \langle a\rangle \sqrt{\xi \kappa} e^{1.0 i {\omega}l t} + \left( -2 \langle {\sigma}{1}^{{22}}\rangle \langle a\rangle ^{2} + \langle a a\rangle \langle {\sigma}{1}^{{22}}\rangle + 2 \langle a {\sigma}{1}^{{22}}\rangle \langle a\rangle \right) \sqrt{\xi \kappa} e^{1.0 i {\omega}l t} \right) -1 i \left( \langle a a\rangle \langle {\sigma}{1}^{{21}}\rangle + 2 \langle a\rangle \langle a {\sigma}{1}^{{21}}\rangle -2 \langle {\sigma}{1}^{{21}}\rangle \langle a\rangle ^{2} \right) g + 1 i g \left( -2 \langle a\rangle \langle {\sigma}{1}^{{12}}\rangle \langle a^\dagger\rangle + \langle {\sigma}{1}^{{12}}\rangle \langle a^\dagger a\rangle + \langle a^\dagger {\sigma}{1}^{{12}}\rangle \langle a\rangle + \langle a {\sigma}{1}^{{12}}\rangle \langle a^\dagger\rangle \right) + \mathrm{pulse}\left( t \right) \langle a\rangle \eta -0.5 \langle a {\sigma}{1}^{{22}}\rangle \kappa -1.0 \langle a {\sigma}{1}^{{22}}\rangle \gamma -1.0 \mathrm{pulse}\left( t \right) \langle a {\sigma}{1}^{{22}}\rangle \eta \ \frac{d}{dt} \langle a {\sigma}{1}^{{12}}\rangle &= -0.5 \mathrm{pulse}\left( t \right) \langle a {\sigma}{1}^{{12}}\rangle \eta -0.5 \langle a {\sigma}{1}^{{12}}\rangle \left( \gamma + \kappa \right) -1 i \left( {\omega}a + {\omega}c \right) \langle a {\sigma}{1}^{{12}}\rangle -1 i \langle a a\rangle g + dW/dt \left( -1 \langle a\rangle \sqrt{\xi \kappa} \langle a {\sigma}{1}^{{12}}\rangle e^{1.0 i {\omega}l t} + \left( -2 \langle a\rangle ^{2} \langle {\sigma}{1}^{{12}}\rangle + \langle a a\rangle \langle {\sigma}{1}^{{12}}\rangle + 2 \langle a\rangle \langle a {\sigma}{1}^{{12}}\rangle \right) \sqrt{\xi \kappa} e^{1.0 i {\omega}l t} + \sqrt{\xi \kappa} e^{-1.0 i {\omega}l t} \left( -2 \langle a\rangle \langle {\sigma}{1}^{{12}}\rangle \langle a^\dagger\rangle + \langle {\sigma}{1}^{{12}}\rangle \langle a^\dagger a\rangle + \langle a^\dagger {\sigma}{1}^{{12}}\rangle \langle a\rangle + \langle a {\sigma}{1}^{{12}}\rangle \langle a^\dagger\rangle \right) -1 \sqrt{\xi \kappa} \langle a {\sigma}{1}^{{12}}\rangle e^{-1.0 i {\omega}l t} \langle a^\dagger\rangle \right) -1 \chi \langle a {\sigma}{1}^{{12}}\rangle -1 i \left( -1 + N \right) g \langle {\sigma}{1}^{{12}} {\sigma}{2}^{{12}}\rangle + 2 i \left( -2 \langle {\sigma}{1}^{{22}}\rangle \langle a\rangle ^{2} + \langle a a\rangle \langle {\sigma}{1}^{{22}}\rangle + 2 \langle a {\sigma}{1}^{{22}}\rangle \langle a\rangle \right) g \ \frac{d}{dt} \langle {\sigma}{1}^{{21}} {\sigma}{2}^{{12}}\rangle &= -1.0 \gamma \langle {\sigma}{1}^{{21}} {\sigma}{2}^{{12}}\rangle -1.0 \mathrm{pulse}\left( t \right) \langle {\sigma}{1}^{{21}} {\sigma}{2}^{{12}}\rangle \eta + dW/dt \left( \sqrt{\xi \kappa} e^{-1.0 i {\omega}l t} \left( \langle a^\dagger {\sigma}{1}^{{12}}\rangle \langle {\sigma}{1}^{{21}}\rangle -2 \langle {\sigma}{1}^{{21}}\rangle \langle {\sigma}{1}^{{12}}\rangle \langle a^\dagger\rangle + \langle {\sigma}{1}^{{21}} {\sigma}{2}^{{12}}\rangle \langle a^\dagger\rangle + \langle a^\dagger {\sigma}{1}^{{21}}\rangle \langle {\sigma}{1}^{{12}}\rangle \right) -1 \sqrt{\xi \kappa} \langle {\sigma}{1}^{{21}} {\sigma}{2}^{{12}}\rangle e^{-1.0 i {\omega}l t} \langle a^\dagger\rangle -1 \sqrt{\xi \kappa} \langle a\rangle \langle {\sigma}{1}^{{21}} {\sigma}{2}^{{12}}\rangle e^{1.0 i {\omega}l t} + \sqrt{\xi \kappa} \left( \langle {\sigma}{1}^{{12}}\rangle \langle a {\sigma}{1}^{{21}}\rangle -2 \langle {\sigma}{1}^{{21}}\rangle \langle a\rangle \langle {\sigma}{1}^{{12}}\rangle + \langle {\sigma}{1}^{{21}}\rangle \langle a {\sigma}{1}^{{12}}\rangle + \langle a\rangle \langle {\sigma}{1}^{{21}} {\sigma}{2}^{{12}}\rangle \right) e^{1.0 i {\omega}l t} \right) -1 i g \langle a {\sigma}{1}^{{21}}\rangle + 2 i \left( \langle {\sigma}{1}^{{22}}\rangle \langle a {\sigma}{1}^{{21}}\rangle + \langle {\sigma}{1}^{{21}}\rangle \langle a {\sigma}{1}^{{22}}\rangle -2 \langle {\sigma}{1}^{{21}}\rangle \langle {\sigma}{1}^{{22}}\rangle \langle a\rangle + {\langle {\sigma}{1}^{{12}} {\sigma}{2}^{{22}}\rangle ^{}} \langle a\rangle \right) g -2 i g \left( \langle a^\dagger {\sigma}{1}^{{12}}\rangle \langle {\sigma}{1}^{{22}}\rangle + \langle {\sigma}{1}^{{12}} {\sigma}{2}^{{22}}\rangle \langle a^\dagger\rangle -2 \langle {\sigma}{1}^{{22}}\rangle \langle {\sigma}{1}^{{12}}\rangle \langle a^\dagger\rangle + \langle a^\dagger {\sigma}{1}^{{22}}\rangle \langle {\sigma}{1}^{{12}}\rangle \right) -2 \chi \langle {\sigma}{1}^{{21}} {\sigma}{2}^{{12}}\rangle + 1 i \langle a^\dagger {\sigma}{1}^{{12}}\rangle g \ \frac{d}{dt} \langle {\sigma}{1}^{{12}} {\sigma}{2}^{{22}}\rangle &= -1 i {\omega}a \langle {\sigma}{1}^{{12}} {\sigma}{2}^{{22}}\rangle -1.5 \gamma \langle {\sigma}{1}^{{12}} {\sigma}{2}^{{22}}\rangle -1.5 \mathrm{pulse}\left( t \right) \eta \langle {\sigma}{1}^{{12}} {\sigma}{2}^{{22}}\rangle -1 i g \left( \langle a {\sigma}{1}^{{22}}\rangle + \langle {\sigma}{1}^{{12}}\rangle \langle a {\sigma}{1}^{{21}}\rangle -2 \langle {\sigma}{1}^{{21}}\rangle \langle a\rangle \langle {\sigma}{1}^{{12}}\rangle + \langle {\sigma}{1}^{{21}}\rangle \langle a {\sigma}{1}^{{12}}\rangle + \langle a\rangle \langle {\sigma}{1}^{{21}} {\sigma}{2}^{{12}}\rangle \right) + \mathrm{pulse}\left( t \right) \langle {\sigma}{1}^{{12}}\rangle \eta -1 \chi \langle {\sigma}{1}^{{12}} {\sigma}{2}^{{22}}\rangle + 1 i \left( 2 \langle a^\dagger {\sigma}{1}^{{12}}\rangle \langle {\sigma}{1}^{{12}}\rangle + \langle {\sigma}{1}^{{12}} {\sigma}{2}^{{12}}\rangle \langle a^\dagger\rangle -2 \langle {\sigma}{1}^{{12}}\rangle ^{2} \langle a^\dagger\rangle \right) g + dW/dt \left( \sqrt{\xi \kappa} e^{1.0 i {\omega}l t} \left( \langle {\sigma}{1}^{{22}}\rangle \langle a {\sigma}{1}^{{12}}\rangle + \langle a\rangle \langle {\sigma}{1}^{{12}} {\sigma}{2}^{{22}}\rangle + \langle a {\sigma}{1}^{{22}}\rangle \langle {\sigma}{1}^{{12}}\rangle -2 \langle {\sigma}{1}^{{22}}\rangle \langle a\rangle \langle {\sigma}{1}^{{12}}\rangle \right) -1 \sqrt{\xi \kappa} \langle a\rangle \langle {\sigma}{1}^{{12}} {\sigma}{2}^{{22}}\rangle e^{1.0 i {\omega}l t} + \sqrt{\xi \kappa} e^{-1.0 i {\omega}l t} \left( \langle a^\dagger {\sigma}{1}^{{12}}\rangle \langle {\sigma}{1}^{{22}}\rangle + \langle {\sigma}{1}^{{12}} {\sigma}{2}^{{22}}\rangle \langle a^\dagger\rangle -2 \langle {\sigma}{1}^{{22}}\rangle \langle {\sigma}{1}^{{12}}\rangle \langle a^\dagger\rangle + \langle a^\dagger {\sigma}{1}^{{22}}\rangle \langle {\sigma}{1}^{{12}}\rangle \right) -1 \sqrt{\xi \kappa} e^{-1.0 i {\omega}l t} \langle {\sigma}{1}^{{12}} {\sigma}{2}^{{22}}\rangle \langle a^\dagger\rangle \right) + 2 i \left( \langle {\sigma}{1}^{{22}} {\sigma}{2}^{{22}}\rangle \langle a\rangle + 2 \langle {\sigma}{1}^{{22}}\rangle \langle a {\sigma}{1}^{{22}}\rangle -2 \langle {\sigma}{1}^{{22}}\rangle ^{2} \langle a\rangle \right) g \ \frac{d}{dt} \langle {\sigma}{1}^{{12}} {\sigma}{2}^{{12}}\rangle &= -1.0 \mathrm{pulse}\left( t \right) \eta \langle {\sigma}{1}^{{12}} {\sigma}{2}^{{12}}\rangle -1.0 \gamma \langle {\sigma}{1}^{{12}} {\sigma}{2}^{{12}}\rangle -2 i {\omega}a \langle {\sigma}{1}^{{12}} {\sigma}{2}^{{12}}\rangle -2 i g \langle a {\sigma}{1}^{{12}}\rangle -2 \chi \langle {\sigma}{1}^{{12}} {\sigma}{2}^{{12}}\rangle + 2 i g \left( 2 \langle {\sigma}{1}^{{22}}\rangle \langle a {\sigma}{1}^{{12}}\rangle + 2 \langle a\rangle \langle {\sigma}{1}^{{12}} {\sigma}{2}^{{22}}\rangle + 2 \langle a {\sigma}{1}^{{22}}\rangle \langle {\sigma}{1}^{{12}}\rangle -4 \langle {\sigma}{1}^{{22}}\rangle \langle a\rangle \langle {\sigma}{1}^{{12}}\rangle \right) + dW/dt \left( -1 \sqrt{\xi \kappa} e^{-1.0 i {\omega}l t} \langle {\sigma}{1}^{{12}} {\sigma}{2}^{{12}}\rangle \langle a^\dagger\rangle + \left( 2 \langle a^\dagger {\sigma}{1}^{{12}}\rangle \langle {\sigma}{1}^{{12}}\rangle + \langle {\sigma}{1}^{{12}} {\sigma}{2}^{{12}}\rangle \langle a^\dagger\rangle -2 \langle {\sigma}{1}^{{12}}\rangle ^{2} \langle a^\dagger\rangle \right) \sqrt{\xi \kappa} e^{-1.0 i {\omega}l t} + \sqrt{\xi \kappa} \left( 2 \langle {\sigma}{1}^{{12}}\rangle \langle a {\sigma}{1}^{{12}}\rangle + \langle a\rangle \langle {\sigma}{1}^{{12}} {\sigma}{2}^{{12}}\rangle -2 \langle a\rangle \langle {\sigma}{1}^{{12}}\rangle ^{2} \right) e^{1.0 i {\omega}l t} -1 \sqrt{\xi \kappa} \langle a\rangle \langle {\sigma}{1}^{{12}} {\sigma}{2}^{{12}}\rangle e^{1.0 i {\omega}l t} \right) \ \frac{d}{dt} \langle {\sigma}{1}^{{22}} {\sigma}{2}^{{22}}\rangle &= 2 \mathrm{pulse}\left( t \right) \langle {\sigma}{1}^{{22}}\rangle \eta -2.0 \langle {\sigma}{1}^{{22}} {\sigma}{2}^{{22}}\rangle \mathrm{pulse}\left( t \right) \eta + dW/dt \left( -1 \langle {\sigma}{1}^{{22}} {\sigma}{2}^{{22}}\rangle \sqrt{\xi \kappa} e^{-1.0 i {\omega}l t} \langle a^\dagger\rangle + \left( 2 \langle {\sigma}{1}^{{22}}\rangle \langle a^\dagger {\sigma}{1}^{{22}}\rangle -2 \langle {\sigma}{1}^{{22}}\rangle ^{2} \langle a^\dagger\rangle + \langle {\sigma}{1}^{{22}} {\sigma}{2}^{{22}}\rangle \langle a^\dagger\rangle \right) \sqrt{\xi \kappa} e^{-1.0 i {\omega}l t} -1 \langle {\sigma}{1}^{{22}} {\sigma}{2}^{{22}}\rangle \sqrt{\xi \kappa} \langle a\rangle e^{1.0 i {\omega}l t} + \left( \langle {\sigma}{1}^{{22}} {\sigma}{2}^{{22}}\rangle \langle a\rangle + 2 \langle {\sigma}{1}^{{22}}\rangle \langle a {\sigma}{1}^{{22}}\rangle -2 \langle {\sigma}{1}^{{22}}\rangle ^{2} \langle a\rangle \right) \sqrt{\xi \kappa} e^{1.0 i {\omega}l t} \right) -2.0 \langle {\sigma}{1}^{{22}} {\sigma}{2}^{{22}}\rangle \gamma + 1 i \left( 2 \langle a^\dagger {\sigma}{1}^{{12}}\rangle \langle {\sigma}{1}^{{22}}\rangle -4 \langle {\sigma}{1}^{{22}}\rangle \langle {\sigma}{1}^{{12}}\rangle \langle a^\dagger\rangle + 2 \langle {\sigma}{1}^{{12}} {\sigma}{2}^{{22}}\rangle \langle a^\dagger\rangle + 2 \langle a^\dagger {\sigma}{1}^{{22}}\rangle \langle {\sigma}{1}^{{12}}\rangle \right) g -1 i g \left( 2 \langle {\sigma}{1}^{{22}}\rangle \langle a {\sigma}{1}^{{21}}\rangle + 2 \langle {\sigma}{1}^{{21}}\rangle \langle a {\sigma}{1}^{{22}}\rangle -4 \langle {\sigma}{1}^{{21}}\rangle \langle {\sigma}{1}^{{22}}\rangle \langle a\rangle + 2 {\langle {\sigma}{1}^{{12}} {\sigma}_{2}^{{22}}\rangle ^{}} \langle a\rangle \right) 
# \end{align}
# ```

# Here we define the actual values for the system parameters. We then show that the deterministic time evolution without the noise terms is still accessible by using the constructor `System` for the stochastic system of equations and the syntax for the simulation of the time evolution is as usual.

## frequencies are in kHz
ωc_ = 0.0
κ_ = 2π*1.13e3
ξ_ = 0.12
N_ = 5e4
ωa_ = 0.0
γ_ = 2π*0.375
η_ = 2π*20
χ_ = 0.016
g_ = 2π*0.73
ωl_ = 2π*1e3
t0 = 0.0
t1 = 20e-3
p = [N, ωa, γ, η, χ, ωc, κ, g, ξ, ωl]
p0 = [N_, ωa_, γ_, η_, χ_, ωc_, κ_, g_, ξ_, ωl_]
T_end = 0.1 # 0.1ms

@named sys = System(scaled_eqs)
u0 = zeros(ComplexF64, length(scaled_eqs))
dict = merge(Dict(unknowns(sys) .=> u0), Dict(p .=> p0))
prob = ODEProblem(sys, dict, (0.0, T_end))
sol_det = solve(prob, RK4(), dt = T_end/2e5)

graph = plot(xlabel = "Time (ms)", ylabel = "Cavity photon number")
plot!(graph, sol_det.t, real(sol_det[a'a]), legend = false)

# The stochastic time evolution is accessible via the constructor `SDESystem`, whose syntax is exactly the same as for the System, but with keyword args as defined in the [SDE tutorial](https://docs.sciml.ai/DiffEqDocs/stable/tutorials/sde_example/). We need to provide a noise process for the measurement. If the noise is white the appropriate noise process is a Wiener process. The `SDEProblem` is constructed just as the `ODEProblem`, but with an additional noise argument.

# We can then make use of the `EnsembleProblem`, which automatically runs multiple instances of the stochastic equations of motion. The number of trajectories can be set in the solve call. See the tutorial linked above for more details of the function calls here.

@named sys_st = SDESystem(scaled_eqs)

Random.seed!(2) # hide
noise = StochasticDiffEq.RealWienerProcess(0.0, 0.0)
prob_st = SDEProblem(sys_st, dict, (0.0, T_end); noise = noise)
sol_test = solve(prob_st, EM(); dt = T_end/2e5);

plot(
    sol_test.t,
    real(sol_test[a'a]),
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
    dt = T_end/2e5,
    save_noise = true,
    trajectories = traj,
    saveat = tspan,
)

# We plot the average of the cavity photon number for the stochastic and deterministic equation of motion with the trajectories in gray in the background. We can see that the dynamics of the system is indeed modified by the measurement backaction.

n_avg_ = zeros(length(tspan)) # average photon number
a_real_avg_ = zeros(length(tspan)) # average field (real part)
traj_succ = zeros(traj) # trajectories can be numerically unstable 
for i = 1:traj
    sol_ = sol.u[i]
    if length(sol_) == length(tspan)
        n_avg_ .+= real(sol_[a'a])
        a_real_avg_ .+= real(sol_[a])
        traj_succ[i] = 1
    end
end

@show sum(traj_succ)
n_avg = n_avg_ / sum(traj_succ)
a_real_avg = a_real_avg_ / sum(traj_succ)
graph1 = plot(xlabel = "Time (ms)", ylabel = "Cavity photon number")
for i = 1:traj
    sol_n = real(sol.u[i][a'a])
    plot!(
        graph1,
        sol.u[i].t,
        real(sol.u[i][a'a]),
        color = :grey,
        alpha = 0.75,
        label = nothing,
    )
end
plot!(graph1, tspan, n_avg, color = :red, label = nothing)

# For the following cavity field amplitude we see that the ensemble average becomes zero, but the trajectories have finite cavity field amplitude

graph2 = plot(xlabel = "Time (ms)", ylabel = "Real part cavity field")
for i = 1:traj
    plot!(graph2, tspan, real(sol.u[i][a']), color = :grey, alpha = 0.75, label = nothing)
end
plot!(graph2, tspan, a_real_avg, color = :red, label = nothing)

# ## Package versions

# These results were obtained using the following versions:

using InteractiveUtils
versioninfo()

using Pkg
Pkg.status(["QuantumCumulants", "OrdinaryDiffEq"], mode = PKGMODE_MANIFEST)
