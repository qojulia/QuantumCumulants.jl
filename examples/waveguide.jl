# # Waveguide Energy Transfer (meanfield setup)

# Collective behaviour of atomic ensembles coupled to a waveguide. The full
# example builds an initial quantum state via `coherentspinstate` and uses
# `initial_values(eqs, ψ)`; that path is being re-implemented in v1.0 and is
# omitted here. The Hamiltonian and equation-set construction work as below.

using QuantumCumulants
using OrdinaryDiffEq, ModelingToolkitBase
using Plots

M = 2

@variables Ω12 Γ12 Γ11 Γ22

h_spin(i) = SpinSpace(Symbol(:spin_, i))
h = tensor([h_spin(i) for i = 1:M]...)

Sx(i) = Spin(h, Symbol(:S, i), 1, i)
Sy(i) = Spin(h, Symbol(:S, i), 2, i)
Sz(i) = Spin(h, Symbol(:S, i), 3, i)

Sm(i) = Sx(i) - 1im * Sy(i)
Sp(i) = Sx(i) + 1im * Sy(i)

H = Ω12 * (Sp(1) * Sm(2) + Sp(2) * Sm(1))

J = [Sm(1), Sm(2)]
rates = [Γ11, Γ22]

ops = [Sz(1), Sz(2)]
eqs = meanfield(ops, H, J; rates = rates, order = 2, simplify = false)
println("Generated ", length(eqs.equations), " base meanfield equations.")
println("First equation: ", eqs.equations[1])
