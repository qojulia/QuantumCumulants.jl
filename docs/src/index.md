```@raw html
---
layout: home

hero:
  name: QuantumCumulants.jl
  text: Generalized mean-field equations for open quantum systems
  tagline: Derive Heisenberg equations for expectation values symbolically, truncate them with the cumulant expansion at any order, and solve the closed system numerically.
  actions:
    - theme: brand
      text: Get started
      link: manual/getting_started/
    - theme: alt
      text: Examples
      link: examples/single-atom-laser-spectrum/
    - theme: alt
      text: Operator algebra (SQA)
      link: https://qojulia.github.io/SecondQuantizedAlgebra.jl/stable/
    - theme: alt
      text: View on GitHub
      link: https://github.com/qojulia/QuantumCumulants.jl
  image:
    src: assets/logo.png
    alt: QuantumCumulants.jl logo

features:
  - icon: "⟨·⟩"
    title: Cumulant expansion
    details: Truncate the hierarchy of moment equations at any order, with mixed orders per subsystem.
  - icon: "∞"
    title: Automatic closure
    details: Find every moment the truncated equations depend on and add its equation until the system closes.
  - icon: "Σ"
    title: Many-body systems
    details: Keep the number of emitters symbolic with indexed operators and sums, or scale identical subsystems to a single representative.
  - icon: "ω"
    title: Correlation functions and spectra
    details: Derive two-time correlation functions and compute emission spectra via the quantum regression theorem.
  - icon: "dW"
    title: Measurement backaction
    details: Derive stochastic mean-field equations for continuously monitored systems.
  - icon:
      src: https://sciml.ai/assets/favicon.png
      alt: SciML logo
      wrap: true
    title: ModelingToolkit bridge
    details: Turn the equations into a ModelingToolkit System and solve them with the SciML differential equation solvers.
---
```

```@meta
CurrentModule = QuantumCumulants
```

QuantumCumulants derives equations of motion for operator averages from the fundamental commutation relations and closes them with the generalized cumulant expansion [(Plankensteiner et al., 2022)](https://doi.org/10.22331/q-2022-01-04-617). The operator algebra comes from [SecondQuantizedAlgebra.jl](@extref SecondQuantizedAlgebra index), which QuantumCumulants re-exports. Its documentation covers how to build a model (Hilbert spaces, operators, indices, sums), and this one covers what QuantumCumulants adds on top of it.

## Quick start

Install the package with Julia's package manager:

```julia-repl
pkg> add QuantumCumulants
```

A first-order mean-field model of a laser with a single atom as gain medium:

```@example home
using QuantumCumulants

h = FockSpace(:cavity) ⊗ NLevelSpace(:atom, (:g, :e))

@variables Δ::Real g::Real κ::Real γ::Real ν::Real
@qnumbers a::Destroy(h)
σ(i, j) = Transition(h, :σ, i, j)

H = Δ * a' * a + g * (a' * σ(:g, :e) + a * σ(:e, :g))
J = [a, σ(:g, :e), σ(:e, :g)]
rates = [κ, γ, ν]

eqs = meanfield([a, σ(:g, :e), σ(:e, :e)], H, J; rates = rates, order = 1)
```

The equations become a ModelingToolkit `System` and are integrated with OrdinaryDiffEq:

```@example home
using ModelingToolkitBase, OrdinaryDiffEqLowOrderRK, Plots

sys = mtkcompile(System(eqs; name = :laser))
p0 = [Δ => 0, g => 1.5, κ => 1, γ => 0.25, ν => 4]
u0 = ComplexF64[1.0e-2, 0, 0]
prob = ODEProblem(sys, merge(initial_values(eqs, u0), Dict(p0)), (0.0, 50.0))
sol = solve(prob, RK4())

ts = range(0.0, 50.0; length = 200)
plot(ts, abs2.(get_solution(sol, a, eqs).(ts)); xlabel = "t", label = "|⟨a⟩|²")
```

The [manual](@ref "Getting started") walks through each step in detail.

## Citing

If you use QuantumCumulants.jl in your research, please cite [Quantum 6, 617 (2022)](https://doi.org/10.22331/q-2022-01-04-617):

```bib
@article{plankensteiner2022quantumcumulants,
  doi = {10.22331/q-2022-01-04-617},
  title = {Quantum{C}umulants.jl: {A} {J}ulia framework for generalized mean-field equations in open quantum systems},
  author = {Plankensteiner, David and Hotter, Christoph and Ritsch, Helmut},
  journal = {{Quantum}},
  volume = {6},
  pages = {617},
  year = {2022}
}
```
