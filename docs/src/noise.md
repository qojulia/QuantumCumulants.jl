# Noise & measurement backaction

The [mean-field pipeline](@ref "The mean-field pipeline") derives the *unconditional* dynamics: the Lindblad master equation averaged over all measurement outcomes, where the white-noise term drops out (see the [theory section](@ref theory)). When a decay channel is *continuously monitored*, that is no longer the whole story. The measurement record feeds back on the state (*measurement backaction*), and the conditional evolution is a stochastic master equation driven by the measurement noise.

QuantumCumulants derives the cumulant-expanded equations for this conditional evolution directly. The result is a [`NoiseMeanfieldEquations`](@ref): each equation carries a deterministic drift *and* a noise drift proportional to a Wiener increment ``dW``.

## Requesting backaction: the `efficiencies` keyword

Noise is opt-in. Pass `efficiencies` to [`meanfield`](@ref), one entry per collapse operator, giving the detection efficiency of that channel. An efficiency of `0` means the channel is dissipative but *not* monitored (it contributes to the drift only); a non-zero efficiency adds the corresponding measurement backaction:

```julia
using QuantumCumulants

# ... build H and the jump operators J ...
J            = [a, σ(1, 2, j), σ(2, 1, j), σ(2, 2, j)]
rates        = [κ, γ, η, χ]
efficiencies = [ξ, 0, 0, 0]   # only the cavity-output channel is measured

eqs = meanfield(ops, H, J; rates=rates, efficiencies=efficiencies, order=2)
```

Because `efficiencies` is supplied, `eqs` is a [`NoiseMeanfieldEquations`](@ref) rather than a [`MeanfieldEquations`](@ref). Each averaged equation now has the schematic form

```math
d\langle \mathcal{O}\rangle = \underbrace{(\dots)}_{\text{drift}}\,dt + \underbrace{(\dots)}_{\text{noise drift}}\,dW,
```

where the noise drift is the measurement-conditioned term. The choice of jump operator sets the measurement: for example a cavity output ``\hat a\, e^{i\omega_l t}`` realises heterodyne detection at local-oscillator frequency ``\omega_l``.

## Closing, scaling, evaluating

[`complete`](@ref), [`scale`](@ref), and [`evaluate`](@ref) operate on a `NoiseMeanfieldEquations` exactly as they do on a deterministic one; the noise column is carried through unchanged:

```julia
eqs_c      = complete(eqs)
scaled_eqs = scale(eqs_c)   # permutation-symmetric ensemble
```

## Solving: deterministic vs. stochastic

Converting to a `ModelingToolkitBase.System` keeps both columns. How you *solve* it selects which you use:

- An `ODEProblem` integrates the **drift only**, recovering the unconditional (ensemble-mean) trajectory.

- An `SDEProblem` integrates **drift + noise**, producing a single conditional trajectory for one measurement realisation. It needs a noise process; for white measurement noise this is a Wiener process from [StochasticDiffEq.jl](https://github.com/SciML/StochasticDiffEq.jl).

```julia
using ModelingToolkitBase, OrdinaryDiffEq, StochasticDiffEq

sys = System(scaled_eqs; name=:sys)
dict = merge(Dict(unknowns(sys) .=> u0), Dict(p .=> p0))

# Unconditional (drift-only) evolution:
sol_det = solve(ODEProblem(sys, dict, (0.0, T)), Tsit5())

# A single conditional trajectory:
noise   = StochasticDiffEq.RealWienerProcess(0.0, 0.0)
sol_traj = solve(SDEProblem(sys, dict, (0.0, T); noise=noise), EM(); dt=T/2e5)
```

Averaging many conditional trajectories (via an `EnsembleProblem`) reproduces the deterministic result, while individual trajectories show the backaction. [`get_solution`](@ref) extracts averages from either kind of solution.

A complete, runnable walkthrough, including the ensemble averaging and plots, is the [heterodyne detection example](@ref "Heterodyne detection of emission from atomic ensemble").

## Retrodiction

Continuous monitoring also supports *retrodiction*: using the measurement record to estimate the state at an *earlier* time (the smoothing / past-quantum-state problem). This is the backward evolution, requested with `direction=Backward()`:

```julia
eqs_backward = meanfield(ops, H, J; rates=rates, efficiencies=efficiencies,
                         direction=Backward(), order=2)
```

`Backward()` flips the sign of the commutator and uses the adjoint Lindblad recycling with the trace-preserving term, producing the backward (effect-operator) equations. The [`EvolutionDirection`](@ref) tag is stored on the equations and threaded through `System`, so the rest of the workflow (complete, scale, solve) is unchanged. Combining a forward filter with a backward filter yields the smoothed estimate.
