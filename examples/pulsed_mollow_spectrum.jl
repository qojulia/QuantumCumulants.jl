# # Mollow Triplet from a Pulsed Drive

# The [Mollow Triplet](@ref) example computes the resonance fluorescence spectrum of a *constantly* driven atom. Here we drive the atom with a **time-dependent** field that is smoothly switched on, and compute the emission spectrum from the two-time correlation function once the atom has settled into a quasi-steady state. This showcases correlation functions of time-dependent Hamiltonians.

# The atom obeys

# $H(t) = \Delta\,\sigma^{ee} + r(t)\,\Omega\,(\sigma^{ge} + \sigma^{eg}),$

# where the dimensionless envelope $r(t)$ ramps the Rabi drive from $0$ to $1$. The emission spectrum is the Fourier transform of the first-order correlation

# $g(t_0,\tau) = \langle \sigma^{eg}(t_0+\tau)\,\sigma^{ge}(t_0)\rangle.$

# By the quantum regression theorem, with a time-dependent generator the $\tau$-evolution is governed by $H(t_0+\tau)$, **not** $H(\tau)$: the drive keeps running during the correlation delay, starting from the absolute time $t_0$ at which the original evolution stopped. **QuantumCumulants.jl** handles this by substituting the original time variable $t \to t_0 + \tau$ in the correlation equations. We pass $t_0$ via the `iv0` keyword and supply its value alongside the other parameters.

using QuantumCumulants
using ModelingToolkitBase, OrdinaryDiffEqTsit5, OrdinaryDiffEqLowOrderRK
using QuantumOptics: timecorrelations
using Plots

# We register the envelope $r(t)$ as a time-dependent function. As in the [Ramsey Spectroscopy](@ref) example, we seed a `meanfield` call first so the envelope is built on the system's own independent variable `t`.

h = NLevelSpace(:atom, (:g, :e))
σ(i, j) = Transition(h, :σ, i, j)

@variables Δ Ω γ
@variables t₀::Real   # absolute time at which the drive evolution stops (t₀)
@register_symbolic r(t)

eqs_seed = meanfield([σ(:e, :e), σ(:e, :g)], Δ * σ(:e, :e), [σ(:g, :e)]; rates = [γ])
t = eqs_seed.iv

H = Δ * σ(:e, :e) + r(t) * Ω * (σ(:g, :e) + σ(:e, :g))
J = [σ(:g, :e)]

eqs = meanfield([σ(:e, :e), σ(:e, :g)], H, J; rates = [γ], iv = t)
complete!(eqs)

# We drive on resonance ($\Delta = 0$) with $\Omega = 3\gamma$, switching the field on smoothly with $r(t) = 1 - e^{-t/t_r}$ so it has effectively saturated long before the stop time $t_0$.

r(t) = 1 - exp(-t / 2.0)

γv, Ωv, Δv = 1.0, 3.0, 0.0
ps = [γ, Ω, Δ]
p0 = [γv, Ωv, Δv]

sys = mtkcompile(System(eqs; name = :atom))
u0 = initial_values(eqs, zeros(ComplexF64, length(eqs)))
tstop = 20.0
prob = ODEProblem(sys, merge(u0, Dict(ps .=> p0)), (0.0, tstop))
sol = solve(prob, Tsit5(); abstol = 1.0e-10, reltol = 1.0e-10)

plot(
    sol.t,
    real.(get_solution(sol, σ(:e, :e), eqs).(sol.t)),
    xlabel = "γt",
    ylabel = "⟨σᵉᵉ⟩",
    label = "excited-state population",
    size = (600, 300),
)

# The population settles onto a plateau: at `t₀ = tstop` the atom is in a quasi-steady state under the (now constant) drive. We build the emission correlation function on this system, passing `iv0 = t₀`.

c = CorrelationFunction(σ(:e, :g), σ(:g, :e), eqs; iv0 = t₀)
nothing # hide

# Without `iv0`, a time-dependent system raises an informative error, since the correlation equations would otherwise contain the orphaned time variable `t`:

try
    CorrelationFunction(σ(:e, :g), σ(:g, :e), eqs)
catch e
    println(e isa ArgumentError ? "ArgumentError: iv0 is required" : e)
end

# We solve the $\tau$-evolution. The steady-state initial values are read off the original solution with [`correlation_u0`](@ref); the parameters, **including** `t₀`, are propagated with [`correlation_p0`](@ref).

csys = mtkcompile(System(c; name = :corr))
u0_c = correlation_u0(c, sol.u[end])
p0_c = correlation_p0(c, sol.u[end], [γ => γv, Ω => Ωv, Δ => Δv, t₀ => tstop])

τ_end = 30.0
prob_c = ODEProblem(csys, merge(u0_c, Dict(p0_c)), (0.0, τ_end))
sol_c = solve(prob_c, Tsit5(); abstol = 1.0e-10, reltol = 1.0e-10, save_idxs = 1)
nothing # hide

# The spectrum is the Fourier transform of the *inelastic* (connected) part of the correlation, $g(t_0,\tau) - |\langle\sigma^{ge}\rangle|^2$, with the coherent (elastically scattered) contribution removed. We borrow the FFT helper from [**QuantumOptics.jl**](https://qojulia.org).

coh = abs2(get_solution(sol, σ(:g, :e), eqs)(tstop))   # |⟨σᵍᵉ⟩|² at t₀ (elastic peak)
τ = collect(range(0.0, τ_end; length = 2001))          # equidistant grid for the FFT
g_inel = sol_c.(τ) .- coh
ω, S = timecorrelations.correlation2spectrum(τ, g_inel)
nothing # hide

# To confirm the result, we compare against the textbook **constant-drive** Mollow triplet: the same atom driven by a time-independent field $\Omega$, evolved to steady state, with its correlation computed the standard way (no `iv0`).

H_const = Δ * σ(:e, :e) + Ω * (σ(:g, :e) + σ(:e, :g))
eqs_const = meanfield([σ(:e, :e), σ(:e, :g)], H_const, J; rates = [γ])
complete!(eqs_const)

sys_const = mtkcompile(System(eqs_const; name = :atom_const))
sol_const = solve(
    ODEProblem(sys_const, merge(initial_values(eqs_const, zeros(ComplexF64, 2)), Dict(ps .=> p0)), (0.0, 40.0)),
    Tsit5(); abstol = 1.0e-10, reltol = 1.0e-10,
)

c_const = CorrelationFunction(σ(:e, :g), σ(:g, :e), eqs_const)
csys_const = mtkcompile(System(c_const; name = :corr_const))
sol_cc = solve(
    ODEProblem(
        csys_const,
        merge(
            correlation_u0(c_const, sol_const.u[end]),
            Dict(correlation_p0(c_const, sol_const.u[end], [γ => γv, Ω => Ωv, Δ => Δv]))
        ),
        (0.0, τ_end),
    ),
    Tsit5(); abstol = 1.0e-10, reltol = 1.0e-10, save_idxs = 1,
)
coh_const = abs2(get_solution(sol_const, σ(:g, :e), eqs_const)(40.0))
_, S_const = timecorrelations.correlation2spectrum(τ, sol_cc.(τ) .- coh_const)
nothing # hide

# The two spectra lie on top of each other: the pulsed drive, sampled from a plateau time $t_0$ via `iv0`, reproduces the steady-state Mollow triplet. The sidebands sit at $\omega \approx \pm 2\Omega$ (the dressed-state splitting at resonance), flanking the central line at $\omega = 0$.

plot(ω, S; xlims = (-12, 12), xlabel = "ω", ylabel = "S(ω)", label = "pulsed drive (iv0 = t₀)", lw = 2, size = (600, 350))
plot!(ω, S_const; xlims = (-12, 12), label = "constant drive (steady state)", ls = :dash, lw = 2)

# ## Package versions

# These results were obtained using the following versions:

using InteractiveUtils
versioninfo()

using Pkg
Pkg.status(
    ["QuantumCumulants", "OrdinaryDiffEqTsit5", "ModelingToolkitBase", "QuantumOptics", "Plots"],
    mode = PKGMODE_MANIFEST,
)
