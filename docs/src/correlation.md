# Correlation functions and spectra

Two quantities that are often of interest in a system are two-time correlation functions and spectral densities. Consider, for example, two operators ``a`` and ``b``. Their two-time correlation function is given by
```math
g(t,\tau) = \langle a(t+\tau) b(t)\rangle.
```
As we can see, we need to take some care here since ``g`` depends on two different times. The spectral density corresponding to ``a`` and ``b`` is given by the Fourier transform,
```math
S(\omega) = 2\text{Re}\left\{\int d\tau e^{-i\omega\tau}g(t,\tau)\right\}.
```
So in any case, we need to treat the two-time correlation function before we can obtain the spectrum.

## Correlation function

In order to compute a correlation function, we first evolve a system of equations up to a time ``t``. Then, we can derive another set of equations that describe the correlation function. This new set of equations is then evolved from time ``t`` up to a time ``t+\tau``. The correlation function is then stored in the first entry of the result. The initial state of the set of equations describing the correlation function will be determined by the state of the original system at time ``t``.

Note, that whenever an instance of a [`CorrelationFunction`](@ref) is created, a set of equations is derived using a custom version of the [`complete`](@ref) function. Depending on the size and order of the considered system, this can take some time. An important distinction that can eventually reduce the computation time by quite a bit is whether or not the original system has been evolved up to steady state, i.e. if the system does not change its state after ``t``. This can be controlled with the keyword argument `steady_state=true` when construction the [`CorrelationFunction`](@ref).

To clarify the usage, consider the simple case of a cavity with resonance frequency $\omega_\mathrm{c}$ that initially has a finite number of photons inside which decay over time at a rate ``\kappa``. We want to compute the two-time correlation function of the field (first-order degree of coherence) given by
```math
g(t,\tau) = \langle a^\dagger(t+\tau)a(t)\rangle.
```
Given the Hamiltonian ``H = \omega_\mathrm{c}a^\dagger a`` and the collapse operator ``a``, it is easy to derive the equation
```math
\frac{d}{d\tau}a^\dagger(t+\tau) = (i\omega_\mathrm{c} - \frac{\kappa}{2})a^\dagger(t+\tau).
```
Since ``a(t)`` is independent of ``\tau`` we can simply multiply the above equation with ``a(t)`` from the right and average to obtain
```math
\frac{d}{d\tau}\langle a^\dagger(t+\tau) a(t)\rangle = (i\omega_\mathrm{c} - \frac{\kappa}{2})\langle a^\dagger(t+\tau) a(t)\rangle.
```
Note that this is generally valid, and will lead to a system of equations that is linear in ``a(t)``. The above is the equation of motion for the correlation function ``g(t,\tau)``. In this very simple case we can solve it analytically and find
```math
g(t,\tau) = \langle a^\dagger(t)a(t)\rangle e^{(i\omega_\mathrm{c} - \kappa/2)\tau}.
```

This is the essential procedure with which correlation functions can be computed within **Qumulants.jl**. In code, the above is just:
```@example correlation
using Qumulants # hide
h = FockSpace(:cavity)
a = Destroy(h,:a)
@parameters ωc κ
H = ωc*a'*a
he = average(heisenberg(a'*a,H,[a];rates=[κ]))

c = CorrelationFunction(a', a, he)
nothing # hide
```
When the [`CorrelationFunction`](@ref) is constructed, an additional Hilbert space is added internally which represents the system at the time ``t``. In our case, this means that another [`FockSpace`](@ref) is added. Note that all operators involved in the correlation function are defined on the [`ProductSpace`](@ref) including this additional Hilbert space.

The equation for ``g(t,\tau)`` is now stored in the first entry of `c.de`. To solve the above numerically, we need to generate code and solve the equations numerically.
```@example correlation
using OrdinaryDiffEq

f = generate_ode(he,(ωc,κ))
n0 = 20.0 # Initial number of photons in the cavity
u0 = [n0]
p0 = (1,1)
prob = ODEProblem(f,u0,(0.0,2.0),p0) # End time not in steady state
sol = solve(prob,RK4())
nothing # hide
```
Numerical computing the correlation function works in the same way. Note, the initial state of the correlation function depends on the final state of the system. However, in general it does not depend on *all* the final values of the system. The correct values can be picked out automatically using the `initial_values` function.
```@example correlation
cf = generate_ode(c, (ωc,κ))
u0_c = initial_values(c, sol.u[end])
prob_c = ODEProblem(cf,u0_c,(0.0,10.0),p0)
sol_c = solve(prob_c,RK4(),save_idxs=1)
nothing # hide
```
Finally, lets check our numerical solution against the analytic one obtained above:
```@example correlation
using Test # hide
g_analytic(τ) = @. sol.u[end] * exp((im*p0[1]-0.5p0[2])*τ)
@test isapprox(sol_c.u, g_analytic(sol_c.t), rtol=1e-4)
```

Note, that this was a very simple case. Usually the system of equations describing the correlation function is much more complex and depends on multiple other correlation functions (see for example [Spectrum of a single atom laser](@ref)).


## Spectrum

There are two possible ways two compute the spectrum given a correlation function:

1. Solving the differential equation needed to obtain ``g(t,\tau)`` and taking the Fourier transform.
2. Taking the (symbolic) Laplace transform of the system of equations describing a correlation function.

On the one hand, the first approach works generally, but is computationally more intense. The second approach, on the other hand, yields a simple linear system of equations which is easy to solve, but only works when the correlation function has been computed starting from the steady state. Both methods can be easily used with **Qumulants.jl**.


### Numerical solution of ``g(t,\tau)``

As mentioned above, this approach works generally, regardless of whether the system is in steady state at time ``t``. However, it has some computational drawbacks. Additionally, the spectrum is not always well defined when not in steady state. This approach is the same as the one used in [**QuantumOptics.jl**](https://qojulia.org), and we can borrow the implemented FFT function from there:
```@example correlation
τ = collect(range(0.0, sol_c.t[end], length=101)) # need equidistant list of times for FFT
using QuantumOptics.timecorrelations: correlation2spectrum
ω, s = correlation2spectrum(τ, sol_c.(τ))
nothing # hide
```
The spectrum obtained in this way roughly has a FWHM of `κ` and is based around the chosen `ωc`. The fact that the FWHM is not *exactly* `κ` illustrates the computation drawback: in order to obtain the correct FWHM we would have to increase the integration time by orders of magnitude. For larger system, this can be computationally expensive.


### Steady state: using the Laplace transform

A useful property of the two-time correlation function is that, if the system is in steady state at time ``t``, then the system of equations is linear, since ``b(t)`` can occur at most once in each product. We can therefore write any system of equations describing the correlation function as
```math
\frac{d}{d\tau} \textbf{y}(\tau) = \textbf{M} \textbf{y}(\tau) + \textbf{c},
```
where ``\textbf{y}(\tau)`` is the vector containing the left-hand-side of the correlation function system (``g(t,\tau) \equiv y_1(\tau)``). The matrix ``\textbf{M}`` contains coefficients consisting of parameters and steady-state values, and is independent of time, and the vector ``\textbf{c}`` includes other constants.

We define ``\textbf{x}(s) = \mathcal{L}\left(\textbf{y}(\tau)\right)``, i.e. ``\textbf{x}(s)`` is the Laplace transform of ``\textbf{y}(\tau)`` with respect to ``\tau``. Applying the Laplace transform to the differential equation above, we obtain
```math
(s - \textbf{M})\textbf{x}(s) = \textbf{y}(0) + \frac{\textbf{c}}{s}.
```
The Laplace transform is equivalent to the Fourier transform at the point ``s=i \omega``, i.e. the spectrum is given by ``S(\omega) = 2\text{Re}\left\{x_1(i\omega)\right\}``. Therefore, we can reduce the task to solving the equation
```math
A\textbf{x} = b,
```
where ``A = i\omega - \textbf{M}`` and ``b = \textbf{y}(0) + \textbf{c}/(i\omega)``. In most cases, solving the above matrix equation is much faster than doing an additional time evolution to obtain the correlation function.

This approach is implemented with the [`Spectrum`](@ref) type, which performs the Laplace transform and generates a function that returns the matrix ``A`` and the vector ``b`` in numerical form depending on the steady-state values and given parameters. This can be used as follows

```@example correlation
c = CorrelationFunction(a', a, he; steady_state=true) # need to specify steady state
S = Spectrum(c,(ωc,κ))
nothing # hide
```

The above performs the Laplace transform on a symbolic level (i.e. it derives the matrix ``A``). To actually compute the spectrum, we can do

```@example correlation
s = S(ω,sol.u[end],p0)
nothing # hide
```

## Examples:

* [Mollow Triplet](@ref)
* [Spectrum of a single atom laser](@ref)
