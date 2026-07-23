# [Convergence and validity of the cumulant expansion](@id convergence)

The [cumulant expansion](@ref theory) closes an open hierarchy of operator
averages by neglecting cumulants above a chosen order ``M``. This raises two
practical questions: does increasing ``M`` improve the result, and when can a
truncation be trusted?

The complete hierarchy is an exact reformulation of the dynamics, but its
truncation is uncontrolled in general. A higher order retains more correlations
without necessarily producing a more accurate solution. Additional structure,
such as Gaussianity or a systematic ``1/N`` expansion, is needed to make the
neglected cumulants parametrically small.

This does not make cumulant closure a bad approximation. It is often extremely
useful, and it is rigorous in identifiable limits. The important distinction is
between accuracy demonstrated for a particular problem and a general guarantee,
which the method does not provide.

## Practical guidance for QuantumCumulants.jl

!!! tip "How to use the truncation order responsibly"
    * Compare consecutive affordable orders for the observables and time interval
      of interest. Agreement is useful evidence, although it is not an error
      bound.
    * When possible, validate a small instance against a full master-equation
      solution, for example with [QuantumOptics.jl](https://qojulia.org/).
    * Check physical bounds and numerical stability. Negative populations,
      expectation values outside their allowed range, and runaway oscillations
      are warning signs.
    * Check the dependence on system size as well as on truncation order. Large
      ``N`` alone does not guarantee that a closure becomes accurate.

## What truncation assumes

The transformation between moments and cumulants is triangular and invertible:
the complete cumulant hierarchy contains the same information as the complete
moment hierarchy. For a finite-dimensional operator algebra, a complete set of
moments determines the state, so retaining the full hierarchy reproduces the
exact dynamics. For bosonic systems with an infinite-dimensional Hilbert space,
the corresponding statement also requires the usual assumptions that make the
[moment problem](https://en.wikipedia.org/wiki/Moment_problem) determinate.

For ``N`` finite-level subsystems, local operator products can be reduced to a
finite operator basis. In an atom-only model of ``N`` two-level systems, for
example, order ``M=N`` can retain every inter-atomic correlation and reproduce
the full dynamics. This does not apply to an untruncated bosonic mode, whose
operator hierarchy remains infinite.

The practical issue is therefore the behavior of affordable intermediate
orders. The order ``M`` is not the power of a coupling constant: unlike a
perturbative expansion, an order-``M`` calculation contains no explicit factor
that forces order-``M+1`` cumulants to be smaller. In a strongly non-Gaussian
state, high-order cumulants need not decay rapidly with order.

QuantumCumulants.jl closes the hierarchy by setting joint cumulants above the
selected order to zero and rewriting the corresponding higher moments through
the moment--cumulant identity. For noncommuting operators, the operator order
within each moment is retained. The precise joint-cumulant definition is given
in [Theoretical background](@ref theory);
[Bianucci and Bologna (2020)](https://arxiv.org/abs/1911.09620) discuss the care
required when extending classical cumulant properties to noncommuting operators.

## Why convergence is not guaranteed

There is no general convergence theorem or a priori applicability criterion for
the dynamical cumulant closure. Its accuracy can depend on the Hamiltonian,
dissipation, initial state, observable, time interval, and the scaling of model
parameters with system size. Agreement for one observable or at one time does
not establish accuracy everywhere.

Successive closures may behave in several ways. They may converge toward a
single result, form an [asymptotic
sequence](https://en.wikipedia.org/wiki/Asymptotic_expansion) that first improves
and then deteriorates, or fail to approach a common limit. Even and odd orders
can form different subsequences. Consequently, order-by-order agreement is an
empirical diagnostic rather than a proof of convergence.

[Kerber, Ritsch and Ostermann (2025)](https://arxiv.org/abs/2511.20115) emphasize
that a general criterion for the applicability and convergence of higher-order
cumulant expansions remains to be found. Rigorous results for special settings
should therefore not be promoted to a guarantee for general driven-dissipative
dynamics.

## Realizability and the special role of second order

Convergence is not the only concern: a closed set of approximate moments may
fail to correspond to any physical state. There is a precise reason why second
order occupies a special position.

For a classical random variable ``X`` with [characteristic
function](https://en.wikipedia.org/wiki/Characteristic_function_%28probability_theory%29)
``\varphi(t)=\langle e^{itX}\rangle``, the cumulant generating function is

```math
K(t)=\log\varphi(t)=\sum_{n\geq 1}\kappa_n\frac{(it)^n}{n!}.
```

[Marcinkiewicz's theorem](https://link.springer.com/article/10.1007/BF01007988)
([original paper in French,
1939](https://eudml.org/doc/168829)) states that if
``K(t)`` is a polynomial, its degree is at most two. A non-Gaussian probability
distribution therefore cannot have an exactly terminating cumulant sequence.
[Rajagopal and Sudarshan
(1974)](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.10.1852) extended
this obstruction to many variables and to commuting and anticommuting
operator-valued variables. In the settings covered by these results, an exactly
terminating hierarchy describes either a [Gaussian
state](https://en.wikipedia.org/wiki/Gaussian_state), with cumulants through
second order, or a special first-order case such as a [coherent
state](https://en.wikipedia.org/wiki/Coherent_state); a genuinely non-Gaussian
state generally has cumulants at arbitrarily high orders.

A second-order closure consequently has a consistent interpretation as a
Gaussian ansatz. A quadratic Hamiltonian with linear dissipation gives closed
equations for first and second moments. Those moments evolve exactly for any
initial state; the entire state is described exactly when the initial state is
Gaussian and the dynamics preserves Gaussianity.

By contrast, setting all cumulants above a finite order ``M>2`` to zero does not,
in general, define the exact cumulants of a physical non-Gaussian state. This is
a structural warning, not a claim that every higher-order numerical solution
must become unphysical. Conversely, second-order closure alone does not
guarantee positivity. A realizability failure may appear as negative
populations, correlations outside their allowed bounds, or unstable and runaway
solutions.

## Regimes with additional control

The following regimes provide either exact closure or an independent reason for
the neglected cumulants to remain small.

| Regime | Source of control | Limitation |
|---|---|---|
| Gaussian states under Gaussian-preserving dynamics | Cumulants above second order vanish exactly | Does not cover dynamics that generates non-Gaussian correlations |
| Near-Gaussian dynamics | Higher cumulants remain small over the relevant time interval | The smallness is generally heuristic and must be checked for the model |
| Large ``N`` or high connectivity | A systematic ``1/N`` expansion can suppress connected correlations | Requires appropriate parameter scaling; large ``N`` alone is insufficient |

Outside such regimes, increasing ``M`` adds correlations but supplies no
independent control of the error.

## What benchmarks show

Benchmarks demonstrate several qualitatively different behaviors. They are
examples, not a classification theorem for Hamiltonians.

### Smooth improvement

For collective radiative decay in a dipole-dipole interacting atom chain, the
orders tested by Kerber, Ritsch and Ostermann (2025) improved smoothly. The gains
beyond second order were modest relative to their computational cost.
Order-by-order comparisons in ordered atomic arrays reached a compatible
conclusion ([Rubies-Bigorda, Ostermann and Yelin,
2023](https://doi.org/10.1103/PhysRevResearch.5.013091)).

### Higher-order artifacts and instability

The same study found very different behavior for a bi-prime factorization
problem formulated as adiabatic annealing. Going beyond mean field did not
improve the useful output. Although the full result was recovered at maximal
order, intermediate closures could agree less well than a lower-order one and,
even for a small system, developed fast oscillations, divergences, and partly
non-physical solutions. The model contains higher-body interaction terms, whose
correlations are not necessarily represented consistently at an intermediate
order.

Kerber, Ritsch and Ostermann summarize these contrasting examples as the
“good,” the “bad,” and the “ugly”: smooth improvement, higher-order artifacts,
and erratic dynamics, respectively.

### Non-uniform convergence

Central-spin models provide a particularly sharp caution
([Fowler-Wright, Arnardóttir, Kirton, Lovett and Keeling,
2023](https://journals.aps.org/prresearch/abstract/10.1103/PhysRevResearch.5.033148)):

1. The common expectation that mean field captures the exact ``N\to\infty``
   limit is false in general. Whether it does so depends on how the model
   parameters scale with ``N``.
2. Convergence can be non-uniform across even and odd orders. The two
   subsequences can approach distinct limits rather than improving monotonically
   toward one result.
3. The error need not decrease monotonically with ``N``. Even when a
   higher-order expansion recovers the correct large-``N`` limit, it can be less
   accurate than mean field at intermediate sizes.

A benchmark of Dicke superradiance gives a compatible but
observable-dependent picture ([Fasser, Genes, Ritsch and Holzinger,
2025](https://doi.org/10.3390/photonics12100996)). The time and magnitude of peak
emission converged reliably with order, whereas long-time populations were
captured only for small emitter numbers. Odd orders showed unphysical behavior,
and at large emitter number neither the individual-spin nor the collective-spin
expansion converged to the exact result. In that study, the individual-spin
expansion was more reliable than the collective-spin variant.

Together, these benchmarks show that reliability must be assessed for the
specific model, observable, and parameter regime.

## Related results and alternative closures

Several related methods use cumulants, moments, or cluster expansions, but they
do not all approximate the same object. Their results should therefore not be
read as convergence guarantees for the dynamical closure implemented by
QuantumCumulants.jl.

### Lessons from classical moment closure

Classical kinetic theories also produce open hierarchies, such as the [BBGKY
hierarchy](https://en.wikipedia.org/wiki/BBGKY_hierarchy). Closing such a
hierarchy by discarding high moments creates the same basic realizability
problem: the retained moments may not correspond to any probability
distribution.

Some classical closures are designed to avoid this problem. [Levermore
(1996)](https://link.springer.com/article/10.1007/BF02179552) constructs closures
that preserve realizability, hyperbolicity, and an entropy. For circular
cumulants of coupled phase oscillators, [Goldobin and Dolmatova
(2019)](https://journals.aps.org/prresearch/abstract/10.1103/PhysRevResearch.1.033139)
show that the Ott-Antonsen ansatz is the only admissible exact finite
truncation in the class they consider. These results do not supply a closure for
quantum operator moments; they show why a useful closure may need structural
constraints beyond a higher truncation order.

### A rigorous result for a different expansion

[Chen et al. (2025)](https://arxiv.org/abs/2512.12010) prove convergence and
polynomial runtime for an expansion of the equilibrium log-partition function
of weakly interacting fermions. The proof relies on an explicit weak-coupling
condition; at any temperature, the resulting algorithm is polynomial in both
system size and precision. It computes a thermodynamic quantity rather than the
time evolution of operator moments.

This result demonstrates that an expansion in connected interaction
contributions can be controlled when the expanded quantity and the small
parameter are specified precisely. It does not establish convergence for
open-system dynamics, bosonic modes, spins, or the closure used by
QuantumCumulants.jl.

### A trajectory-based extension

The truncated-cumulant-trajectory method applies a cumulant closure separately
along stochastic quantum trajectories and then averages over the trajectories
([Verstraelen, Huybrechts, Roscilde and Wouters,
2023](https://journals.aps.org/prxquantum/abstract/10.1103/PRXQuantum.4.030304)).
The correlator hierarchy is cut at an order ``k_c`` on each trajectory.
Fluctuations between trajectories can then generate correlations that are
absent from a single deterministic closure.

The method performed well for the dissipative spin lattices studied by the
authors, but it remains approximate and can develop unstable trajectories. It
does not provide a general positivity or convergence guarantee and is not
currently implemented by QuantumCumulants.jl.

### Other closure principles

A [maximum-entropy
closure](https://en.wikipedia.org/wiki/Principle_of_maximum_entropy) reconstructs
a state or distribution by maximizing its entropy subject to the retained
moments. Reduced-density-matrix closures instead evolve low-body density
matrices and approximate the higher-body matrices, often while enforcing
representability conditions. Model-specific closures may use a controlled
``1/N`` expansion, a symmetry, or another physical constraint.

These approaches can preserve properties that direct cumulant truncation does
not, but they require different variables, constraints, or numerical steps.
They are not interchangeable variants of the current `order` keyword and are
not currently implemented by QuantumCumulants.jl.

## Open problems

Important unresolved questions include:

* developing predictive criteria and rigorous error bounds that determine when
  dynamical cumulant closures converge, are asymptotic, or fail in
  driven-dissipative bosonic and spin systems; and
* constructing useful quantum closures that preserve realizability or
  positivity while retaining the computational advantages of moment equations.

## References

* J. Marcinkiewicz.
  [Sur une propriété de la loi de Gauß. *Mathematische Zeitschrift* **44**, 612--618 (1939)](https://eudml.org/doc/168829).
* R. Kubo.
  [Generalized cumulant expansion method. *Journal of the Physical Society of Japan* **17**, 1100 (1962)](https://www.jstage.jst.go.jp/article/jpsj1946/17/7/17_7_1100/_article/-char/ja/).
* A. K. Rajagopal and E. C. G. Sudarshan.
  [Some generalizations of the Marcinkiewicz theorem and its implications to certain approximation schemes in many-particle physics. *Physical Review A* **10**, 1852 (1974)](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.10.1852).
* P. Hänggi and P. Talkner.
  [A remark on truncation schemes of cumulant hierarchies. *Journal of Statistical Physics* **22**, 65--67 (1980)](https://link.springer.com/article/10.1007/BF01007988).
* M. Kira and S. W. Koch.
  [Cluster-expansion representation in quantum optics. *Physical Review A* **78**, 022102 (2008)](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.78.022102).
* M. Bianucci and M. Bologna.
  [About the foundation of the Kubo generalized cumulants theory: a revisited and corrected approach. *Journal of Statistical Mechanics* (2020) 043405](https://arxiv.org/abs/1911.09620).
* D. Plankensteiner, C. Hotter and H. Ritsch.
  [QuantumCumulants.jl: A Julia framework for generalized mean-field equations in open quantum systems. *Quantum* **6**, 617 (2022)](https://doi.org/10.22331/q-2022-01-04-617).
* O. Rubies-Bigorda, S. Ostermann and S. F. Yelin.
  [Characterizing superradiant dynamics in atomic arrays via a cumulant expansion approach. *Physical Review Research* **5**, 013091 (2023)](https://doi.org/10.1103/PhysRevResearch.5.013091).
* M. Fasser, C. Genes, H. Ritsch and R. Holzinger.
  [Benchmarking the cumulant expansion method using Dicke superradiance. *Photonics* **12**, 996 (2025)](https://doi.org/10.3390/photonics12100996).
* P. Fowler-Wright, K. B. Arnardóttir, P. Kirton, B. W. Lovett and J. Keeling.
  [Determining the validity of cumulant expansions for central spin models. *Physical Review Research* **5**, 033148 (2023)](https://journals.aps.org/prresearch/abstract/10.1103/PhysRevResearch.5.033148).
* W. Verstraelen, D. Huybrechts, T. Roscilde and M. Wouters.
  [Quantum and classical correlations in open quantum-spin lattices via truncated-cumulant trajectories. *PRX Quantum* **4**, 030304 (2023)](https://journals.aps.org/prxquantum/abstract/10.1103/PRXQuantum.4.030304).
* C. D. Levermore.
  [Moment closure hierarchies for kinetic theories. *Journal of Statistical Physics* **83**, 1021 (1996)](https://link.springer.com/article/10.1007/BF02179552).
* D. S. Goldobin and A. V. Dolmatova.
  [Ott-Antonsen ansatz truncation of a circular cumulant series. *Physical Review Research* **1**, 033139 (2019)](https://journals.aps.org/prresearch/abstract/10.1103/PhysRevResearch.1.033139).
* J. Kerber, H. Ritsch and L. Ostermann.
  [The cumulants expansion approach: the good, the bad and the ugly. arXiv:2511.20115 (2025)](https://arxiv.org/abs/2511.20115).
* H. Chen et al.
  [Convergence of the cumulant expansion and polynomial-time algorithm for weakly interacting fermions. arXiv:2512.12010 (2025)](https://arxiv.org/abs/2512.12010).
