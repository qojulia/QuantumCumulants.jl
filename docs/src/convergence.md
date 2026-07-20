# [Convergence and validity of the cumulant expansion](@id convergence)

The [cumulant expansion](@ref theory) turns an open hierarchy of operator
averages into a finite set of equations by neglecting cumulants above a chosen
order ``M``. This raises two practical questions: does increasing ``M`` improve
the result, and when can a truncation be trusted?

The short answer is that cumulants provide an exact reformulation when the full
hierarchy is retained, but truncation is uncontrolled in general. A higher order
includes more correlations, yet is not guaranteed to be more accurate. The
method becomes controlled only when an additional argument—such as Gaussianity
or a well-defined ``1/N`` expansion—makes the neglected cumulants small.

This does not make cumulant closure a bad approximation. It is often extremely
useful and becomes rigorous in identifiable limits. The important distinction is
between demonstrated accuracy for a particular problem and a general guarantee
that the method does not possess.

## Practical guidance for QuantumCumulants.jl

!!! tip "How to use the truncation order responsibly"
    * **Check convergence empirically.** Solve at consecutive affordable orders
      and compare the observables of interest. Agreement between successive
      orders is useful evidence, although it is not an error bound.
    * **Validate against an exact reference when possible.** Compare a small
      instance to a full master-equation solution, for example with
      [QuantumOptics.jl](https://qojulia.org/), before scaling up.
    * **Watch for unphysical results.** Negative populations, expectation values
      outside their allowed range, or runaway oscillations may indicate that the
      truncated equations no longer describe a physical state.
    * **Do not assume that higher is better.** Increasing the order can degrade
      accuracy or destabilize the equations.
    * **Use second order when Gaussian correlations dominate.** It is exact for
      Gaussian states under Gaussian-preserving dynamics and is often a good
      accuracy-for-cost compromise, but it is not automatically physical or
      accurate for general dynamics.
    * **Check the system-size scaling.** Large-``N`` collective models are an
      important controlled regime when their parameter scaling produces a
      genuine ``1/N`` expansion. Large ``N`` alone is not sufficient.

## The expansion and its truncation

The transformation between moments and cumulants is triangular and invertible:
the complete cumulant hierarchy contains the same information as the complete
moment hierarchy. For a finite-dimensional operator algebra, a complete set of
moments determines the state, so retaining the full hierarchy reproduces the
exact dynamics. For bosonic systems with an infinite-dimensional Hilbert space,
the corresponding statement also requires the usual assumptions that make the
moment problem determinate.

For ``N`` finite-level subsystems, local operator products can be reduced to a
finite operator basis. In an atom-only model of ``N`` two-level systems, for
example, order ``M=N`` can retain every inter-atomic correlation and reproduce
the full dynamics. This does not apply to an untruncated bosonic mode, whose
operator hierarchy remains infinite.

The practical question is therefore not whether the full hierarchy is exact,
but whether the affordable intermediate orders approach the desired result.
There is no general affirmative answer.

## The order ``M`` is not a small parameter

The truncation order is not the power of a coupling constant. Unlike a
perturbative expansion, an order-``M`` calculation does not contain an explicit
factor that forces the omitted order-``M+1`` cumulants to be smaller. For a
strongly non-Gaussian state, high-order cumulants need not decay rapidly with
order.

Increasing ``M`` therefore adds correlations without, by itself, guaranteeing a
smaller error. A convergence argument needs some independent smallness
mechanism, such as weak non-Gaussianity or suitable scaling with system size.

## Gaussian closure and realizability

There is a precise reason why second order is special. For a classical random
variable ``X`` with characteristic function
``\varphi(t)=\langle e^{itX}\rangle``, the cumulant generating function is

```math
K(t)=\log\varphi(t)=\sum_{n\geq 1}\kappa_n\frac{(it)^n}{n!}.
```

Marcinkiewicz's theorem states that if ``K(t)`` is a polynomial, its degree is at
most two. Thus a non-Gaussian probability distribution cannot have an exactly
terminating cumulant sequence. Rajagopal and Sudarshan extended this obstruction
to many variables and to commuting and anticommuting operator-valued variables.
In the settings covered by these results, an exactly terminating hierarchy
describes either a Gaussian object, with cumulants through second order, or a
special first-order case such as a coherent state; a genuinely non-Gaussian
state generally has cumulants at arbitrarily high orders.

This result has two consequences for cumulant closure:

1. A second-order closure has a consistent interpretation as a Gaussian ansatz.
   A quadratic Hamiltonian with linear dissipation gives closed equations for
   first and second moments. Those moments evolve exactly for any initial state;
   the entire state is described exactly when the initial state is Gaussian and
   the dynamics preserves Gaussianity.

2. Setting all cumulants above a finite order ``M>2`` to zero does not, in
   general, define the exact cumulants of a physical non-Gaussian state. This is
   a structural warning about realizability, not a statement that every
   higher-order numerical solution must become unphysical. Conversely, a
   second-order closure is not by itself a guarantee of positivity.

The practical symptoms of a realizability failure can include negative
populations, correlations outside their allowed bounds, and unstable or runaway
solutions.

!!! note "Operator cumulants"
    For noncommuting operators, generalized cumulants require care about operator
    ordering and about which properties of classical cumulants survive. The
    closure used by QuantumCumulants.jl is based on the joint-cumulant definition
    given in [Theoretical background](@ref theory), rather than on naive moment
    factorization. Bianucci and Bologna (2020) discuss corrections and
    qualifications to Kubo's original operator treatment.

## No general convergence theorem

Mean field is the lowest-order cumulant closure, but neither it nor the sequence
of higher closures has a general convergence theorem. Applicability depends on
the Hamiltonian, dissipation, initial state, observables, and the scaling of
parameters with system size. In particular, agreement at one time or for one
observable does not guarantee agreement elsewhere.

As Kerber, Ritsch and Ostermann (2025) emphasize, a general criterion for the
applicability and convergence of higher-order cumulant expansions remains to be
found. The rigorous results that do exist apply to specific regimes and should
not be promoted to a guarantee for general driven-dissipative dynamics.

### Asymptotic or convergent?

Even a sequence that initially improves may be asymptotic rather than
convergent. It can improve up to some order and then deteriorate; even and odd
orders can also behave differently. Unless a model supplies a genuine small
parameter or a separate convergence result, order-by-order agreement should be
treated as an empirical diagnostic rather than a proof.

The central-spin results discussed below give a concrete warning: if even and
odd orders approach different limits, the sequence has no single limit at all.
It is therefore safest to regard the expansion as potentially asymptotic unless
a convergence argument applies to the model at hand.

## What benchmarks show

A systematic benchmark by Kerber, Ritsch and Ostermann (2025) found behavior
spanning the “good,” the “bad,” and the “ugly”:

* **Good.** For collective radiative decay in a dipole-dipole interacting atom
  chain, the tested truncations improved smoothly with order. The gains beyond
  second order were modest relative to their computational cost. Order-by-order
  comparisons in ordered atomic arrays reached a compatible conclusion
  (Rubies-Bigorda, Ostermann and Yelin, 2023).
* **Bad.** In a model with higher-body interaction terms, an intermediate
  truncation could agree less well with the exact result than a lower-order one.
* **Ugly.** For a bi-prime factorization problem formulated as adiabatic
  annealing, going beyond mean field did not improve the useful output. Even for
  a small system, intermediate-order equations developed fast oscillations,
  divergences, and partly non-physical solutions.

These labels describe the examples in that benchmark, not a classification
theorem for Hamiltonians. They illustrate why the behavior of one model cannot
be transferred automatically to another.

### Non-monotonic and non-uniform convergence: the central-spin lesson

Central-spin models provide a particularly sharp caution (Fowler-Wright,
Arnardóttir, Kirton, Lovett and Keeling, 2023):

1. The common expectation that mean field captures the exact ``N\to\infty``
   limit is false in general. Whether it does so depends on how the model
   parameters scale with ``N``.

2. Convergence can be non-uniform across even and odd orders. The two
   subsequences can approach distinct limits rather than improving monotonically
   toward one result.

3. The error need not decrease monotonically with ``N``. Even when a
   higher-order expansion recovers the correct large-``N`` limit, it can be less
   accurate than mean field at intermediate sizes.

A benchmark of Dicke superradiance against an exact solution gives a compatible
but observable-dependent picture (Fasser, Genes, Ritsch and Holzinger, 2025).
The time and magnitude of peak emission converged reliably with order, whereas
long-time populations were captured only for small emitter numbers. Odd orders
showed unphysical behavior, and at large emitter number neither the
individual-spin nor the collective-spin expansion converged to the exact result.
In that study, the individual-spin expansion was more reliable than the
collective-spin variant.

The practical conclusion is deliberately limited: higher order is not
guaranteed to improve accuracy, either as a function of truncation order or as a
function of system size.

## Regimes with additional control

| Regime | Source of control | Limitation |
|---|---|---|
| Gaussian states under Gaussian-preserving dynamics | Cumulants above second order vanish exactly | Does not cover dynamics that generates non-Gaussian correlations |
| Near-Gaussian dynamics | Higher cumulants remain small over the relevant time interval | Smallness must be justified or checked for the model |
| Large ``N`` or high connectivity | A systematic ``1/N`` expansion may suppress connected correlations | Requires appropriate parameter scaling; large ``N`` alone is insufficient |

There are rigorous convergence results for other kinds of cumulant expansions.
For example, Chen et al. (2025) prove convergence and polynomial runtime for a
linked-cluster expansion of the equilibrium log-partition function of weakly
interacting fermions. Their randomized algorithm is polynomial in both system
size and precision under an explicit weak-coupling condition at any temperature.
That is an important related result, but it concerns an equilibrium
linked-cluster expansion and is not a proof for the driven-dissipative moment
closure implemented by QuantumCumulants.jl.

## Related approaches

The operator hierarchy is the quantum analogue of the classical BBGKY hierarchy,
and closing it by neglecting high cumulants is an instance of the broader
moment-closure problem. Two closely related issues arise in the classical
setting:

* **Realizability.** A finite set of approximate moments or cumulants need not
  correspond to any probability distribution. Levermore's moment-closure
  hierarchies (1996) were designed to preserve realizability, hyperbolicity, and
  an entropy rather than closing the hierarchy by direct truncation.
* **Admissible closures can be highly restricted.** For circular cumulants of
  coupled phase oscillators, Goldobin and Dolmatova (2019) showed that the
  Ott-Antonsen ansatz is the only admissible exact finite truncation in the class
  they considered.

These classical results do not directly determine the validity of a quantum
closure. They demonstrate that retaining more moments is not, by itself, the
same as preserving the set of realizable states.

### Truncated-cumulant trajectories

For open quantum-spin systems, truncated-cumulant trajectories combine a
cumulant closure along stochastic quantum trajectories with an ensemble average
(Verstraelen, Huybrechts, Roscilde and Wouters, 2023). The correlator hierarchy
is truncated at an order ``k_c`` on each trajectory, while fluctuations between
trajectories contribute additional correlations to the ensemble result.

This construction can recover correlations that are absent from a single
deterministic closure and performed well for the dissipative spin lattices
studied by the authors. It should not, however, be interpreted as a general proof
of positivity or convergence: the method remains approximate and can develop
unstable trajectories. It is not currently part of QuantumCumulants.jl, but is a
useful example of how cumulant closure can be combined with a different
representation of open-system dynamics.

## Open problems

Important unresolved questions include:

* finding an a priori criterion, based on interaction structure, symmetry,
  near-Gaussianity, or parameter scaling with ``N``, that predicts whether a
  given dynamical closure converges, is merely asymptotic, or fails;
* determining whether rigorous convergence techniques for equilibrium
  linked-cluster expansions can be adapted to non-equilibrium and
  driven-dissipative dynamics, and to bosonic and spin systems; and
* constructing useful quantum closures that preserve realizability or
  positivity while retaining the computational advantages of moment equations.

## References

* J. Marcinkiewicz.
  [Sur une propriété de la loi de Gauß. *Mathematische Zeitschrift* **44**, 612 (1939)](https://link.springer.com/article/10.1007/BF01210677).
* R. Kubo.
  [Generalized cumulant expansion method. *Journal of the Physical Society of Japan* **17**, 1100 (1962)](https://www.jstage.jst.go.jp/article/jpsj1946/17/7/17_7_1100/_article/-char/ja/).
* A. K. Rajagopal and E. C. G. Sudarshan.
  [Some generalizations of the Marcinkiewicz theorem and its implications to certain approximation schemes in many-particle physics. *Physical Review A* **10**, 1852 (1974)](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.10.1852).
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
