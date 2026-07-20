# [Convergence and validity of the cumulant expansion](@id convergence)

The [cumulant expansion](@ref theory) turns the infinite hierarchy of operator
averages into a finite, closed set by discarding cumulants above a chosen order
``M``. A natural question follows immediately: as we raise ``M``, do the results
converge to the exact dynamics, and how trustworthy is a given truncation? This
page collects what is known.

The short answer is that the expansion is *exact in principle* and *often
excellent in practice*, but as a general-purpose method it is *uncontrolled*.
There is no general theorem guaranteeing that a higher order is more accurate,
and no a priori criterion that tells you in advance whether a given truncation
can be trusted. The method becomes provably controlled only in special regimes.
Calling it a "bad truncation" is too blunt. A more accurate statement is
"uncontrolled in general, very useful in practice, and rigorous in identifiable
limits".

## The expansion is exact; truncation is the only approximation

The map between ordinary averages (moments) and cumulants is a triangular
polynomial transform, hence invertible: the full set of cumulants carries exactly
the same information as the full set of averages, and therefore as the density
matrix itself. Keeping *all* clusters reproduces the exact dynamics; the only
approximation is the act of dropping cumulants above order ``M``. In this sense
the cumulant (or *cluster*) expansion is an exact reformulation of the many-body
problem, not an approximation in itself.

A useful consequence: for a finite system the hierarchy terminates, so a
sufficiently large ``M`` is exact. For ``N`` two-level atoms, for instance,
truncating at ``M = N`` drops nothing. The practical question is therefore never
"does ``M`` large enough work" (it does, trivially), but "do the *intermediate*,
affordable orders improve monotonically toward the answer". That is the question
with no general answer.

## The order ``M`` is not a small parameter

The crucial conceptual point is that the truncation order ``M`` is *not* the power
of any coupling constant. Unlike a perturbative expansion, where the ``n``-th term
carries an explicit factor of a small parameter, nothing forces cumulants of
order ``M+1`` to be smaller than those of order ``M``. For states that are far
from Gaussian the higher cumulants need not decay at all. Increasing ``M`` adds
more correlations, but without a smallness mechanism it does not, by itself,
guarantee a better answer. This is the structural reason the expansion is
uncontrolled in general.

## The Marcinkiewicz obstruction

There is a precise theorem behind the special role of second order.

For a random variable ``X`` with characteristic function
``\varphi(t) = \langle e^{itX}\rangle``, the cumulant generating function is
``K(t) = \log\varphi(t) = \sum_{n\ge 1}\kappa_n\,(it)^n/n!``, whose Taylor
coefficients ``\kappa_n`` are the cumulants. **Marcinkiewicz's theorem** (1939)
states that if ``K(t)`` is a polynomial, then its degree is at most ``2``.
Equivalently, the Gaussian is the only distribution with finitely many nonzero
cumulants; every non-Gaussian distribution has infinitely many nonzero cumulants.

Rajagopal and Sudarshan (1974) generalized this result to many variables and to
both commuting (Bose) and anticommuting (Fermi) operator-valued variables. For
quantum states the upshot is that a physical state has either

* only a nonzero first-order cumulant (a coherent state),
* nonzero cumulants only up to second order (a Gaussian state), or
* nonzero cumulants to all orders.

No physical state has cumulants that vanish identically above some finite order
``M > 2``. Two consequences follow directly.

1. **Order two is special.** The Gaussian is the unique consistent finite
   truncation. This is why a purely quadratic Hamiltonian with linear dissipation
   closes exactly at order ``2`` (see [Theoretical background](@ref theory)), and
   why a second-order expansion is exact for any Gaussian-preserving dynamics.

2. **Finite truncation breaks positivity.** Setting all cumulants above some
   ``M > 2`` to zero is structurally inconsistent with the existence of a genuine
   state. Rajagopal and Sudarshan already noted that such truncation schemes can
   exhibit non-positive behavior, meaning the truncated description need not
   correspond to a valid density matrix. This is the root of the unphysical
   outputs that a truncated expansion can produce, such as negative populations,
   correlations outside their allowed range, or runaway oscillations.

!!! note "Why the code uses cumulants, not just factorized moments"
    Because operators do not commute, there is in general no closed-form formula
    giving cumulants from moments (only the reverse direction has one), so the
    expansion is built from the joint cumulant definition rather than naive
    factorization. See Bianucci and Bologna (2020) for the corrected operator
    treatment of Kubo's original formulation.

## No general convergence theorem and no a priori criterion

Despite the ubiquity of the method (mean field is simply its lowest order), there
is no general convergence theorem and no general criterion for when a higher order
helps. As the most recent dedicated study puts it: "a general criterion for
applicability and convergence properties of higher order cumulants expansions
remains to be found" (Kerber, Ritsch and Ostermann, 2025). Convergence depends on
the system, the state, and on how parameters scale with system size. The
rigorous results that *do* exist are narrowly regime-specific, which by itself
confirms that no general criterion is available.

## Asymptotic or convergent?

Even where low orders improve, it is rarely known whether the series *converges*
or is merely *asymptotic*, that is, improves up to some order and then degrades.
The only setting with a proven convergent series is weakly interacting fermions
(see below); in general there is no such guarantee. The central-spin result is a
concrete warning: when even and odd orders approach different limits, the sequence
of truncations has no single limit at all, so "convergent" is simply the wrong
word for it. A safe default is to treat the expansion as asymptotic unless a
convergence proof or a genuine small parameter (such as ``1/N`` or weak coupling)
applies.

## What happens in practice: the Good, the Bad and the Ugly

A systematic benchmark (Kerber, Ritsch and Ostermann, 2025) found the observed
behavior spans a spectrum:

* **Good.** For the collective radiative decay of a dipole-dipole interacting
  atom chain, the truncation converges smoothly, with monotone improvement as the
  order increases. Even here the gains beyond second order are modest, so second
  order is usually the best accuracy-for-cost choice. Order-by-order comparisons
  in ordered atomic arrays reach the same conclusion (Rubies-Bigorda, Ostermann
  and Yelin, 2023).
* **Bad.** Higher-order interaction terms in the Hamiltonian can make a
  higher-order truncation agree *less* with the exact solution than a lower-order
  one.
* **Ugly.** For a bi-prime factorization problem solved by adiabatic annealing,
  going beyond mean field is useless, and even at small system size the equations
  become numerically pathological and partly non-physical.

## Non-monotonic and non-uniform convergence: the central-spin lesson

The sharpest cautionary results come from central-spin (many-to-one connectivity)
models (Fowler-Wright, Arnardóttir, Kirton, Lovett and Keeling, 2023). Three
findings overturn common expectations:

1. The widespread assumption that mean field captures the exact ``N\to\infty``
   limit, and that higher cumulant orders converge to that same limit while
   improving finite-``N`` accuracy, is **false in general**. Whether mean field is
   even exact at large ``N`` depends on how the model parameters scale with ``N``.

2. **Convergence can be non-uniform across even and odd orders.** Even-order and
   odd-order truncations can approach *distinct* limits, so the expansion does not
   converge monotonically order by order.

3. **The error need not be monotonic in ``N``.** Even when a higher-order
   expansion does recover the correct ``N\to\infty`` limit, its error as a
   function of ``N`` is non-monotonic and can exceed that of plain mean field.

The practical takeaway is blunt: a higher order is not guaranteed to be more
accurate, neither in the truncation order nor in the system size.

A benchmark of Dicke superradiance against an exact solution (Fasser, Genes,
Ritsch and Holzinger, 2025) reaches a compatible verdict in the collective
setting: time and magnitude of the peak emission converge reliably with order,
but the long-time populations are only captured for small emitter numbers, odd
orders show unphysical behavior, and at large emitter number both an
individual-spin and a collective-spin cumulant expansion fail to converge to the
exact result. There too an individual-spin expansion proves more reliable than
the collective-spin variant.

## Connection to classical moment closure and realizability

None of this is unique to quantum mechanics. The operator hierarchy is the
quantum analog of the classical BBGKY hierarchy, and closing it by discarding high
cumulants is exactly the *moment-closure* problem of kinetic theory and
statistical physics. The same two difficulties appear there:

* **Realizability.** A truncated set of moments or cumulants may not correspond to
  *any* genuine probability distribution. This is the classical mirror of the
  positivity violation above, and it is precisely the content of the Marcinkiewicz
  theorem: among finite truncations only the Gaussian closure (second order) is
  guaranteed realizable. Levermore's moment-closure hierarchies (1996) were
  constructed specifically to preserve realizability, hyperbolicity, and an
  entropy, rather than truncating naively.
* **Only special truncations are admissible.** In some systems the set of valid
  finite truncations is severely restricted. For the circular-cumulant series of
  coupled phase oscillators, Goldobin and Dolmatova (2019) proved that the
  Ott-Antonsen ansatz is the *only* admissible truncation, with every other finite
  truncation forbidden, much as Marcinkiewicz forbids finite non-Gaussian cumulant
  truncations.

The lesson carried over from the classical literature is that a good closure is
one that is realizability-preserving by construction, not one that merely keeps
more terms.

## When the expansion is actually controlled

The method is trustworthy, and in some cases provably so, in the following
regimes.

| Regime | Why it is controlled | Caveat |
|---|---|---|
| Near-Gaussian dynamics, order ``2`` | Gaussian is the unique consistent finite truncation; exact for quadratic ``H`` with linear dissipation | Fails once genuinely non-Gaussian correlations build up |
| Weak coupling | Rigorous convergence proof exists for weakly interacting fermions | Proven for the equilibrium log-partition function, not yet for general driven-dissipative dynamics |
| Large ``N`` / high connectivity | Mean field plus a ``1/N`` expansion; ``1/N`` is a genuine small parameter | Only when the parameter scaling with ``N`` cooperates (see the central-spin caveat above) |

The one fully rigorous "converges and runs in polynomial time" guarantee proven
to date is for weakly interacting fermions at equilibrium (Chen et al., 2025):
the cumulant series for the log-partition function converges, with a randomized
algorithm whose runtime is polynomial in both system size and precision, valid
under an explicit weak-coupling (intensive interaction strength) condition at any
temperature. This is narrower than the general open-system operator truncation
used here, but it is the kind of theorem the field had been missing.

## Restoring positivity: truncated-cumulant trajectories

Because a single deterministic truncation can leave the realizable set, one remedy
is to avoid imposing one global closure. The truncated-cumulant-trajectory method
(Verstraelen, Huybrechts, Roscilde and Wouters, 2023) combines stochastic quantum
trajectories with a cumulant truncation applied *along each trajectory*: the
correlator hierarchy is cut at a finite cumulant order ``k_c`` on every trajectory,
and the ensemble average over trajectories restores correlations beyond that order.
Because the full state is never forced into a fixed density-matrix form, positivity
is respected by the dynamics rather than imposed by hand, which mitigates the
unphysical solutions that a single deterministic truncation can produce. The
authors invoke the Marcinkiewicz theorem to motivate the Gaussian (second-order)
level as the natural per-trajectory truncation. This scheme is not part of
QuantumCumulants.jl, but it is a useful template for what to do when a plain
truncation turns non-physical.

## Practical guidance for QuantumCumulants.jl

!!! tip "How to use the truncation order responsibly"
    * **Check convergence empirically.** Derive and solve at consecutive orders
      (for example ``M = 1, 2, 3``) and compare. Treat agreement between
      successive orders as the diagnostic, since no a priori bound is available.
    * **Validate against an exact reference when possible.** For a small instance,
      compare to a full master-equation solution (for example via
      [QuantumOptics.jl](https://qojulia.org/)) before scaling up.
    * **Watch for unphysical results.** Negative populations, expectation values
      outside their allowed range, or runaway oscillations are red flags that the
      truncation has broken positivity, not numerical noise to be smoothed over.
    * **Second order is privileged.** It is the highest order consistent with a
      genuine quantum state, it is exact for Gaussian-preserving dynamics, and it
      is usually the best accuracy-for-cost tradeoff.
    * **Do not assume higher is better.** Increasing the order can degrade
      accuracy or destabilize the equations. Verify, do not extrapolate.
    * **The natural home of the method is large ``N`` with collective coupling**,
      where ``1/N`` provides genuine control, but still check how your parameters
      scale with ``N``.

## Open problems

* A general a priori criterion (in terms of interaction structure, symmetry,
  parameter scaling with ``N``, or near-Gaussianity) that predicts whether a
  given expansion will converge, be merely asymptotic, or fail, rather than
  diagnosing it model by model.
* Whether the rigorous weak-coupling fermion convergence proof extends to
  non-equilibrium and driven-dissipative dynamics, to bosonic and spin systems,
  and beyond the weak-coupling condition.
* A unified realizability condition (a valid-state or positivity criterion)
  linking the even/odd and non-monotonic-in-``N`` behavior of quantum truncations
  to the classical moment-closure literature, and ideally a positivity-preserving
  truncation that exploits the Marcinkiewicz constraint.

## References

* J. Marcinkiewicz.
  [Sur une propriété de la loi de Gauß. *Mathematische Zeitschrift* **44**, 612 (1939)](https://link.springer.com/article/10.1007/BF01210677).
* R. Kubo.
  [Generalized cumulant expansion method. *Journal of the Physical Society of Japan* **17**, 1100 (1962)](https://www.jstage.jst.go.jp/article/jpsj1946/17/7/17_7_1100/_article/-char/ja/).
* A. K. Rajagopal and E. C. G. Sudarshan.
  [Some generalizations of the Marcinkiewicz theorem and its implications to certain approximation schemes in many-particle physics. *Physical Review A* **10**, 1852 (1974)](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.10.1852).
* M. Kira and S. W. Koch.
  [Cluster-expansion representation in quantum optics. *Physical Review A* **78**, 022102 (2008)](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.78.022102).
  See also M. Kira and S. W. Koch, *Semiconductor Quantum Optics* (Cambridge
  University Press, 2011).
* M. Bianucci and M. Bologna.
  [About the foundation of the Kubo Generalized Cumulants theory. A revisited and corrected approach. *Journal of Statistical Mechanics* (2020) 043405](https://arxiv.org/abs/1911.09620).
* D. Plankensteiner, C. Hotter and H. Ritsch.
  [QuantumCumulants.jl: A Julia framework for generalized mean-field equations in open quantum systems. *Quantum* **6**, 617 (2022)](https://doi.org/10.22331/q-2022-01-04-617).
* O. Rubies-Bigorda, S. Ostermann and S. F. Yelin.
  [Characterizing superradiant dynamics in atomic arrays via a cumulant expansion approach. *Physical Review Research* **5**, 013091 (2023)](https://doi.org/10.1103/PhysRevResearch.5.013091).
* M. Fasser, C. Genes, H. Ritsch and R. Holzinger.
  [Benchmarking the Cumulant Expansion Method Using Dicke Superradiance. *Photonics* **12**, 996 (2025)](https://doi.org/10.3390/photonics12100996).
* P. Fowler-Wright, K. B. Arnardóttir, P. Kirton, B. W. Lovett and J. Keeling.
  [Determining the validity of cumulant expansions for central spin models. *Physical Review Research* **5**, 033148 (2023)](https://journals.aps.org/prresearch/abstract/10.1103/PhysRevResearch.5.033148).
* W. Verstraelen, D. Huybrechts, T. Roscilde and M. Wouters.
  [Quantum and classical correlations in open quantum-spin lattices via truncated-cumulant trajectories. *PRX Quantum* **4**, 030304 (2023)](https://journals.aps.org/prxquantum/abstract/10.1103/PRXQuantum.4.030304).
* C. D. Levermore.
  [Moment closure hierarchies for kinetic theories. *Journal of Statistical Physics* **83**, 1021 (1996)](https://link.springer.com/article/10.1007/BF02179552).
* D. S. Goldobin and A. V. Dolmatova.
  [Ott-Antonsen ansatz truncation of a circular cumulant series. *Physical Review Research* **1**, 033139 (2019)](https://journals.aps.org/prresearch/abstract/10.1103/PhysRevResearch.1.033139).
* J. Kerber, H. Ritsch and L. Ostermann.
  [The Cumulants Expansion Approach: The Good, The Bad and The Ugly. arXiv:2511.20115 (2025)](https://arxiv.org/abs/2511.20115).
* Chen et al.
  [Convergence of the Cumulant Expansion and Polynomial-Time Algorithm for Weakly Interacting Fermions. arXiv:2512.12010 (2025)](https://arxiv.org/abs/2512.12010).
