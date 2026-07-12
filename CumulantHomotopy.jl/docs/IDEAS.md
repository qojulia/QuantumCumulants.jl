# Controlling stationary roots and incorporating structure preserving reductions

## 1. Objective

At truncation order (n), QuantumCumulants produces a finite nonlinear stationary system

[
F_n(x)=0,
]

where (x) contains the retained moments or cumulants.

Even when the original Liouvillian has a unique stationary density operator, the truncated system may have many algebraic roots. CumulantHomotopy should not attempt to enumerate every root at high order. Instead, it should construct, track, validate, and interpret a bounded set of physically motivated stationary candidates.

The central objective is therefore not

[
#{x:F_n(x)=0}=1.
]

The relevant objective is that the set of physically admissible candidates collapses in its low order observables:

[
\operatorname{diam}
\left(
\pi_k\mathcal P_n
\right)
\longrightarrow0
]

for every fixed projection order (k), where (\mathcal P_n) denotes the physically plausible stationary candidates at order (n).

This allows formal algebraic roots to remain while requiring all physically admissible candidates to produce the same low order stationary observables.

## 2. Controlling the number of tracked roots

Let

[
\mathcal R_n
============

{x:F_n(x)=0}
]

denote the complete algebraic root set.

At high hierarchy order, CumulantHomotopy should maintain only a bounded active candidate set

[
\mathcal C_n\subseteq\mathcal R_n,
\qquad
|\mathcal C_n|\leq B,
]

where (B) is a user configurable branch budget.

The branch budget prevents the method from scaling with the Bézout or mixed volume root count.

### 2.1 Primary branch

The first candidate should be connected to a known physical reference state.

Possible references include:

1. A noninteracting or Gaussian limit.
2. A weak drive or weak coupling limit.
3. A finite Hilbert space stationary state.
4. A maximum entropy reconstruction.
5. A phase space stationary approximation.
6. A reset regularized Liouvillian.

A useful reset path is

[
\mathcal L_{\lambda,\epsilon}
=============================

\lambda\mathcal L_{\mathrm{target}}
+
\epsilon\mathcal R_\sigma,
]

with

[
\mathcal R_\sigma(\rho)
=======================

\operatorname{Tr}(\rho)\sigma-\rho.
]

At (\lambda=0), the reference state is exactly (\sigma). The target Liouvillian is activated while the reset term regularizes the stationary problem, after which the reset strength is removed.

### 2.2 Branch spawning

Additional branches should only be generated when the current path becomes ambiguous.

Possible branch spawning conditions include:

[
\sigma_{\min}(J_n)<\varepsilon_J,
]

large Newton corrections, repeated corrector failures, abrupt changes in physical diagnostics, disagreement between different continuation paths, or failure of forward and reverse tracking.

Near an ambiguous point, new candidates can be generated using:

1. Singular vectors of the Jacobian.
2. Deflated Newton corrections.
3. Small perturbations normal to the current branch tangent.
4. Alternative physical parameter paths.
5. Alternative hierarchy lift predictors.
6. Phase space identified metastable states.

### 2.3 Candidate pruning

A branch should be removed from the active set when it has one or more of the following properties:

[
\begin{aligned}
&\text{large stationary residual},\
&\text{large adjoint inconsistency},\
&\text{strong moment matrix negativity},\
&\text{violation of symmetry or conserved quantities},\
&\text{unphysical occupation or energy growth},\
&\text{large held out Liouvillian residual},\
&\text{large correction under physical purification},\
&\text{failure to continue to the next hierarchy order}.
\end{aligned}
]

Pruning should be lexicographic rather than based on a single weighted score. Strong physical violations should remove a candidate before weaker ranking criteria are considered.

### 2.4 Root clustering

Algebraically distinct roots may be physically indistinguishable in low order observables.

For candidates (x_i) and (x_j), define

[
d_k(x_i,x_j)
============

\left|
S_k
\left(
\pi_kx_i-\pi_kx_j
\right)
\right|.
]

Candidates satisfying

[
d_k(x_i,x_j)<\varepsilon_k
]

should be placed in the same physical cluster.

The important quantity is therefore

[
N_{n,k}^{\mathrm{clusters}},
]

rather than the raw number of polynomial roots.

A practical definition of hierarchy collapse is

[
N_{n,k}^{\mathrm{clusters}}=1
]

for all requested low order projections (k), together with decreasing order lift defects and improving physical diagnostics.

### 2.5 Local root certification

Once a candidate root is obtained, interval Newton or Krawczyk methods can certify that exactly one root lies in a neighborhood of the numerical solution.

This does not prove global uniqueness, but it prevents numerical duplicates from being counted as independent branches and provides a stronger endpoint certificate.

## 3. Metastable branches and stationary mixtures

Multiple physically plausible cumulant roots need not all be spurious.

In systems with semiclassical bistability or multistability, different roots may approximate distinct metastable basins. The exact stationary state can still be unique because rare fluctuations generate transitions between these basins.

The stationary state may then be approximated by

[
\rho_{\mathrm{ss}}
\approx
\sum_{j=1}^{m}w_j\rho_j,
\qquad
w_j\geq0,
\qquad
\sum_jw_j=1,
]

where (\rho_j) are metastable branch states.

The weights satisfy a reduced stationary Markov equation

[
Q^\mathsf T w=0,
\qquad
\sum_jw_j=1,
]

where (Q_{ij}) contains switching rates between metastable branches.

Raw moments combine linearly:

[
\langle O\rangle_{\mathrm{ss}}
==============================

\sum_jw_j\langle O\rangle_j.
]

Cumulants should therefore be reconstructed only after the mixture has been assembled in raw moment coordinates.

CumulantHomotopy should distinguish two successful outcomes:

1. A single stationary branch.
2. A metastable branch mixture with uniquely determined weights.

This distinction is important because forcing several metastable branches to collapse into one local root may be physically incorrect.

## 4. Maximum entropy reconstruction

Given retained moments

[
m_i=\operatorname{Tr}(\rho O_i),
]

a maximum entropy state can be defined as

[
\rho_{\mathrm{ME}}(m)
=====================

\frac{
\exp
\left(
\log\rho_0+\sum_i\lambda_iO_i
\right)
}{
Z(\lambda)
},
]

where (\lambda_i) are chosen so that the state reproduces the retained moments.

The missing higher moments are then reconstructed as

[
\Phi_n^{\mathrm{ME}}(m)
=======================

\operatorname{Tr}
\left[
\rho_{\mathrm{ME}}(m)O_{>n}
\right].
]

Maximum entropy should be used in CumulantHomotopy in four ways.

### 4.1 Order lift predictor

When moving from order (n) to order (n+1), newly retained moments may be predicted using both the ordinary cumulant closure and the maximum entropy reconstruction:

[
z_{\mathrm{cum}}
================

\Phi_n^{\mathrm{cum}}(x),
]

[
z_{\mathrm{ME}}
===============

\Phi_n^{\mathrm{ME}}(x).
]

The disagreement

[
\Delta z
========

z_{\mathrm{cum}}-z_{\mathrm{ME}}
]

provides an estimate of closure uncertainty.

### 4.2 Realizable reconstruction

Maximum entropy gives a positive normalized state whenever the retained moments lie in the interior of the realizable moment set.

This allows a candidate cumulant root to be mapped to a physical state reconstruction and compared with the original closure.

### 4.3 Unique selector

If a finite stationary feasible set contains many admissible pseudomoment sequences, maximum entropy or minimum relative entropy can select one unique representative.

### 4.4 Closure homotopy

A homotopy between maximum entropy and cumulant closure can be defined by

[
\Phi_n^{(\lambda)}
==================

(1-\lambda)\Phi_n^{\mathrm{ME}}
+
\lambda\Phi_n^{\mathrm{cum}}.
]

This allows the method to begin from a realizable reconstruction and track toward the requested cumulant closure.

Maximum entropy does not by itself guarantee that the stationary reduced equations have one root. It should therefore be treated as a reconstruction, predictor, regularizer, and selector.

## 5. BBGKY style purification

Truncated BBGKY hierarchies can violate compatibility and positivity even when their equations are solved accurately.

A similar purification step can be used for quantum moment hierarchies.

Given a candidate moment vector (y), define the purified vector

[
\widetilde y
============

\arg\min_z
|W(z-y)|^2
]

subject to

[
M_q(z)\succeq0,
]

[
Cz=0,
]

and

[
A_{\mathrm{stat}}z=0.
]

Here:

[
M_q(z)\succeq0
]

enforces moment matrix positivity,

[
Cz=0
]

contains adjoint, symmetry, compatibility, and conservation constraints, and

[
A_{\mathrm{stat}}z=0
]

contains retained linear stationary identities.

The purification distance

[
d_{\mathrm{pur}}(y)
===================

|\widetilde y-y|
]

is a useful physicality diagnostic.

A small value indicates that the algebraic root lies close to a stationary realizable moment sequence.

A large value indicates that substantial correction is required and that the root is likely a closure artifact.

The first implementation should apply purification only at selected continuation endpoints. A later constrained continuation method may alternate between Newton correction and physical projection.

## 6. Mori Zwanzig and stationary resolvent closure

The exact raw moment hierarchy is linear.

Partition the moments into retained variables (x) and unresolved variables (y):

[
\begin{pmatrix}
\dot x\
\dot y
\end{pmatrix}
=============

\begin{pmatrix}
A_{xx}&A_{xy}\
A_{yx}&A_{yy}
\end{pmatrix}
\begin{pmatrix}
x\
y
\end{pmatrix}.
]

At stationarity,

[
A_{yx}x+A_{yy}y=0.
]

If the unresolved block is invertible on the relevant subspace,

[
y
=

-A_{yy}^{-1}A_{yx}x.
]

Substitution gives the exact reduced stationary equation

[
\left(
A_{xx}
------

A_{xy}A_{yy}^{-1}A_{yx}
\right)x=0.
]

This is the zero frequency Mori Zwanzig memory correction or stationary Schur complement.

The reduced problem remains linear and therefore avoids the nonlinear algebraic root explosion produced by ordinary cumulant closure.

The unresolved inverse can be approximated using:

1. Krylov subspaces.
2. Block Lanczos methods.
3. Continued fraction expansions.
4. Rational approximations.
5. Auxiliary memory variables.
6. Finite unresolved hierarchy sections.
7. Reset regularization.

A reset regularized inverse is

[
(A_{yy}-\epsilon I)^{-1},
]

followed by continuation toward

[
\epsilon\rightarrow0.
]

CumulantHomotopy should include a future stationary backend such as

```julia
ResolventClosure()
```

which provides a linear or rational stationary approximation.

This approximation can be used as:

1. A reference solution.
2. An order lift predictor.
3. A comparison against nonlinear cumulant branches.
4. A possible uniqueness preserving alternative to hard cumulant factorization.

## 7. Koopman invariance and adaptive basis construction

An exact finite hierarchy exists when the retained operator space is invariant under the adjoint Liouvillian.

For an operator basis

[
\mathcal B_n
============

{O_1,\ldots,O_m},
]

exact closure requires

[
\mathcal L^\dagger
\operatorname{span}(\mathcal B_n)
\subseteq
\operatorname{span}(\mathcal B_n).
]

Then

[
\mathcal L^\dagger O_i
======================

\sum_jK_{ij}O_j,
]

and the retained expectations obey an exact finite linear system.

Generic interacting systems do not possess small exact invariant spaces, but approximate invariance can still guide basis selection.

Define an invariance defect

[
\eta(\mathcal B_n)
==================

\left|
(I-P_{\mathcal B_n})
\mathcal L^\dagger
P_{\mathcal B_n}
\right|.
]

QuantumCumulants already records which new moments are generated when applying (\mathcal L^\dagger). This information can be used to add the operators responsible for the largest invariance defect.

The basis should therefore be selectable by more than polynomial order.

Possible adaptive criteria include:

1. Largest contribution to the invariance defect.
2. Largest predicted stationary amplitude.
3. Largest effect on requested observables.
4. Strongest coupling in the moment graph.
5. Largest contribution to the order lift defect.
6. Symmetry sector relevance.

An adaptive approximately invariant basis may use fewer variables than a complete order based hierarchy and may reduce the number of artificial stationary roots.

## 8. Lumpability and metastable coarse graining

Lumpability is the Markov chain analogue of an exact closed observable space.

A coarse partition is lumpable when the coarse probabilities evolve autonomously as a Markov process.

For CumulantHomotopy, lumpability provides two useful concepts.

First, it gives a criterion for identifying exactly or approximately closed groups of observables.

Second, it provides the mathematical model for combining metastable cumulant branches.

If several roots represent metastable basins, the branch set can be treated as a coarse Markov state space with switching generator (Q).

The stationary branch weights then follow from

[
Q^\mathsf T w=0.
]

Approximate lumpability diagnostics can determine whether the metastable description is self consistent.

## 9. Finite state projection and Fock space truncation

Finite state projection truncates the original state space while preserving a finite linear generator.

For bosonic quantum systems, the analogue is a finite Fock space cutoff.

The resulting stationary problem is

[
\mathcal L_N\rho_N=0,
\qquad
\operatorname{Tr}\rho_N=1.
]

When the truncated generator remains primitive, it has one stationary density matrix.

Finite state projection and Fock space truncation should be used as:

1. Small system reference solvers.
2. Sources of initial stationary moments.
3. Sources of occupation and energy tail estimates.
4. Benchmarks for hierarchy convergence.
5. Benchmarks for metastable switching.
6. Hybrid solvers for strongly quantum subsystems.

A hybrid decomposition may use an explicit Hilbert space for a small strongly quantum subsystem and a cumulant description for the remaining modes or collective degrees of freedom.

Finite state projection does not scale to all systems, but it preserves the original linear stationary structure and is therefore an important reference method.

## 10. Phase space methods

Phase space formulations are particularly valuable for the quantum to classical limit and for identifying metastable structure.

The exact Wigner equation can be written schematically as

[
\partial_tW
===========

\mathcal L_{\mathrm{cl}}W
+
\hbar^2\mathcal L_2W
+
\hbar^4\mathcal L_4W
+\cdots.
]

The truncated Wigner approximation discards derivative terms above second order, producing a Fokker Planck equation or an equivalent stochastic differential equation.

For nonlinear bosonic systems, phase space methods should be incorporated into CumulantHomotopy as predictors, reference models, and metastability estimators.

### 10.1 Semiclassical stationary predictor

A truncated Wigner stationary distribution can provide approximate moments

[
\langle(a^\dagger)^pa^q\rangle_{\mathrm{TW}}.
]

These moments can initialize the cumulant stationary problem.

This is superior to initializing from a deterministic mean field fixed point when noise connects several semiclassical basins.

### 10.2 Metastable branch weights

Truncated Wigner, positive (P), gauge (P), or quantum trajectory simulations can estimate switching rates between metastable branches.

These rates define the reduced generator (Q), from which the unique stationary branch mixture can be constructed.

### 10.3 Phase space order lift

Samples from a phase space distribution can estimate newly required moments:

[
z_{\mathrm{PS}}
===============

\Phi_n^{\mathrm{PS}}(x).
]

The method can compare three order lift predictors:

[
z_{\mathrm{cum}},
\qquad
z_{\mathrm{ME}},
\qquad
z_{\mathrm{PS}}.
]

Agreement among the predictors provides evidence that the hierarchy is locally well controlled.

Large disagreement identifies a difficult nonGaussian or strongly quantum regime.

### 10.4 Variable scaling

Phase space fluctuations provide natural scales for high order moments and cumulants.

For a large occupation scale (N), standardized cumulants may be defined by

[
\widetilde\kappa_r
==================

\frac{\kappa_r}{N^{r/2}}.
]

Such scaling can improve conditioning and provide an order independent interpretation of cumulant magnitude.

### 10.5 Quantum correction homotopy

A quantum to classical homotopy can be defined through the Wigner generator:

[
\mathcal L_W(\lambda)
=====================

\mathcal L_{\mathrm{TW}}
+
\lambda\mathcal L_{\mathrm{quantum}}.
]

At

[
\lambda=0,
]

the dynamics is the truncated Wigner approximation.

At

[
\lambda=1,
]

the complete Wigner generator is recovered.

A finite spectral or moment representation of this generator could be used to continue stationary solutions from the semiclassical problem toward the quantum problem.

### 10.6 Limitations

Truncated Wigner is not a general stationary state certificate.

The neglected quantum terms may accumulate over long times, and the approximation can approach an incorrect classical stationary distribution.

A Wigner function may also be negative and therefore cannot serve as a positivity certificate.

Positive (P) can represent many bosonic systems exactly, but it may suffer from large sampling variance, boundary terms, and trajectory instabilities.

Phase space methods should therefore guide and validate the cumulant method rather than replace moment realizability tests.

## 11. Proposed architecture

CumulantHomotopy should be organized into four layers.

### 11.1 Basis construction

The basis layer chooses retained observables using:

1. Cumulant order.
2. Symmetry.
3. Adjoint relations.
4. Moment graph connectivity.
5. Koopman invariance defect.
6. User requested observables.
7. Occupation and energy bounds.

A future interface could be:

```julia
basis = AdaptiveMomentBasis(
    max_order = 8,
    invariance_tolerance = 1e-4,
    symmetries = :auto,
)
```

### 11.2 Candidate generation

Candidate stationary branches may be generated using:

1. Physical parameter continuation.
2. Closure consistent order continuation.
3. Maximum entropy prediction.
4. Phase space prediction.
5. Resolvent closure.
6. Finite state projection references.
7. Branch spawning near singularities.

```julia
candidates = stationary_candidates(
    hierarchy;
    branch_budget = 4,
    predictors = (
        CumulantPredictor(),
        MaximumEntropyPredictor(),
        PhaseSpacePredictor(),
        ResolventPredictor(),
    ),
)
```

### 11.3 Physical correction and validation

Each candidate is evaluated using:

1. Stationary residual.
2. Adjoint consistency.
3. Held out Liouvillian equations.
4. Moment matrix positivity.
5. BBGKY style purification distance.
6. Order lift defect.
7. Shared moment convergence.
8. Path dependence.
9. Optional local root certification.
10. Optional semidefinite observable bounds.

```julia
validated = validate_candidates(
    candidates;
    purification = MinimalMomentPurification(),
    certification = KrawczykCertification(),
)
```

### 11.4 Stationary interpretation

The result should distinguish between:

```julia
SingleStationaryBranch(...)
```

and

```julia
MetastableMixture(
    branches,
    weights,
    switching_generator,
)
```

The returned status should indicate whether the physical result is supported by one converged branch or by a weighted mixture of metastable branches.

## 12. Recommended first combined prototype

The first prototype should combine:

[
\text{cumulant continuation}
+
\text{phase space metastable analysis}
+
\text{moment realizability validation}.
]

A driven Kerr resonator provides a suitable initial benchmark.

The experiment should:

1. Generate all stationary roots at low hierarchy order.
2. Track a bounded branch set through parameter space.
3. Continue each candidate through increasing hierarchy order.
4. Identify dim and bright metastable branches where present.
5. Estimate switching rates using truncated Wigner, positive (P), or trajectories.
6. Construct the weighted stationary raw moments.
7. Compare against the exact or high accuracy stationary state.
8. Apply moment matrix positivity tests.
9. Evaluate held out stationary equations.
10. Compare cumulant, maximum entropy, phase space, and resolvent order lift predictors.
11. Study convergence as occupation increases toward the semiclassical limit.
12. Record whether the result is a single branch or a metastable mixture.

The purpose of the benchmark is to determine whether physical roots literally merge with increasing order or instead organize into metastable components whose unique weighted mixture converges to the exact stationary state.
