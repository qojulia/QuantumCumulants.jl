# Technical specification

## 1. Purpose

CumulantHomotopy.jl is a hierarchy aware framework for computing, comparing, and validating stationary solutions of truncated quantum cumulant and moment equations.

The project addresses two coupled problems:

1. How to generate stationary candidates of a nonlinear truncated hierarchy without relying on uncontrolled random initial conditions.
2. How to determine whether the physically admissible predictions of successive finite hierarchies collapse toward the unique stationary state of the original Liouvillian.

The project is not defined by one numerical backend. Homotopy continuation, nonlinear solves, time integration, pseudoarclength continuation, and semidefinite optimization are complementary tools inside a larger method.

## 2. Exact problem and source of artificial roots

### 2.1 Infinite linear hierarchy

Let `\mathcal A` be an operator algebra and define

```math
y(O)=\operatorname{Tr}(\rho O).
```

The adjoint Lindblad generator produces the exact moment hierarchy

```math
\frac{d}{dt}y(O)=y(\mathcal L^\dagger O).
```

After choosing an ordered operator basis, this can be written formally as

```math
\dot y = A y,
```

where `A` is generally infinite dimensional. The physical state space is not the full vector space of sequences. A physical sequence must define a normalized positive linear functional:

```math
y(1)=1,
\qquad
y(X^\dagger X)\geq0
```

for every admissible operator polynomial `X`.

If the Liouvillian is primitive or otherwise known to possess one stationary state, then

```math
\ker A
\cap
\{\text{normalized positive functionals}\}
=
\{y_{\mathrm{ss}}\}.
```

The exact uniqueness statement is therefore a statement about the intersection of a linear kernel with a positive cone.

### 2.2 Nonlinear finite closure

At truncation order `n`, let `P_n` project onto retained moments and let

```math
S_n:x\mapsto y
```

reconstruct unresolved moments using a cumulant closure. The finite equations are

```math
\dot x=f_n(x)=P_n A S_n(x).
```

The stationary polynomial system is

```math
F_n(x)=f_n(x)=0.
```

Additional roots are introduced by the nonlinear reconstruction manifold `S_n(x)` and by projection. They are not stationary states of the original linear Liouvillian merely because they solve `F_n(x)=0`.

### 2.3 Root sets and physical sets

Define the algebraic root set

```math
\mathcal R_n=\{x:F_n(x)=0\}.
```

The total number of isolated complex roots may be controlled only by algebraic bounds such as total degree or mixed volume bounds. These counts do not encode Hermiticity, positivity, moment extendibility, stability, or consistency with neighbouring orders.

Define a hierarchy of filtered candidate sets

```math
\mathcal P_n^{(q,s,\varepsilon)}
\subseteq
\mathcal R_n,
```

where the filters may include:

1. Reality or adjoint consistency.
2. Symmetry constraints.
3. Moment matrix positivity through level `q`.
4. Positive extension through `s` additional moment levels.
5. Stationary and held out residual tolerances `\varepsilon`.
6. Energy or number growth bounds.
7. Continuity from a physical reference state.

The project studies how these filtered sets evolve with hierarchy order.

## 3. Definition of hierarchy collapse

The desired phenomenon is not necessarily that all algebraic roots merge.

For a fixed set of low order observables represented by a projection `\pi_k`, define the candidate diameter

```math
D_{n,k}
=
\operatorname{diam}\!\left(\pi_k\mathcal P_n\right).
```

The hierarchy exhibits physical collapse on these observables when

```math
D_{n,k}\longrightarrow0.
```

A second notion uses a convex stationary moment envelope. Let

```math
K_{n,q}
=
\left\{
 y:
 y(1)=1,
 \;y(\mathcal L^\dagger O)=0
 \text{ for retained }O,
 \;M_q(y)\succeq0,
 \;\text{algebraic and growth constraints hold}
\right\}.
```

For an observable `O`, define

```math
\underline O_{n,q}=\inf_{y\in K_{n,q}}y(O),
\qquad
\overline O_{n,q}=\sup_{y\in K_{n,q}}y(O),
```

and interval width

```math
W_{n,q}(O)=\overline O_{n,q}-\underline O_{n,q}.
```

Collapse is supported when both

```math
D_{n,k}\to0
```

and

```math
W_{n,q}(O)\to0
```

for the observables of interest.

The candidate roots and the convex envelope answer different questions. Roots describe the selected nonlinear closure. The envelope describes moment sequences compatible with finite stationary and positivity information from the original Liouvillian.

## 4. Method architecture

CumulantHomotopy consists of three interacting layers.

### 4.1 Layer A: candidate generation

Generate one or more roots of `F_n(x)=0` using physically informed methods.

Supported strategies should include:

1. Newton or nonlinear least squares from user supplied or automatically constructed initial conditions.
2. Long time integration followed by root polishing.
3. Physical parameter continuation.
4. Reset regularized continuation.
5. Continuation between hierarchy orders.
6. Pseudoarclength continuation near real folds.
7. Optional polynomial homotopy continuation.
8. Low order root enumeration for validation experiments.

No candidate generator by itself proves physicality.

### 4.2 Layer B: physical validation and collapse analysis

For every candidate, compute stationarity, representability, cross order, and path consistency diagnostics. Construct stationary moment envelopes and observable bounds where practical.

This layer produces the evidence needed to distinguish:

1. A root of the closure.
2. A dynamically selected fixed point of the closure.
3. A physically plausible approximation to the Liouvillian stationary state.
4. A cluster of branches that are indistinguishable on selected low order observables.

### 4.3 Layer C: closure design

Use defects discovered by Layer B to improve the reconstruction `S_n`.

Research directions include:

1. Maximum entropy and exponential family reconstruction.
2. Variational moment closure.
3. Invariant manifold corrections.
4. Representability purification.
5. Adaptive operator basis selection.
6. Memory augmented closures.
7. Convex projection of a cumulant prediction onto a stationary moment envelope.

Layer C is a research track. It should not block implementation of Layers A and B.

## 5. Stationary hierarchy representation

### 5.1 Core abstraction

```julia
struct StationaryHierarchy
    equations
    variables
    parameters
    moments
    operator_orders
    adjoint_pairs
    symmetry_constraints
    closure_metadata
    residual_evaluator
    jacobian_evaluator
end
```

The representation must preserve a bidirectional map between numerical coordinates and symbolic QuantumCumulants expressions.

### 5.2 Coordinate policies

Realification is optional, not mandatory.

#### Complexified coordinates

Treat `m_O` and `m_{O^\dagger}` as independent complex variables. This is natural for polynomial homotopy backends. Physical adjoint consistency is tested at the endpoint:

```math
r_\dagger
=
\max_O|m_{O^\dagger}-m_O^*|.
```

#### Reduced real coordinates

For an adjoint pair, use

```math
m_O=q_O+ip_O,
\qquad
m_{O^\dagger}=q_O-ip_O.
```

Hermitian moments receive one real variable. This representation enforces Hermiticity exactly and is natural for native physical solvers.

#### Policy interface

```julia
abstract type CoordinatePolicy end
struct ComplexifiedCoordinates <: CoordinatePolicy end
struct ReducedRealCoordinates <: CoordinatePolicy end
```

Backends select a compatible policy without exposing the representation to high level users.

### 5.3 Deterministic reduction

The exporter must remove redundant equations and variables deterministically using:

1. Adjoint pairs.
2. Hermitian observables.
3. Operator identities.
4. Symmetry sectors.
5. Conserved quantities where suitable.

Random squaring of an overdetermined system should be a fallback for algebraic geometry experiments, not the default package representation.

## 6. Candidate generation

### 6.1 Baseline stationary solves

Before continuation, the package should support baseline methods:

1. `NonlinearSolve.jl` root solves.
2. Nonlinear least squares when the system is overdetermined.
3. ODE integration to a stationary point.
4. Multiple physically motivated starts.

These baselines are required to measure the actual benefit of continuation.

### 6.2 Physical parameter continuation

For a parameter path `p(t)`, solve

```math
H(x,t)=F_n(x;p(t))=0.
```

The tangent satisfies

```math
J_x\dot x=-\partial_tH.
```

A predictor corrector implementation should support adaptive step size, equation scaling, variable scaling, sparse Jacobians, and conditioning estimates.

The physical path should be chosen so that the start state is known or easily computed. Examples include:

1. Zero interaction.
2. Zero coherent coupling.
3. Independent damping or thermalization.
4. A Gaussian reference model.
5. A reset regularized Liouvillian.

### 6.3 Reset regularized path

Choose a full rank reference state `\sigma` and define

```math
\mathcal R_\sigma(\rho)
=
\operatorname{Tr}(\rho)\sigma-\rho.
```

Use

```math
\mathcal L_{\lambda,\epsilon}
=
\lambda\mathcal L_{\mathrm{target}}
+
\epsilon\mathcal R_\sigma.
```

A default path is

```math
(\lambda,\epsilon):
(0,\epsilon_0)
\rightarrow
(1,\epsilon_0)
\rightarrow
(1,0).
```

At `\lambda=0`, the exact stationary state is `\sigma`. For `\epsilon>0`, the reset term regularizes the exact generator on the traceless subspace. This does not prevent the truncated closure from developing artificial branches, so path ambiguity must still be tested.

### 6.4 Continuation between hierarchy orders

Let `x` contain variables retained at order `n` and `z` variables first retained at order `n+1`. The lower order cumulant closure defines

```math
z=\Phi_n(x).
```

Introduce closure deviation coordinates

```math
u=z-\Phi_n(x).
```

At `u=0`, all newly retained moments equal their lower order closure values.

A padded start system is

```math
G_0(x,u)
=
\begin{pmatrix}
F_n(x)\\
D u
\end{pmatrix},
```

where `D` is a scaling or damping matrix. The target is

```math
G_1(x,u)
=
F_{n+1}\!\left(x,\Phi_n(x)+u\right).
```

Possible order transitions include:

1. Direct Newton correction from `(x_n,0)`.
2. A linear start target homotopy.
3. Separate activation of corrections to old equations and equations for new moments.
4. A pseudoarclength path when the transition has folds.
5. A least squares predictor for `u` before continuation.

### 6.5 Linearized new moment predictor

At the closure consistent lift, let

```math
r=F_{n+1}\!\left(x_n,\Phi_n(x_n)\right).
```

With `J_u` the Jacobian with respect to new closure deviations, compute

```math
\Delta u
=
\arg\min_v
\|W(r+J_uv)\|_2^2
+
\gamma\|Sv\|_2^2.
```

This predictor uses the higher order equations to estimate how the newly retained moments should depart from the lower order factorization.

### 6.6 HomotopyContinuation.jl backend

HomotopyContinuation.jl is useful for:

1. Reference path tracking at low and moderate dimensions.
2. Parameter homotopies with fixed polynomial support.
3. Low order root enumeration.
4. Comparison with a native tracker.
5. Singular and difficult polynomial test cases.

It is not required as the first project milestone and should not dictate the public API. The package should first establish the hierarchy representation, closure lift, diagnostics, and benchmark harness. A homotopy backend can then be added as one candidate generator.

### 6.7 Pseudoarclength and branch continuation

A real branch may pass through a fold where direct parameter continuation fails. A pseudoarclength method augments the equations with

```math
\tau_x^T(x-x_0)+\tau_p(p-p_0)-\Delta s=0.
```

BifurcationKit.jl is an appropriate reference backend for this problem. A native implementation can be considered after profiling.

## 7. Physical diagnostics

### 7.1 Stationary residual

```math
r_n=\|W_nF_n(x_n)\|.
```

All accepted candidates require a small residual, but a small residual alone has no physical meaning beyond solving the closure.

### 7.2 Adjoint and symmetry residuals

When not enforced by coordinates, report

```math
r_\dagger
=
\max_O|m_{O^\dagger}-m_O^*|.
```

Evaluate all known symmetry and conservation constraints separately.

### 7.3 Held out stationary equations

Select operators `O` not used to construct the closed stationary system and evaluate

```math
h_O=y(\mathcal L^\dagger O).
```

Held out equations are out of sample tests of the closure. Operator selection may use the existing moment graph or repeated application of `\mathcal L^\dagger` to observables of interest.

### 7.4 Order lift defect

The closure consistent next order defect is

```math
\delta_n
=
\left\|
W_{n+1}F_{n+1}\!\left(x_n,\Phi_n(x_n)\right)
\right\|.
```

Residuals should be separated into:

1. Equations inherited from order `n`.
2. Equations first introduced at order `n+1`.

### 7.5 Shared moment convergence

For a fixed projection order `k`,

```math
e_{n,k}
=
\left\|
S_k\left(\pi_kx_n-\pi_kx_{n-1}\right)
\right\|.
```

The package should inspect several consecutive orders and report odd and even subsequences separately when parity effects occur.

### 7.6 Stability metadata

For the truncated dynamical Jacobian,

```math
\alpha_n
=
\max\operatorname{Re}\operatorname{eig}J_n(x_n).
```

Stability under the closed dynamics is useful metadata. It is not a necessary or sufficient criterion for physicality.

### 7.7 Path dependence

Compare endpoints obtained from:

1. Different physical parameter paths.
2. Different reset states.
3. Forward and reverse continuation.
4. Closed loops in parameter space.
5. Different candidate generators.

Path dependence indicates branch ambiguity in the truncated problem and must be reported rather than hidden.

## 8. Moment realizability and stationary envelopes

### 8.1 Moment matrices

Given an operator basis `B=(B_1,\ldots,B_m)`, define

```math
[M_B(y)]_{ij}=y(B_i^\dagger B_j).
```

A necessary physical condition is

```math
M_B(y)\succeq0.
```

Report:

```math
\mu_B=\lambda_{\min}(M_B)
```

and the total negative eigenvalue weight. Passing a finite family of tests is evidence, not proof of a full density operator.

### 8.2 Positive extension tests

Given candidate moments through order `k`, introduce unknown higher moments and test whether they can satisfy:

```math
A_ny=0,
\qquad
y(1)=1,
\qquad
M_q(y)\succeq0,
```

plus operator identities and growth bounds.

This is a semidefinite feasibility problem. Increasing the extension depth provides a stronger filter.

### 8.3 Stationary moment envelope

The set `K_{n,q}` is convex because its stationarity and algebraic constraints are linear and its positivity constraints are semidefinite.

Operations should include:

1. Lower and upper bounds for selected observables.
2. Feasibility of a candidate low order moment vector.
3. Distance from a candidate to the feasible set.
4. Projection of a candidate onto the feasible set.
5. A unique central estimate selected by a strictly convex objective.

### 8.4 Projection of a cumulant candidate

Given raw moments `y_c` obtained from a cumulant root, solve

```math
\begin{aligned}
\min_y\quad&
\frac12\|W(y-y_c)\|_2^2,\\
\text{subject to}\quad&
y\in K_{n,q}.
\end{aligned}
```

This gives a unique projection when the feasible set is nonempty and closed and the objective is strictly convex.

The projection has three uses:

1. Quantify the physical defect of a cumulant root.
2. Produce a realizability corrected stationary estimate.
3. Determine whether distinct cumulant branches collapse to the same physical pseudomoment vector.

### 8.5 Unique finite order selector

More generally, define

```math
y_{n,q}^\star
=
\arg\min_{y\in K_{n,q}}R_n(y),
```

where `R_n` is strictly convex. Choices include:

1. Distance to the previous order estimate.
2. Distance to a cumulant prediction.
3. A regularized log determinant center.
4. A maximum entropy or relative entropy objective where tractable.

The selector provides one reproducible finite order result. The interval widths of `K_{n,q}` remain the measure of unresolved physical ambiguity.

### 8.6 Bosonic growth control

For infinite dimensional bosonic systems, moment positivity alone is insufficient for useful compactness. The hierarchy should support constraints such as

```math
\langle N^r\rangle\leq C_r
```

or Lyapunov inequalities derived from `\mathcal L^\dagger V`.

When no occupation or energy control is available, the package must report that moment tightness is unverified.

## 9. Collapse metrics and branch clustering

### 9.1 Algebraic count

Record the number of algebraic roots only when enumeration is feasible. This is a diagnostic of the closure geometry, not the primary convergence measure.

### 9.2 Physical survivor count

For specified validation levels, define

```math
N_{n,q,s}
=
\#\mathcal P_n^{(q,s,\varepsilon)}.
```

For fixed closure order `n`, increasing the extension depth or validation strength should not increase the number of survivors when filters are nested.

No monotonicity is assumed across closure order because the systems `F_n` are not nested.

### 9.3 Projected branch diameter

For selected observables,

```math
D_{n,k}
=
\max_{x,x'\in\mathcal P_n}
\|\pi_kx-\pi_kx'\|.
```

This is the principal root based collapse metric.

### 9.4 Envelope width

For each observable,

```math
W_{n,q}(O)
=
\overline O_{n,q}-\underline O_{n,q}.
```

This is the principal convex collapse metric.

### 9.5 Distance to physical envelope

For a candidate `x`, convert to raw moments and define

```math
d_{n,q}(x)
=
\operatorname{dist}(y(x),K_{n,q}).
```

A large distance is evidence that the closure root is incompatible with retained stationary positivity information.

### 9.6 Collapse status

Suggested statuses:

```julia
:algebraic_only
:physically_plausible
:multiple_physical_clusters
:collapsed_on_requested_observables
:envelope_not_collapsed
:positivity_violation
:no_positive_extension
:path_dependent
:insufficient_growth_control
:tracking_failed
```

## 10. Structure preserving closure research

### 10.1 Maximum entropy and exponential families

Use a positive state manifold

```math
\rho_\lambda
=
\frac{
\exp\!\left(\log\rho_0+\sum_j\lambda_jO_j\right)
}{Z(\lambda)}.
```

This enforces positivity and normalization. The projected dynamics may be obtained using an information geometric metric or variational principle.

Advantages:

1. Realizability by construction.
2. A smooth finite dimensional manifold.
3. Natural entropy and covariance geometry.

Limitations:

1. Evaluating the log partition function can be expensive.
2. The reduced dynamics can still have several fixed points.
3. Generic nonequilibrium Lindbladians need not admit a strict convex potential on the manifold.

### 10.2 Variational moment closure

Interpret the closure as projection of the full dynamics onto a parametric family. Minimize an instantaneous residual or divergence rather than imposing zero high cumulants without justification.

The variational formulation should be benchmarked against conventional cumulant factorization using the same observables and dimensions.

### 10.3 Invariant manifold correction

For a state manifold `M_n`, define the invariance defect

```math
\Delta(\rho)
=
\left(I-P_{T_\rho M_n}\right)\mathcal L(\rho).
```

Correct the manifold or closure map by reducing `\|\Delta\|`. Newton type invariant manifold methods provide a model for this procedure.

A hierarchy specific version may use the moment graph to identify which omitted moments dominate the defect.

### 10.4 Thermodynamic projection

For detailed balance Lindbladians with a relative entropy gradient structure, construct a projection that preserves entropy production. This may give a finite dynamics with stronger uniqueness and stability properties than an arbitrary moment projection.

This approach is restricted to systems with suitable thermodynamic structure and should be treated separately from generic driven open systems.

### 10.5 Representability correction

Inspired by BBGKY purification, minimally modify either:

1. The candidate moments.
2. The closure reconstruction.
3. The projected equations.

so that selected moment matrices remain positive while disturbing the original closure as little as possible.

Corrections must be monitored for conservation, symmetry, and stationarity defects.

### 10.6 Adaptive operator spaces

Let `V_n` be the retained observable space. Measure its invariance defect under `\mathcal L^\dagger`:

```math
E_n=(I-P_{V_n})\mathcal L^\dagger P_{V_n}.
```

Add operators responsible for the largest components of `E_n`, rather than increasing polynomial order uniformly. This connects the project to invariant subspace, Koopman, and adaptive moment methods.

### 10.7 Memory augmented reduction

Exact elimination of unresolved variables generally produces memory. A finite autonomous closure may therefore need auxiliary variables approximating a memory kernel.

A stationary solver may not require the full transient memory model, but memory based analysis can explain why a Markovian cumulant closure creates artificial fixed points. Auxiliary variables should be investigated when stationary errors correlate with slow unresolved modes.

## 11. Backend independent API

### 11.1 Candidate methods

```julia
abstract type AbstractCandidateMethod end

struct NewtonCandidates <: AbstractCandidateMethod end
struct TimeIntegrationCandidates <: AbstractCandidateMethod end
struct ParameterContinuation <: AbstractCandidateMethod end
struct OrderContinuation <: AbstractCandidateMethod end
struct PolynomialHomotopy <: AbstractCandidateMethod end
struct PseudoArclengthContinuation <: AbstractCandidateMethod end
```

### 11.2 Validation methods

```julia
abstract type AbstractValidator end

struct AdjointValidator <: AbstractValidator end
struct HeldoutStationarity <: AbstractValidator end
struct MomentMatrixValidator <: AbstractValidator end
struct PositiveExtensionValidator <: AbstractValidator end
struct StationaryEnvelopeValidator <: AbstractValidator end
```

### 11.3 High level entry points

```julia
stationary_candidates(hierarchy; methods, kwargs...)
stationary_state(hierarchy; methods, validation, selector, kwargs...)
stationary_sequence(hierarchies; methods, validation, collapse, kwargs...)
continue_parameters(hierarchy, path, solution; kwargs...)
continue_order(lower, higher, solution; kwargs...)
stationary_envelope(hierarchy; level, observables, kwargs...)
validate_stationary(candidate, hierarchy; validators, kwargs...)
```

### 11.4 Result types

```julia
struct StationaryCandidate
    values
    moments
    parameters
    order
    source_method
    branch_id
    residuals
    diagnostics
    status
end

struct StationarySequence
    orders
    candidates
    selected
    envelope_results
    branch_clusters
    collapse_metrics
    status
end
```

## 12. Backend strategy

### 12.1 Initial dependencies

The first implementation should use the existing SciML stack where possible:

1. ModelingToolkit generated residuals and Jacobians.
2. NonlinearSolve.jl for baseline roots.
3. OrdinaryDiffEq.jl for integration based candidates.
4. LinearSolve.jl and SparseArrays for correctors.
5. An optional semidefinite backend selected later.

### 12.2 Optional continuation backends

1. HomotopyContinuation.jl for polynomial tracking and low order enumeration.
2. BifurcationKit.jl for pseudoarclength continuation and stability.
3. A native hierarchy aware predictor corrector after profiling.

### 12.3 Custom tracker criterion

A custom low level continuation engine is justified when profiling shows that a general backend cannot exploit:

1. Sparse hierarchy Jacobians.
2. Matrix free Jacobian vector products.
3. Block structure by moment order.
4. Factorization reuse between nearby path points.
5. Physical branch beam management.
6. Mixed precision only on difficult segments.

The hierarchy formulation and validation layer should be developed before replacing mature continuation machinery.

## 13. Verification strategy

### 13.1 Exact finite dimensional references

For finite spin or truncated bosonic systems, construct the full Liouvillian and solve

```math
\mathcal L\rho_{\mathrm{ss}}=0,
\qquad
\operatorname{Tr}\rho_{\mathrm{ss}}=1.
```

These models provide exact stationary moments and a known terminal operator basis.

### 13.2 Benchmark models

1. Driven Kerr resonator with exact or high accuracy stationary moments.
2. Parametric oscillator and dissipative cat model.
3. Finite spin models with exact Liouvillian solutions.
4. Permutation symmetric light matter model.
5. A model with metastability and a small Liouvillian gap.
6. A classical chemical master equation with known closure pathologies.
7. A BBGKY style finite bosonic model with representability issues.

### 13.3 Comparative methods

Compare:

1. Random start Newton.
2. Physically initialized Newton.
3. Long time integration.
4. Parameter continuation.
5. Order continuation.
6. Polynomial homotopy continuation.
7. Pseudoarclength continuation.
8. Stationary moment envelope selectors.
9. Exact finite dimensional stationary states.

### 13.4 Metrics

1. Error in selected observables.
2. Runtime and memory.
3. Candidate and physical survivor counts.
4. Branch diameter on low order observables.
5. Stationary envelope widths.
6. Distance to the stationary envelope.
7. Order lift defect.
8. Held out residuals.
9. Positivity and extension margins.
10. Path dependence.
11. Failure classification quality.
12. Scaling with hierarchy order.

## 14. Research hypotheses

### Hypothesis A: physical collapse

For a primitive open quantum system and a suitable sequence of physical filters, the projected physically admissible candidate sets contract on fixed low order observables over a useful parameter regime.

### Hypothesis B: closure consistent prediction

The lift `z=\Phi_n(x_n)` and its linearized correction provide substantially better initial conditions for order `n+1` than zero initialization or unrelated random starts.

### Hypothesis C: stationary envelope discrimination

Moment positivity, positive extension, and held out stationary constraints reject a significant fraction of algebraically valid but physically spurious cumulant roots.

### Hypothesis D: branch clustering

Several formal roots may survive while their low order observables and projections onto the stationary envelope converge to one cluster.

### Hypothesis E: adaptive bases

Selecting new observables by the invariance defect under `\mathcal L^\dagger` reaches a given stationary accuracy with fewer variables than uniform polynomial order truncation.

### Hypothesis F: corrected closures

A realizability or invariance corrected reconstruction improves stationary convergence without requiring full state space truncation.

These are research hypotheses and must not be presented as general theorems without additional assumptions.

## 15. Mathematical targets

### 15.1 Finite dimensional target

For a finite operator algebra, a complete retained basis reduces the stationary problem to the original finite linear Liouvillian system. The method should recover exact uniqueness at the terminal level.

### 15.2 Nested envelope target

Let `K_n` be a nested sequence of stationary moment relaxations. Under compactness, completeness, and uniqueness assumptions, aim to establish

```math
\bigcap_{n\geq k}\pi_kK_n
=
\{\pi_ky_{\mathrm{ss}}\}.
```

Then

```math
\operatorname{diam}(\pi_kK_n)\to0,
```

and every continuous unique selector on `K_n` converges on fixed low order moments.

### 15.3 Root based target

A weaker empirical target is

```math
\operatorname{diam}(\pi_k\mathcal P_n)\to0,
```

where `\mathcal P_n` is defined by explicit validation levels. This target does not assume monotonicity of the raw algebraic root count.

### 15.4 Infinite dimensional requirements

For bosonic systems, convergence statements require additional control such as:

1. Existence of the relevant moments.
2. Tightness or energy bounds.
3. Moment determinacy or an appropriate operator algebra topology.
4. Completeness of the selected observable sequence.

The package must make these assumptions visible in result metadata.

## 16. Risks and mitigations

### Risk: continuation selects a path dependent root

Mitigation: multiple candidate generators, alternate paths, branch beams, loop tests, and explicit path dependence status.

### Risk: positivity tests pass at low order but no state exists

Mitigation: increasing moment matrix levels, positive extension tests, growth constraints, and conservative wording.

### Risk: the stationary envelope is too loose

Mitigation: adaptive operator selection, additional held out stationarity constraints, symmetry reduction, and higher moment levels.

### Risk: semidefinite programs dominate runtime

Mitigation: run envelope and extension tests at selected endpoints, use small adaptive operator bases, and keep them optional.

### Risk: higher cumulant order is nonmonotonic

Mitigation: report several consecutive orders, parity subsequences, branch clusters, and envelope widths rather than one difference.

### Risk: bosonic moments escape to infinity

Mitigation: derive or request energy and occupation bounds and flag unverified tightness.

### Risk: a finite autonomous closure cannot reproduce slow memory

Mitigation: diagnose slow unresolved modes, investigate auxiliary memory variables, and separate stationary accuracy from transient accuracy.

### Risk: implementation starts with an unnecessarily restrictive backend

Mitigation: backend independent hierarchy and validation interfaces, baseline nonlinear methods first, optional homotopy and pseudoarclength backends later.

## 17. Literature map

The following areas directly inform the method.

### 17.1 Stationary finite state projection

Gupta, Mikelson, and Khammash develop a finite state projection algorithm for stationary chemical master equations, with convergence and error control under stability conditions. This is a model for preserving a finite linear stationary problem while expanding the approximation space.

* A. Gupta, J. Mikelson, M. Khammash, *A finite state projection algorithm for the stationary solution of the chemical master equation*, 2017. <https://arxiv.org/abs/1704.07259>

### 17.2 Rigorous stationary moment bounds

Kuntz, Thomas, Stan, and Barahona use semidefinite and linear programming to compute convergent stationary moment and distribution bounds for chemical master equations. In the unique case, the bounds provide a computational uniqueness test.

* J. Kuntz, P. Thomas, G.-B. Stan, M. Barahona, *Rigorous bounds on the stationary distributions of the chemical master equation via mathematical programming*, 2017. <https://arxiv.org/abs/1702.05468>

### 17.3 Bosonic Lindblad bootstrap

Robichon and Tilloy construct a hierarchy of semidefinite relaxations giving upper and lower bounds on stationary observables of bosonic open quantum systems.

* G. Robichon, A. Tilloy, *Bootstrapping the stationary state of bosonic open quantum systems*, 2024. <https://arxiv.org/abs/2410.07384>

### 17.4 Quantum moment problems

Noncommutative moment hierarchies provide positive extension and infeasibility certificates for operator moment data.

* A. C. Doherty, Y.-C. Liang, B. Toner, S. Wehner, *The quantum moment problem and bounds on entangled multi-prover games*, 2008. <https://arxiv.org/abs/0803.4373>

### 17.5 BBGKY representability correction

Krönke and Schmelcher study high order truncated bosonic BBGKY hierarchies, representability violations, dynamical purification, and minimally invasive equation corrections.

* S. Krönke, P. Schmelcher, *The BBGKY hierarchy for ultracold bosonic systems*, 2017. <https://arxiv.org/abs/1712.00819>

### 17.6 Variational moment closure

Bronstein and Koeppl derive moment closures variationally and connect them to entropic matching, clarifying closure failure modes and principled reconstruction families.

* L. Bronstein, H. Koeppl, *A variational approach to moment-closure approximations for the kinetics of biomolecular reaction networks*, 2017. <https://arxiv.org/abs/1709.02963>

### 17.7 Entropy based realizability

Entropy based moment methods use a realizable ansatz and convex optimization. Regularization controls failures near the realizability boundary.

* G. W. Alldredge, M. Frank, C. D. Hauck, *A regularized entropy-based moment method for kinetic equations*, 2018. <https://arxiv.org/abs/1804.05447>

### 17.8 Invariant manifold methods

Gorban and Karlin formulate closure as the construction of an approximately invariant low dimensional manifold and solve the invariance equation using Newton corrections and thermodynamic projection.

* A. N. Gorban, I. V. Karlin, *Method of invariant manifold for chemical kinetics*, 2002. <https://arxiv.org/abs/cond-mat/0207231>

### 17.9 Quantum exponential projection

Exponential quantum projection filters demonstrate finite dimensional projection onto positive exponential families of density operators.

* Q. Gao, G. Zhang, I. R. Petersen, *An Exponential Quantum Projection Filter for Open Quantum Systems*, 2017. <https://arxiv.org/abs/1705.09114>

### 17.10 Entropic gradient structure of Lindblad dynamics

For detailed balance Lindblad generators, relative entropy gradient formulations provide a basis for thermodynamically consistent projection.

* M. Mittnenzweig, A. Mielke, *An entropic gradient structure for Lindblad equations and couplings of quantum systems to macroscopic models*, 2016. <https://arxiv.org/abs/1609.05765>

### 17.11 Carleman finite sections

Carleman finite section analysis shows how truncations of infinite linear hierarchies can converge under analyticity and stability assumptions, including some infinite time results. The setting is the reverse of cumulant closure but the distinction between finite time and stationary convergence is directly relevant.

* A. Amini, C. Zheng, Q. Sun, N. Motee, *Carleman Linearization of Nonlinear Systems and Its Finite-Section Approximations*, 2022. <https://arxiv.org/abs/2207.07755>

### 17.12 Numerical continuation

General polynomial homotopy and pseudoarclength continuation remain important candidate generation tools. They do not supply physical selection by themselves. HomotopyContinuation.jl and BifurcationKit.jl should be treated as optional reference backends behind package interfaces.

### 17.13 Mori-Zwanzig and optimal prediction

Projection of an exact high dimensional dynamics onto a finite set of resolved observables generally produces memory and unresolved forcing. This explains why an autonomous finite cumulant closure can change the stationary structure and motivates auxiliary memory variables or invariant manifold corrections.

* A. J. Chorin, O. H. Hald, R. Kupferman, *Optimal prediction and the Mori-Zwanzig representation of irreversible processes*, 2000. <https://doi.org/10.1073/pnas.97.7.2968>
* A. J. Chorin, O. H. Hald, R. Kupferman, *Optimal prediction with memory*, 2002. <https://doi.org/10.1016/S0167-2789(02)00446-3>

### 17.14 Koopman invariant observable spaces

Finite exact closure of observables requires an invariant observable subspace, a restrictive condition for generic nonlinear dynamics. This motivates measuring and reducing the invariance defect of the retained operator space.

* S. L. Brunton, B. W. Brunton, J. L. Proctor, J. N. Kutz, *Koopman invariant subspaces and finite linear representations of nonlinear dynamical systems for control*, 2015. <https://arxiv.org/abs/1510.03007>

### 17.15 Markov lumpability

Lumpability is the classical exact reduction condition under which projected variables retain a closed Markov description. It is a useful analogue for identifying rare exact finite reductions of Liouvillian or moment dynamics.

* B. C. Geiger, C. Temmel, *Lumpings of Markov chains, entropy rate preservation, and higher-order lumpability*, 2012. <https://arxiv.org/abs/1212.4375>

## 18. Acceptance criteria for the first research release

A first research release is successful when it can:

1. Import stationary hierarchies from QuantumCumulants.jl.
2. Preserve symbolic mappings, adjoints, symmetries, and closure metadata.
3. Generate stationary candidates with at least two independent methods.
4. Lift a candidate through at least three successive hierarchy orders using `\Phi_n`.
5. Report stationary residuals, held out residuals, shared moment errors, and order lift defects.
6. Evaluate configurable moment matrix positivity tests.
7. Cluster physically plausible branches on selected observables.
8. Construct at least one stationary moment envelope or positive extension test on a benchmark.
9. Compare with exact or high accuracy stationary observables on at least two models.
10. Demonstrate one case where algebraic roots remain multiple but physical predictions collapse.
11. Demonstrate one case where the hierarchy does not collapse and report the failure honestly.
