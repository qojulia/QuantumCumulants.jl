# Technical specification

## 1. Problem statement

Given a stationary truncated hierarchy generated from an open quantum system,

```math
F_n(x; p) = 0,
```

compute one or more candidate roots that are connected to a known physical reference state and rank them using information inherited from the original Liouvillian and from neighbouring hierarchy orders.

The solver must distinguish three concepts:

1. Algebraic solution: a root of the finite polynomial system.
2. Dynamically selected solution: a stable fixed point of the truncated equations or a root reached by a chosen continuation path.
3. Physically plausible solution: a root that satisfies available adjoint, symmetry, realizability, and cross order consistency tests.

These concepts need not coincide at finite order.

## 2. Design goals

### 2.1 Primary goals

1. Track a physically initialized stationary branch across model parameters.
2. Lift a stationary solution from truncation order `n` to order `n + 1` using the exact closure map.
3. Generate a minimal real polynomial system from QuantumCumulants expressions.
4. Report branch ambiguity rather than hiding it.
5. Provide physical diagnostics and optional certification interfaces.
6. Exploit sparse hierarchy structure and generated Jacobians.
7. Remain useful when complete root enumeration is impossible.

### 2.2 Non goals for the first release

1. Complete enumeration of all complex roots at high order.
2. Replacement of HomotopyContinuation.jl.
3. A general proof that cumulant closures converge.
4. A claim that finite order positivity tests prove existence of a density operator.
5. Full semidefinite optimization support in the core package.

## 3. Mathematical model

### 3.1 Stationary hierarchy

At truncation order `n`, the hierarchy is represented as

```math
F_n : \mathbb R^{d_n} \times \mathbb R^q \rightarrow \mathbb R^{d_n}.
```

The real formulation is preferred for the package interface because physical moment variables obey conjugation relations. A complex polynomial backend may be supported through an internal complexification.

### 3.2 Realification

For every non Hermitian moment pair,

```math
m_O = q_O + i p_O,
\qquad
m_{O^\dagger} = q_O - i p_O.
```

Hermitian moments receive one real variable. Equations are split into real and imaginary parts. Redundant adjoint equations are removed deterministically.

The realification layer must preserve a bidirectional map between numerical coordinates and symbolic moments.

### 3.3 Closure map

Let `x` contain variables retained at order `n`, and let `z` contain variables first retained at order `n + 1`. QuantumCumulants provides or implies a closure map

```math
\Phi_n : x \mapsto z_{
\mathrm{closure}}.
```

The natural new coordinates are

```math
u = z - \Phi_n(x).
```

At `u = 0`, the newly retained moments agree exactly with the lower order closure.

## 4. Continuation methods

### 4.1 Parameter continuation

For a path `p(t)`, solve

```math
H(x,t) = F_n(x; p(t)) = 0.
```

A predictor corrector step uses

```math
J_x \dot{x} = -\partial_t H
```

followed by damped Newton correction at the new path parameter.

Required features:

1. Adaptive step size.
2. Sparse direct linear solves initially.
3. Jacobian condition estimates.
4. Automatic variable and equation scaling.
5. Mixed precision fallback as a later feature.
6. Optional pseudoarclength continuation.

### 4.2 Reset regularized path

Define

```math
\mathcal L_{\lambda,\epsilon}
= \lambda \mathcal L_{\mathrm{target}}
+ \epsilon \mathcal R_\sigma,
```

with

```math
\mathcal R_\sigma(\rho)
= \operatorname{Tr}(\rho)\sigma - \rho.
```

A default two stage path is

```math
(\lambda,\epsilon):
(0,\epsilon_0)
\rightarrow
(1,\epsilon_0)
\rightarrow
(1,0).
```

This gives an exact reference state and regularizes the exact Liouvillian. It does not guarantee absence of closure induced singularities, so the tracker must still monitor conditioning and branch ambiguity.

### 4.3 Order continuation

Use the padded lower order start system

```math
G_0(x,u)
=
\begin{pmatrix}
F_n(x)\\
D u
\end{pmatrix},
```

where `D` is a scaling or damping matrix.

Write the order `n + 1` target equations in closure deviation coordinates,

```math
G_1(x,u)
=
F_{n+1}(x, \Phi_n(x) + u).
```

A first implementation may use

```math
H(x,u,t)
= (1-t)G_0(x,u) + tG_1(x,u).
```

Alternative structured homotopies should be benchmarked. In particular, the old equations and equations for new variables may be activated separately.

### 4.4 Predictor for newly retained moments

The closure consistent point is

```math
u = 0.
```

A stronger predictor solves a regularized linearized least squares problem for `Δu`,

```math
\Delta u
= \arg\min_v
\left\|W(r + J_u v)\right\|_2^2
+ \gamma\|S v\|_2^2.
```

This corrects the newly introduced moments before full continuation begins.

## 5. Branch management

Tracking one path does not establish that the endpoint is canonical. The package must represent ambiguity explicitly.

### 5.1 Branch beam

Maintain a small set of candidate paths. Spawn alternatives when:

1. The smallest singular value of the Jacobian drops below a threshold.
2. Newton correction requires excessive iterations.
3. The predictor error grows sharply.
4. Physical diagnostics change discontinuously.
5. Different admissible parameter paths yield different endpoints.

### 5.2 Path dependence tests

Optional tests include:

1. Two different physical parameter paths.
2. Two different reset reference states.
3. A closed loop in parameter space.
4. Forward and reverse tracking.
5. A small complex parameter detour when supported by the backend.

The result status must distinguish a unique observed branch from an unresolved branch selection problem.

## 6. Physical diagnostics

### 6.1 Stationary residual

```math
r_n = \|W_n F_n(x_n)\|.
```

### 6.2 Adjoint consistency

When adjoint relations are not eliminated exactly,

```math
r_{\dagger}
= \max_O
\left|m_{O^\dagger} - m_O^*\right|.
```

### 6.3 Order lift defect

```math
\delta_n
= \left\|W_{n+1}F_{n+1}(x_n,\Phi_n(x_n))\right\|.
```

Residuals should also be separated into equations inherited from order `n` and equations first introduced at order `n + 1`.

### 6.4 Shared moment convergence

For fixed projection order `k`,

```math
e_{n,k}
= \left\|S_k\left(\pi_kx_n-\pi_kx_{n-1}\right)\right\|.
```

The package should report convergence separately for odd and even truncation sequences when parity effects are present.

### 6.5 Held out stationary equations

Evaluate

```math
\langle\mathcal L^\dagger O\rangle
```

for operators not used to define the closed stationary system. These provide out of sample stationarity tests.

### 6.6 Moment matrix positivity

For an operator basis `B`, construct

```math
M_B(x)_{ij}
= \langle B_i^\dagger B_j\rangle_x.
```

Report

```math
\mu_B = \lambda_{\min}(M_B)
```

and total negative eigenvalue weight. Passing a finite collection of tests is necessary evidence, not a proof of full realizability.

### 6.7 Stability metadata

For the truncated dynamical Jacobian,

```math
\alpha_n
= \max \operatorname{Re}\operatorname{eig}J_n(x_n).
```

This classifies local stability of the closed dynamics. It is not a mandatory physicality criterion.

## 7. Optional semidefinite validation

A validator may construct a convex feasibility problem using raw moments `y`:

```math
A_n y = 0,
\qquad
y(1)=1,
\qquad
M_q(y) \succeq 0.
```

Possible operations:

1. Bound selected observables.
2. Test whether candidate low order moments admit a positive extension.
3. Project a slightly inconsistent candidate onto the feasible pseudomoment set.
4. Reject a branch outside certified observable intervals.

This functionality should be optional and backend agnostic.

## 8. Scaling and sparse structure

### 8.1 Variable scaling

Each variable receives a scale based on one or more of:

1. Operator order.
2. Reference state magnitude.
3. Previous hierarchy order.
4. Estimated occupation number.
5. User supplied scales.

### 8.2 Equation scaling

Each residual row receives a scale based on:

1. Coefficient norm.
2. Local Jacobian row norm.
3. Physical damping rate.
4. User supplied weights.

### 8.3 Hierarchy block structure

Variables and equations are grouped by operator or cumulant order. The interface should expose this block structure for sparse direct solvers and future block preconditioners.

## 9. Proposed API

```julia
abstract type AbstractStationaryMethod end

struct PhysicalContinuation <: AbstractStationaryMethod
    tracker
    scaler
    branch_policy
    validators
end
```

Primary entry points:

```julia
stationary_state(hierarchy; kwargs...)
stationary_sequence(hierarchies; kwargs...)
continue_parameters(problem, path; kwargs...)
continue_order(lower, higher, solution; kwargs...)
validate_stationary(solution, hierarchy; kwargs...)
```

Core types:

```julia
StationaryPolynomialSystem
RealifiedCoordinates
ClosureMap
ParameterPath
ResetPath
OrderLift
TrackedBranch
BranchEnsemble
PhysicalDiagnostics
StationarySolution
StationarySequence
```

### 9.1 Stationary solution result

```julia
struct StationarySolution
    values
    moments
    parameters
    order
    residual
    adjoint_residual
    order_lift_defect
    positivity_margins
    stability_abscissa
    path_condition
    branch_id
    status
end
```

Suggested status values:

```julia
:converged
:converged_with_warnings
:ambiguous_branch
:positivity_violation
:singular_endpoint
:tracking_failed
:insufficient_validation
```

## 10. Backend strategy

### 10.1 Initial backends

1. HomotopyContinuation.jl for polynomial path tracking and low order root enumeration.
2. BifurcationKit.jl for pseudoarclength continuation and stability.
3. SparseArrays and LinearSolve.jl for custom predictor corrector prototypes.

### 10.2 Custom tracker criterion

A custom low level tracker should only be implemented after profiling identifies general tracker overhead as a dominant cost. The initial novelty should remain in hierarchy construction, prediction, branch policy, and physical validation.

## 11. Verification strategy

### 11.1 Unit tests

1. Correct realification of adjoint pairs.
2. Correct reconstruction of complex moments.
3. Closure consistent order lift.
4. Jacobian agreement with finite differences.
5. Predictor corrector convergence on known polynomial systems.
6. Moment matrix construction.
7. Invariance under deterministic variable reordering.

### 11.2 Reference models

1. Driven Kerr resonator with known stationary moments.
2. Finite spin models with exact Liouvillian nullspace solutions.
3. Parametric oscillator or dissipative cat model.
4. Permutation symmetric light matter model.
5. A classical stochastic reaction model with known moment closure pathologies.

### 11.3 Comparative methods

Compare against:

1. Random start Newton solves.
2. Long time integration of the truncated ODE.
3. HomotopyContinuation.jl parameter homotopy.
4. BifurcationKit continuation.
5. Exact finite Hilbert space solutions where available.
6. Semidefinite observable bounds where available.

### 11.4 Metrics

1. Error in selected observables.
2. Runtime and memory.
3. Number of tracked branches.
4. Corrector iterations and rejected steps.
5. Minimum Jacobian singular value.
6. Positivity margin.
7. Order lift defect.
8. Shared moment convergence.
9. Path dependence.
10. Failure classification quality.

## 12. Research hypotheses

### Hypothesis A

For a primitive open quantum system and a suitable physical reference path, a hierarchy aware continuation method identifies a sequence of stationary roots whose fixed low order moments converge over a useful parameter regime.

### Hypothesis B

Closure consistent order prediction substantially reduces continuation failures compared with zero initialization or independent random start solves.

### Hypothesis C

Moment positivity and held out stationarity tests reject a significant fraction of algebraically valid but physically spurious roots.

### Hypothesis D

Branch ambiguity correlates with small Jacobian singular values, metastability, and regions where the cumulant hierarchy converges nonuniformly.

These are empirical hypotheses. They should not be presented as general theorems without additional assumptions.

## 13. Main risks and mitigations

### Risk: path dependent endpoint

Mitigation: branch beams, alternate paths, loop tests, and explicit ambiguity reporting.

### Risk: closure artifacts remain stable and positive at tested order

Mitigation: next order defects, held out equations, higher moment extension tests, and comparison across several orders.

### Risk: poor scaling at high order

Mitigation: automatic order aware scaling, sparse factorization reuse, and block preconditioning.

### Risk: singular or nearly singular branch

Mitigation: adaptive steps, pseudoarclength continuation, precision escalation, and multiple predictors.

### Risk: positivity tests become too expensive

Mitigation: use operator bases selected from the moment graph and run SDP validation only at selected endpoints.

### Risk: infinite dimensional moments are not tight

Mitigation: support user supplied energy or number bounds and report when no growth control is available.

## 14. Acceptance criteria for a first release

A first release is successful when it can:

1. Import stationary hierarchies from QuantumCumulants.jl.
2. Realify them without redundant equations.
3. Track a known stationary solution through a physical parameter path.
4. Lift that solution through at least three successive hierarchy orders.
5. Report stationary residuals, shared moment errors, and order lift defects.
6. Evaluate at least one configurable moment matrix positivity test.
7. Reproduce exact or high accuracy reference observables on two benchmark models.
8. Detect and report at least one deliberately constructed ambiguous or nonphysical case.
