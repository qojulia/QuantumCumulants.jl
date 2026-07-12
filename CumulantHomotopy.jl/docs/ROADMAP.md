# Roadmap

The roadmap is organized around the scientific problem rather than around one numerical backend. The first milestones establish the hierarchy representation, physical diagnostics, and benchmark evidence. Continuation and homotopy methods are then added as interchangeable candidate generators.

## Phase 0: formal problem and benchmark harness

1. Define algebraic roots, physically filtered roots, stationary moment envelopes, and collapse metrics.
2. Select benchmark observables and tolerances before implementing branch selection.
3. Build exact or high accuracy reference solvers for a driven Kerr resonator and one finite spin model.
4. Generate stationary cumulant equations at several orders.
5. Record the number and structure of variables and equations, closure maps, and operator orders.
6. Implement a common experiment output format for candidates, diagnostics, and reference errors.

Deliverable: reproducible notebooks showing the stationary root structure and exact reference observables without assuming a preferred solver.

## Phase 1: hierarchy interoperability

1. Define `StationaryHierarchy` and deterministic variable ordering.
2. Preserve mappings to QuantumCumulants moments and operators.
3. Expose adjoint pairs, Hermitian variables, symmetry sectors, and operator identities.
4. Extract or reconstruct the closure map `\Phi_n` between adjacent hierarchy orders.
5. Expose sparse residuals, Jacobians, and block structure by operator order.
6. Support both complexified and reduced real coordinate policies.

Deliverable: a backend independent stationary hierarchy representation.

## Phase 2: physical diagnostics and collapse metrics

1. Implement stationary residuals and equation scaling.
2. Implement adjoint, symmetry, and conservation diagnostics.
3. Generate held out stationary equations from the moment graph and `\mathcal L^\dagger`.
4. Implement closure consistent next order defects.
5. Implement shared moment convergence across several orders.
6. Implement odd and even subsequence diagnostics.
7. Define branch clustering and projected branch diameter.
8. Define result statuses for physical plausibility, ambiguity, and collapse.

Deliverable: physical evidence can be attached to any stationary candidate, independently of how it was generated.

## Phase 3: baseline candidate generators

1. Implement physically initialized Newton solves with NonlinearSolve.jl.
2. Implement nonlinear least squares for overdetermined systems.
3. Implement long time integration followed by root polishing.
4. Add multiple physically motivated initial conditions.
5. Compare candidate sets from independent methods.
6. Establish failure classifications and benchmark runtime.

Deliverable: `stationary_candidates` with at least two independent candidate generation methods.

## Phase 4: continuation in parameters and hierarchy order

1. Implement physical parameter paths.
2. Implement reset regularized paths.
3. Construct closure deviation coordinates `u = z - \Phi_n(x)`.
4. Implement direct higher order correction from the closure consistent lift.
5. Implement the regularized linearized predictor for newly retained moments.
6. Implement a native sparse predictor corrector prototype.
7. Add optional HomotopyContinuation.jl tracking as a reference backend.
8. Add path reversal and alternate path tests.

Deliverable: `stationary_sequence` that follows candidates through parameters and at least three successive hierarchy orders.

## Phase 5: branch management

1. Maintain a small ensemble of candidate branches.
2. Detect poor conditioning, large predictor errors, and corrector instability.
3. Spawn additional candidates near ambiguous segments.
4. Add pseudoarclength continuation through BifurcationKit.jl or a native backend.
5. Compare different reset states and physical parameter paths.
6. Perform forward, reverse, and loop consistency tests.
7. Cluster endpoints by selected low order observables.

Deliverable: explicit reporting of unique observed branches, multiple physical clusters, and unresolved path dependence.

## Phase 6: moment realizability

1. Construct moment matrices from configurable operator bases.
2. Select operator bases adaptively from observables and the moment graph.
3. Report minimum eigenvalues and negative eigenvalue weights.
4. Add bosonic energy and occupation constraints.
5. Implement positive extension feasibility for selected endpoints.
6. Compute distance from a cumulant candidate to a realizable pseudomoment set.

Deliverable: physically implausible roots can be rejected or flagged using increasingly strong tests.

## Phase 7: stationary moment envelope

1. Define the convex set `K_{n,q}` from stationary equations, operator identities, positivity, and growth bounds.
2. Compute lower and upper bounds on selected observables.
3. Implement a strictly convex finite order selector.
4. Project cumulant candidates onto the stationary envelope.
5. Compare projections of distinct cumulant branches.
6. Measure envelope widths alongside branch diameters.
7. Add an optional semidefinite backend after selecting a suitable Julia optimization stack.

Deliverable: one reproducible finite order stationary estimate plus quantified remaining ambiguity.

## Phase 8: structure preserving closure experiments

1. Implement a maximum entropy or exponential family reconstruction for small finite models.
2. Implement a variational closure objective and compare it with zero higher cumulants.
3. Compute the invariance defect of the retained observable space under `\mathcal L^\dagger`.
4. Add operators adaptively based on the invariance defect.
5. Test minimal representability corrections inspired by BBGKY purification.
6. Investigate thermodynamic projection for detailed balance models.
7. Investigate auxiliary memory variables for models with slow unresolved modes.

Deliverable: at least one closure that improves stationary collapse relative to the standard cumulant truncation at comparable dimension.

## Phase 9: performance and native solvers

1. Profile equation generation, residual evaluation, Jacobian construction, linear solves, continuation overhead, and semidefinite validation separately.
2. Reuse sparse factorizations between nearby path points.
3. Add matrix free Jacobian vector products.
4. Add hierarchy block preconditioners.
5. Add mixed precision only for difficult segments.
6. Compare the native tracker with HomotopyContinuation.jl and BifurcationKit.jl on the same paths.
7. Replace general backends only where profiling demonstrates a meaningful advantage.

Deliverable: high order stationary sequences that outperform independent root solves while retaining physical diagnostics.

## Phase 10: theory

1. Prove finite dimensional termination when the retained operator basis becomes complete.
2. Formalize nested stationary moment envelopes and their projections.
3. Identify compactness and moment determinacy assumptions for bosonic systems.
4. Establish conditions under which envelope widths converge to zero for a unique stationary state.
5. Study whether physical branch diameters can be bounded by envelope widths and closure defects.
6. Relate order lift defects to errors in selected stationary observables.
7. Characterize when a projected detailed balance dynamics retains a unique entropy minimizer.

Deliverable: precise theorems for restricted model classes and clearly stated conjectures for the general method.

## Proposed first issues

1. Define `StationaryHierarchy` and the candidate result schema.
2. Extract closure maps between adjacent hierarchy orders.
3. Build exact stationary references for Kerr and a finite spin model.
4. Implement held out stationary equation generation.
5. Implement order lift defect reporting.
6. Implement shared moment and branch diameter metrics.
7. Add physically initialized Newton and time integration baselines.
8. Implement closure consistent higher order correction.
9. Implement basic moment matrix positivity checks.
10. Design `K_{n,q}` and an initial small semidefinite feasibility experiment.
11. Add an optional HomotopyContinuation.jl reference backend.
12. Design collapse and ambiguity statuses.

## First publication quality experiment

The first complete study should answer the following questions on at least three models:

1. How many stationary roots are found at each cumulant order?
2. How many survive adjoint, symmetry, positivity, and held out stationarity tests?
3. Do surviving roots form one or several clusters on low order observables?
4. Does the branch diameter decrease with order?
5. Do stationary moment envelope widths decrease with order?
6. Does closure consistent order continuation find the same cluster as independent candidate generators?
7. In which parameter regimes does collapse fail?
8. Which diagnostics predict failure earliest?

The central success criterion is evidence that the finite hierarchy produces a shrinking set of physically admissible stationary predictions, not merely that one numerical path can be followed.
