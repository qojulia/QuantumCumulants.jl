# Roadmap

## Phase 0: validation notebook

1. Select a driven Kerr benchmark.
2. Generate stationary equations at orders two through six.
3. Export a minimal real polynomial system.
4. Track from a linear reference point with HomotopyContinuation.jl.
5. Compare against exact stationary moments.
6. Record order lift defects and positivity margins.

Deliverable: a reproducible notebook that demonstrates the full concept without a stable public API.

## Phase 1: hierarchy representation

1. Define `StationaryPolynomialSystem`.
2. Implement deterministic variable ordering.
3. Implement realification and reconstruction.
4. Expose sparse residuals and Jacobians.
5. Preserve mappings to QuantumCumulants moments.

Deliverable: reliable conversion from a QuantumCumulants hierarchy to a stationary numerical system.

## Phase 2: parameter continuation

1. Implement a HomotopyContinuation.jl backend.
2. Add physical parameter paths.
3. Add reset regularized paths.
4. Add scaling and diagnostics.
5. Add path reversal tests.

Deliverable: `stationary_state` for one hierarchy order.

## Phase 3: order continuation

1. Extract the closure map `Φ_n`.
2. Construct closure deviation coordinates.
3. Implement closure consistent order lifts.
4. Add linearized predictors for new moments.
5. Add shared moment convergence and order lift diagnostics.

Deliverable: `stationary_sequence` across several truncation orders.

## Phase 4: branch management

1. Track a small branch beam.
2. Detect poor conditioning and corrector instability.
3. Add alternate parameter paths.
4. Integrate BifurcationKit pseudoarclength continuation.
5. Report branch ambiguity explicitly.

Deliverable: robust operation near folds and closure induced bifurcations.

## Phase 5: physical validation

1. Construct moment matrices from selected operator bases.
2. Add positivity margins and negative eigenvalue diagnostics.
3. Add held out stationary equations.
4. Add optional extension feasibility interfaces.
5. Add optional semidefinite backends.

Deliverable: physical evidence attached to every stationary candidate.

## Phase 6: performance

1. Profile hierarchy generation, residual evaluation, Jacobian construction, and linear solves separately.
2. Reuse sparse factorizations between nearby continuation steps.
3. Add matrix free Jacobian vector products.
4. Add hierarchy block preconditioners.
5. Add precision escalation only for difficult path segments.

Deliverable: high order benchmarks that outperform generic independent stationary solves.

## Proposed first issues

1. Define the stationary system interchange format.
2. Implement realification for adjoint moment pairs.
3. Extract closure maps between adjacent hierarchy orders.
4. Build a Kerr validation notebook.
5. Implement the HomotopyContinuation.jl parameter tracker backend.
6. Implement order lift defect reporting.
7. Implement basic moment matrix positivity checks.
8. Design branch ambiguity statuses and diagnostics.
