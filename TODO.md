# TODO

## P1 graph and closure maintainability

- Make correlation and spectrum graph-native. The core pipeline (`complete`, `scale`,
  `evaluate`, `simplify`, `modify_equations`, `cumulant_expansion`) now reads and rewrites
  the wrapper's `MomentGraph` directly, but `CorrelationFunction` still rebuilds a graph
  from the assembled ancilla equations via `_graph_from_eqs`, and ambient/source averages
  are not yet graph-tracked. Make these graph-native too and drop `_graph_from_eqs`.

- Cache canonicalization. `canonical_rep`, `_treatment_key`, `symmetric_min`,
  `qadd_order_key`, `acts_on`, `get_order`, and conjugate representatives are repeatedly
  recomputed across closure, scale/evaluate, `System`, `get_solution`, correlation,
  and spectrum. Add per-graph caches keyed by operator identity / structural hash.

- Refresh the stored `ctx` vocabulary as closure mints indices. `MomentGraph.ctx` is
  frozen at `meanfield` time, so consumers needing the post-completion index set (e.g.
  `scale`'s coefficient `sym_to_space`) must currently recover it from the graph's moments.
  Either update the ctx as new per-subspace indices are minted, or make that recovery the
  documented contract.

- Add closure diagnostics. Users need tooling to explain why a hierarchy grows:
  state-count per iteration, which RHS leaf introduced each new moment, filter hits,
  conjugate folding, and per-subspace treatment effects.

## P2 feature coverage

- Add multi-channel stochastic examples and tests. Current examples mostly monitor
  one nonzero-efficiency channel. Add tests where two channels are monitored and
  verify independent-noise structure and deterministic ensemble consistency.

- Add allocation-focused tests or benchmarks for symbolic hot spots. Track
  allocations in cumulant expansion, canonicalization, `evaluate`, `scale`,
  `System`, and `_spectrum_kernel`.

- Update docs that drifted from implementation.

- Document the MTK bridge invariants. Explain why averages are not currently used
  directly as MTK unknowns, how conjugate folding is represented, and how users
  should reason about `get_adjoints=false`.

- Add public introspection helpers:
  - `states(eqs)` / `operators(eqs)` / `moments(eqs)` style accessors;
  - `moment_variable_map(eqs)` for the MTK bridge;
  - `closure_report(eqs)` for missing and generated moments;
  - `noise_channels(eqs)` once multi-channel noise is represented.

## P2 API and dependency hygiene

- Improve error messages for unsupported cases: backward retrodiction with
  nondiagonal measurements, unsupported indexed coefficient arity, unclosed systems
  at `System` construction, and symbolic limits missing from `evaluate`.

- Keep in-place and non-mutating APIs semantically identical. Add tests for every
  `foo` / `foo!` pair that compare equations, states, operators, treatments, and
  subsequent `System` behavior.
