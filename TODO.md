# TODO



## P1 Symbolics / MTK architecture

- Use CSE deliberately. `Symbolics.build_function(...; cse=true)` creates temporary
  bindings for repeated subexpressions, but the package currently leaves most RHSs
  raw and only simplifies on request. Benchmark CSE for large cumulant-expanded
  systems, especially indexed/evaluated systems, before enabling by default.

- Revisit parameter collection and arrayization with Symbolics APIs. `_collect_params!`,
  `_build_callable_to_array_sub`, `_parse_slot`, and `parameter_map(eqs, pairs)` are
  hand-rolled around Symbolics variable naming and `getindex`. Replace name parsing
  with metadata or Symbolics-supported array variable introspection where possible.

## P1 graph and closure maintainability

- Make `MomentGraph` the single source of semantic truth. Transformations currently
  rebuild graphs from stored equations in some paths and rederive in others. Preserve
  user RHS modifications, filter substitutions, treatment state, source averages,
  and generated variables in one graph-backed representation.

- Cache canonicalization. `canonical_rep`, `_treatment_key`, `symmetric_min`,
  `qadd_order_key`, `acts_on`, `get_order`, and conjugate representatives are repeatedly
  recomputed across closure, scale/evaluate, `System`, `get_solution`, correlation,
  and spectrum. Add per-graph caches keyed by operator identity / structural hash.

- Preserve and validate user-modified equations. `modify_equations` changes stored
  RHSs, but `complete` can rebuild from operators and discard those modifications.
  Either make modifications graph-native or document and test the supported order of
  operations.

- Add closure diagnostics. Users need tooling to explain why a hierarchy grows:
  state-count per iteration, which RHS leaf introduced each new moment, filter hits,
  conjugate folding, and per-subspace treatment effects.

## P2 feature coverage

- Add multi-channel stochastic examples and tests. Current examples mostly monitor
  one nonzero-efficiency channel. Add tests where two channels are monitored and
  verify independent-noise structure and deterministic ensemble consistency.

- Add performance baselines to the repo or CI artifacts. Benchmarks exist, but a
  maintainer should be able to see current medians/allocation trends for meanfield,
  completion, scale, evaluate, System construction, correlation and spectrum.

- Add allocation-focused tests or benchmarks for symbolic hot spots. Track
  allocations in cumulant expansion, canonicalization, `evaluate`, `scale`,
  `System`, and `_spectrum_kernel`.

- Update docs that drifted from implementation. `docs/src/correlation.md` still
  describes the spectrum path as generating numerical functions with
  `Symbolics.build_function`; the implementation currently classifies leaves and
  substitutes numerically. Either update the docs or implement the generated-kernel
  path.

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
