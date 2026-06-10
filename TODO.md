# TODO

This is the active package audit backlog. The package goal is to be a performant,
maintainable symbolic moment-layer for deriving, closing, reducing, compiling and
solving generalized mean-field equations for open quantum systems.

Resolved historical investigations (collective indexed dissipation / NE-constraint
flow, det/stoch closure asymmetry, conjugate dedup, ...) are archived in
[INVESTIGATION_NOTES.md](INVESTIGATION_NOTES.md).

## P0 correctness and reproducibility

- Fix `scale!` treatment metadata. `scale(eqs)` returns an equation set whose
  `treatments` mark scaled subspaces as `Scaled`, but `scale!(eqs)` only replaces
  vector fields through `_replace_contents!` and leaves `eqs.treatments` unchanged
  (`Free`). This makes later canonical lookup / `System` / `get_solution` operate
  with stale semantic state. Add a regression test comparing `scale(eqs).treatments`
  and `scale!(copy).treatments`.

- Fix scalar `order::Int` normalization. `_normalize_order(order::Int, eqs)` sizes
  the order vector from `eqs.hamiltonian` only. If `ops` or jumps touch subspaces
  not touched by `H`, `meanfield(...; order=2)` can build an order vector that is
  too short and later throw a `BoundsError` in `cumulant_expansion`. Normalize over
  the full system inputs: observables, Hamiltonian, jumps and adjoint jumps.


- Rework stochastic systems to support independent measurement channels. The API
  accepts `efficiencies` per jump, but the current SDE bridge creates one Brownian
  `_qc_dW` and sums all noise drifts into a single column. Multiple monitored
  channels should generally produce independent Wiener processes. Represent
  `noise_equations` as channel columns or store per-channel noise drifts in the graph,
  then generate one Brownian per nonzero independent channel.

- Enforce `Spectrum` preconditions. The docs say the Laplace-transform spectrum is
  the steady-state method, but `Spectrum(c::CorrelationFunction)` does not reject or
  warn when `c.steady_state == false`. Add an explicit check or a separate
  non-steady-state API.

## P1 Symbolics / MTK architecture

- Replace string-based moment serialization as the identity bridge. `average(op)` is
  already a `BasicSymbolic` average leaf (`operation == avg`, `symtype == AvgSym`),
  but `_state_registry` throws that identity away and mints fresh variables named
  by `serialize(op)`, e.g. `avg_a(t)`. Directly using `average(a)` as an MTK
  unknown currently fails because it is not a function of `t`, so the answer is not
  to pass average leaves blindly to MTK. Instead, create a first-class
  `MomentVariableRegistry` that stores:
  - the canonical average leaf / operator representative;
  - the generated time-dependent Symbolics variable;
  - conjugation side information;
  - treatment state (`Free` / `Scaled` / `Concrete`);
  - stable metadata linking the MTK variable back to the average.
  Keep generated variable names purely as display/debug names, never as identity.

- Investigate an SQA / Symbolics extension for time-dependent averages. A stronger
  long-term design would make averaged operators usable as dependent variables, e.g.
  a symbolic `avg(op)(t)` or metadata-enhanced variable whose source is the average
  leaf. This likely belongs at the SQA/QC boundary. Prototype before committing,
  because SQA's current `AvgSym` is deliberately an average term, not an MTK unknown.

- Centralize tree rewriting. The package currently has several ad hoc walkers:
  `mapleaves`, `_filter_expr`, `_scale_expr`, `_materialise_walk`,
  `_apply_callable_walk`, `_collect_params!`, spectrum substitution, etc. Replace
  them with one small expression-rewrite utility layer that preserves metadata,
  handles `complex(re, im)` consistently, and makes "do not descend into average
  leaves" an explicit option.

- Use `Symbolics.build_function` for generated numerical kernels. The local
  Symbolics version supports `build_function(...; cse=true)` and emits both
  out-of-place and in-place functions. Candidate uses:
  - spectrum kernels (`A`, `b`, ambient substitutions) instead of repeated
    `SymbolicUtils.substitute` and `_scalarize`;
  - RHS evaluation benchmarks independent of MTK compilation;
  - optional fast paths for closed systems where users want a raw ODEFunction-like
    kernel.

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
  `_serialize`, `acts_on`, `get_order`, and conjugate representatives are repeatedly
  recomputed across closure, scale/evaluate, `System`, `get_solution`, correlation,
  and spectrum. Add per-graph caches keyed by operator identity / structural hash.

- Replace string serialization in canonical comparisons where possible. `_serialize`
  currently compares stringified operators to pick representatives. This is simple
  but expensive and fragile. Prefer a structural ordering/key from SQA or a cached
  normalized tuple representation.

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
