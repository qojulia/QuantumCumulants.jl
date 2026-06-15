# TODO

## P2 feature coverage

- Add allocation-focused tests or benchmarks for symbolic hot spots. Track
  allocations in cumulant expansion, canonicalization, `evaluate`, `scale`,
  `System`, and `_spectrum_kernel`.

- Update docs that drifted from implementation.

- Document the MTK bridge invariants. Explain why averages are not currently used
  directly as MTK unknowns, how conjugate folding is represented, and how users
  should reason about `get_adjoints=false`.

## P2 API and dependency hygiene

- Improve error messages for unsupported cases: backward retrodiction with
  nondiagonal measurements, unsupported indexed coefficient arity, unclosed systems
  at `System` construction, and symbolic limits missing from `evaluate`.

- Keep in-place and non-mutating APIs semantically identical. Add tests for every
  `foo` / `foo!` pair that compare equations, states, operators, treatments, and
  subsequent `System` behavior.
