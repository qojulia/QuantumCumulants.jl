# TODO

## P2 API and dependency hygiene

- Improve error messages for unsupported cases: backward retrodiction with
  nondiagonal measurements, unsupported indexed coefficient arity, unclosed systems
  at `System` construction, and symbolic limits missing from `evaluate`.

- Keep in-place and non-mutating APIs semantically identical. Add tests for every
  `foo` / `foo!` pair that compare equations, states, operators, treatments, and
  subsequent `System` behavior.
