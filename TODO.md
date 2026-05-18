# v1 regressions discovered during test porting

These were identified while porting the master test suite. Each item is a
behavior the master branch handled but the v1 rewrite either drops, weakens,
or gets wrong. Listed in rough priority order.

## `scale(::NoiseMeanFieldEquations)` not defined

**Symptom:** master scaled indexed noise systems via `scale(stoch_eqs)`.
v1 only has `scale(::MeanFieldEquations)`.

**Status:** missing. Would need to scale both `eqs.equations` and
`eqs.noise_equations` consistently.

## Deferred from CHANGELOG (not regressions per se, just not yet ported)

- `evaluate(eqs; limits=(N=>n,))` to unroll indexed sums to fixed-N
- `initial_values(eqs, ψ::Ket)` for spin-coherent / Fock state initialisation
- Rate-matrix collective decay (`J = [[J1, J2]]; rates = [[Γ11 Γ12; Γ21 Γ22]]`)
- `meanfield_backward(...)` (use `meanfield(...; direction=Backward())`)
- `modify_equations(eqs, f)` / `translate_W_to_Y(eqs)` (retrodiction helpers)
- `find_operators(h, order)` returning a deduplicated, canonically-ordered
  operator set so closure tests can assert exact counts
