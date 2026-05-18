# v1 regressions discovered during test porting

These were identified while porting the master test suite. Each item is a
behavior the master branch handled but the v1 rewrite either drops, weakens,
or gets wrong. Listed in rough priority order.

## `scale(::NoiseMeanFieldEquations)` not defined

**Symptom:** master scaled indexed noise systems via `scale(stoch_eqs)`.
v1 only has `scale(::MeanFieldEquations)`.

**Status:** missing. Would need to scale both `eqs.equations` and
`eqs.noise_equations` consistently.

## Example regressions surfaced during 1:1 port from master

Each was a working master example; ported syntax loads fine but fails at a
later stage. The `complex(re, im)` literal bug that initially masked these
is fixed (see `_rewrite_complex_literals` in [src/meanfield.jl](src/meanfield.jl)
and the defensive case in `_substitute_conj_avgs` in [src/mtk.jl](src/mtk.jl),
both using `Symbolics.IM` to keep the additive form in SymReal land).

- **[examples/heterodyne_detection.jl](examples/heterodyne_detection.jl)** —
  `TypeError: in Complex, in T, expected T<:Real, got Type{Number}`. Type
  promotion blows up somewhere in the SDE construction. Likely SQA-side.
- **[examples/retrodiction_homodyne.jl](examples/retrodiction_homodyne.jl)** —
  `SDEProblem` rejects the system with "A completed system is required.
  Call `complete` or `mtkcompile` on the system before creating a
  `SDEProblem`." Example needs a `mtkcompile` step (and possibly QC-side
  too: `complete(::NoiseMeanFieldEquations{...,Backward})` may not be
  fully wired).
- **[examples/excitation-transport-chain.jl](examples/excitation-transport-chain.jl)** —
  `UndefVarError: EnsembleProblem`. Example needs `using
  DifferentialEquations` (or `using SciMLBase`); separately, the
  `EnsembleProblem` pathway should be exercised against the rewritten
  noise machinery.
- **[examples/optomechanical-cooling.jl](examples/optomechanical-cooling.jl)** —
  `MethodError: objects of type SecondQuantizedAlgebra.AvgFunc are not
  callable` at `solve()`. Looks like an `AvgFunc` (sum/index-parametrised
  average) survives into the numerics where a scalar `Num` is expected.
- **[examples/single-atom-laser-spectrum.jl](examples/single-atom-laser-spectrum.jl)** —
  `CorrelationFunction` is fundamentally broken. Current code does
  `new_op = op1 * op2; meanfield([new_op], ...)` which derives equations
  for `⟨a†(τ)a(τ)⟩` (same time), NOT the two-time correlation
  `⟨a†(τ)a(0)⟩`. For `op1=a†, op2=a` this produces 14 states instead of
  the expected 2. The fix is the **Quantum Regression Theorem via an
  ancilla subspace**: create `op2_anc` on a fresh subspace index
  `aon_anc = max(existing aons) + 1` using `_embed_on` helpers, form
  `new_op = op1 * op2_anc`, and call `_complete_ancilla!` which only
  derives equations for averages touching the ancilla. A partial attempt
  at this fix already lives in `src/correlation.jl` (written but not
  tested) — read it, it may be close but have bugs. SQA operator API:
  `Destroy(name, space_index, index)`, `SQA.acts_on(op)` returns subspace
  indices. After fixing `CorrelationFunction`, the `Spectrum` Laplace
  computation also needs fixing: `Symbolics.derivative` on the substituted
  RHS fails with 'array expression'/'symtype' errors. Use a
  finite-difference approach instead: `A[i,j] = f(e_j) - f(0)` by
  substituting all states=0, then setting state_j=1.
- **[examples/unique_squeezing.jl](examples/unique_squeezing.jl)** —
  `ArgumentError: Cannot add arguments of different sizes - encountered
  shapes UnitRange{Int64}[] and SymbolicUtils.Unknown(2)`. `scale` shape
  mismatch on this system.
- **[examples/superradiant_laser_indexed.jl](examples/superradiant_laser_indexed.jl)** —
  `MethodError: no method matching *(::Type{Number}, ::Type{Real},
  ::Type{Any})` at `to_system` after `scale`. Type promotion on the
  scaled-equation RHS.

## Deferred from CHANGELOG (not regressions per se, just not yet ported)

- `evaluate(eqs; limits=(N=>n,))` to unroll indexed sums to fixed-N
- `initial_values(eqs, ψ::Ket)` for spin-coherent / Fock state initialisation
- Rate-matrix collective decay (`J = [[J1, J2]]; rates = [[Γ11 Γ12; Γ21 Γ22]]`)
- `meanfield_backward(...)` (use `meanfield(...; direction=Backward())`)
- `modify_equations(eqs, f)` / `translate_W_to_Y(eqs)` (retrodiction helpers)
- `find_operators(h, order)` returning a deduplicated, canonically-ordered
  operator set so closure tests can assert exact counts
