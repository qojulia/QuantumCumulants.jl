# TODO

Open work items for the v1 rewrite. Each entry names the failure mode and where it
shows up so the fix can be verified end to end.

## Numerics / MTK bridge

### `parameter_map` does not unroll `DoubleIndexedVariable` keys

Mapping `Γ(i, j) => Γmatrix` (and `Ω(i, j) => Ωmatrix`) at the `parameter_map`
level does not produce MTK array parameters that `ODEProblem` recognises. The
resulting build fails with:

```
AssertionError: Expected an `Initial` parameter to exist for variable `Γ(i, j)`,
but did not find one.
```

The single-indexed path (`g(i) => gvec`) works; only the doubly-indexed case is
broken. After `evaluate(eqs; limits=N=>N_)` the MTK system surfaces the
unrolled scalars as `i_1, i_2, ...` rather than as an `M_by_M` array, so
`parameter_map` cannot find a matching `arr_by_name[:Γ]`.

Surfaces in: `examples/cavity_antiresonance_indexed.jl` (transmission sweep).

Fix sketch: extend `_collect_named_params!` / `_param_name` in `src/mtk.jl` to
recognise the per-pair scalars synthesised by `evaluate` from a
`DoubleIndexedVariable` and group them under the user-visible name (`:Γ`).
`parameter_map` then needs to flatten a `Matrix` value to the right
acts-on order.

### `IndexedVariable(::Symbol, ::Int)` is not a constructor

After `scale` (or `evaluate`), calling the user closure `g(1)` to refer to the
collapsed coupling errors with:

```
MethodError: no method matching IndexedVariable(::Symbol, ::Int64)
```

SQA only defines `IndexedVariable(::Symbol, ::Index)`. Master accepted integer
indices on scaled systems; v1 forces users into either a hand-rolled
`SymbolicUtils.Sym{...}(:g; type=Real)` or the `parameter_map(eqs, Dict(g(i) => v))`
detour with the closure as the key.

Surfaces in: `examples/superradiant_laser_indexed.jl`,
`examples/filter-cavity_indexed.jl`.

Fix sketch: add an SQA method `IndexedVariable(name, k::Integer)` returning the
scalar `Num` that `scale`/`evaluate` produced, or document the canonical
v1 pattern (only `parameter_map(eqs, Dict(g(i) => v))`) and adjust the
examples accordingly.

### `scale(::NoiseMeanFieldEquations)` not implemented

`scale(eqs_c)` on stochastic mean-field equations bottoms out in
`IndexedVariable(::Symbol, ::Int)` (same symptom as above) and never produces
a scaled equation set. Master propagated indexed couplings transparently.

Surfaces in: `examples/heterodyne_detection.jl`.

Fix sketch: thread the scaling pass through `NoiseMeanFieldEquations` the same
way it works for `MeanFieldEquations`. Likely a missing method dispatch in
`src/scaling.jl`.

### Matrix-valued `rates` rejected by `meanfield`

Passing a full decay-rate matrix to `meanfield`:

```julia
rates = [Γp(c1, c2) for c1 in 1:M, c2 in 1:M]
meanfield(ops, H, J; rates = rates, order = 2)
```

errors with:

```
TypeError: in keyword argument rates, expected AbstractVector, got a value of type Matrix{Num}
```

Master supported this for collective dissipation (rate-matrix Lindblad).

Surfaces in: `examples/waveguide.jl`.

Fix sketch: extend `meanfield` to accept `AbstractMatrix` `rates` and lower it
to the symmetric Lindblad sum
`(1/2) Σ_ij R_ij (2 J_i ρ J_j† . J_i† J_j ρ . ρ J_i† J_j)`.

## Examples currently broken (downstream of the items above)

- `examples/cavity_antiresonance_indexed.jl`: transmission sweep, blocked on
  `parameter_map` for `DoubleIndexedVariable`.
- `examples/filter-cavity_indexed.jl`: `ODEProblem` build, blocked on the
  integer `IndexedVariable` constructor (`[δ(i) for i in 1:M_]`).
- `examples/heterodyne_detection.jl`: deterministic and stochastic solves,
  blocked on `scale(::NoiseMeanFieldEquations)`.
- `examples/superradiant_laser_indexed.jl`: numerical sweep and spectrum,
  blocked on the integer `IndexedVariable` constructor (`g(1)`).
- `examples/waveguide.jl`: derivation step itself, blocked on matrix-valued
  `rates` argument to `meanfield`.
- `examples/retrodiction_homodyne.jl`: depends on `modify_equations` and the
  `Backward()` SDE path. Already commented out of the docs build; relies on
  the same family of fixes above.
