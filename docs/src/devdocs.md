# Developer Documentation

This page describes the internal architecture of **QuantumCumulants.jl**. It is aimed at
contributors and at users who want to understand what happens between
[`meanfield`](@ref) and a numerical solution. The operator algebra itself
(`QSym`/`QAdd`, normal ordering, indices, sums) lives in
[SecondQuantizedAlgebra.jl](https://github.com/qojulia/SecondQuantizedAlgebra.jl); see its
[developer documentation](https://qojulia.github.io/SecondQuantizedAlgebra.jl/stable/devdocs/)
for that layer. QuantumCumulants is the *moment* layer on top: it turns the algebra into a
closed set of c-number differential equations.

## Pipeline overview

The source files are loaded in dependency order (`src/QuantumCumulants.jl`):

| Layer | Files | Role |
|-------|-------|------|
| Containers & identity | `equations.jl`, `tree.jl`, `canonical.jl` | equation structs, leaf traversal, moment keys |
| Algebra to moments | `operator_drift.jl`, `cumulant.jl`, `moments.jl` | Heisenberg RHS, cumulant expansion, node derivation |
| Hierarchy | `graph.jl` | the coupled moment equations as a graph |
| Orchestration | `meanfield.jl`, `completion.jl`, `scaling.jl`, `evaluate.jl`, `mtk.jl`, `correlation.jl` | public passes over the graph |
| Display | `printing.jl` | plain-text and LaTeX |

Every public entry point is a pass that builds, extends, or rewrites one central data
structure: the **moment graph**.

## The moment graph

`MomentGraph` (`graph.jl`) is the in-memory representation of a (possibly incomplete)
cumulant hierarchy:

```julia
struct MomentGraph{S <: SystemSpec}
    nodes::OrderedDict{QAdd, NodeData}        # one entry per tracked moment
    sys::S                                    # Hamiltonian, jumps, rates, order, direction, …
    ctx::CanonCtx                             # index vocabulary + subspace treatments
    treatments::Dict{Int, SubspaceTreatment} # per-subspace: Free / Scaled / Concrete
end
```

Each node is keyed by a canonical-representative `QAdd` (the moment's operator) and holds a
`NodeData` (`moments.jl`) carrying the drift (the averaged, cumulant-expanded right-hand
side), the operator-level RHS, and, when measurement backaction is requested, the noise
drift. The `SystemSpec` (`graph.jl`) is the immutable bundle of everything the derivation
needs: `hamiltonian`, `jumps`, `jumps_dagger`, `rates`, `efficiencies` (`nothing` means
deterministic), `iv`, `order`, `mix_choice`, and `direction`.

## Deriving a node

`derive(op, sys, ctx)` (`moments.jl`) produces the equation of motion for a single moment:

1. `_operator_rhs` (`operator_drift.jl`) builds the operator-level RHS. `Forward` uses
   ``+i[H, \cdot]`` plus the Lindblad recycling
   ``\sum_i r_i\,(J_i^\dagger \cdot J_i - \tfrac12\{J_i^\dagger J_i, \cdot\})``; `Backward`
   flips the commutator sign, swaps ``J \leftrightarrow J^\dagger`` in the recycling, and
   adds the trace-preserving term (retrodiction).
2. `average_and_truncate` averages the RHS and applies the [`cumulant_expansion`](@ref) to
   the requested `order`, redistributing the SQA sum-scope metadata
   (`SumIndices`/`SumNonEqual`) onto each factorised product so a later `scale` can recover
   the per-index range prefactor.

The leaves of the resulting drift (`eachleaf`, `tree.jl`) are the moments this node couples
to, the input that drives hierarchy closure.

## Cumulant expansion and truncation order

[`cumulant_expansion`](@ref) (`cumulant.jl`) factorises an average of order ``n`` into
products of lower-order averages by setting the joint cumulant of the truncation order to
zero (see the [theory section](@ref theory)). [`get_order`](@ref) measures a term's order
(numbers 0, a single operator 1, a product the number of factors).

`order` is normalised to a `TruncOrder` (a `Vector{Int}`, one cap per Hilbert subspace, or
`nothing`). A term acting on several subspaces combines its per-subspace caps through
`mix_choice` (default `maximum`).

## Moment keys and canonicalisation

The heart of the moment layer is `canonical.jl`. A *key* assigns each average a canonical
label that is identical for two averages exactly when they are the same expectation value,
so the hierarchy closes on a non-redundant set. Two averages can coincide through:

- relabelling of free atom indices,
- the permutation symmetry ``S_n`` of identical atoms, and
- Hermitian conjugation, ``\langle A^\dagger\rangle = \langle A\rangle^*``.

`CanonCtx` holds the index vocabulary used for relabelling (always derived from the user's
declared indices, never freshly minted prefixes) and the set of symmetric/selected
subspaces. Each subspace carries a `SubspaceTreatment` (`equations.jl`):

- `Free`: its atom index stays a free symbolic index (the default after `meanfield`).
- `Scaled`: reduced under permutation symmetry (after [`scale`](@ref)).
- `Concrete`: pinned to fixed sites ``1..M`` (after [`evaluate`](@ref)).

`canon_key` keys at the `Free`/`Concrete` level; `scaled_key` additionally quotients by
permutation symmetry. `canonical_rep` returns both the representative and which side of a
conjugate pair an operator sits on, which is how a single stored moment stands in for both
``\langle A\rangle`` and ``\langle A^\dagger\rangle``.

## Closing the hierarchy

`seed` (`graph.jl`) derives the user's requested operators into the initial nodes.
`closure!` then repeatedly pops a node, walks its drift leaves, and derives any moment not
yet tracked, until no new moments appear. A moment and its conjugate are one physical
unknown: if the conjugate partner is already tracked, the new moment is covered
(`get_adjoints=false` keeps a single representative per pair and recovers the partner via
`conj` at code generation). The `max_iter` guard is a runaway backstop, not a truncation
limiter; hitting it errors rather than silently returning a non-closed system.

[`complete!`](@ref) is `closure!` exposed publicly: it rebuilds a graph from an existing
equation set (`_graph_from_eqs`) and extends it. A `filter_func` can drop unwanted moments
(e.g. phase-invariant or ancilla terms).

## Passes over the graph

All three reductions are graph rewrites that re-key the nodes under a new treatment:

- [`scale`](@ref) (`scaling.jl`) calls `quotient`, marking the selected symmetric subspaces
  `Scaled` and merging permutation-equivalent moments. The `h::Vector{Int}` keyword
  restricts the quotient to specific `space_index` values; empty means all symmetric
  subspaces.
- [`evaluate`](@ref) (`evaluate.jl`) calls `specialize`, unrolling each targeted symbolic
  range to its concrete integer size (from `limits`) and pinning indices to `Concrete`
  sites.
- Because both are keyed by treatment rather than a hardcoded representative, a system can
  be partially scaled and partially evaluated, in any order (hybrid systems).

## Noise and evolution direction

When `efficiencies` is supplied, `meanfield` returns a [`NoiseMeanfieldEquations`](@ref):
each node additionally stores a noise drift built per evolution direction. `Forward` and
`Backward` are singleton tags (`equations.jl`) dispatched at compile time, selecting the
forward measurement-backaction term or the backward (retrodiction) one. Everything
downstream (`complete`, `scale`, `evaluate`) carries the noise column through unchanged. See
[Noise & measurement backaction](@ref) for the user-facing workflow.

## Building the numerical system

`mtk.jl` bridges to ModelingToolkitBase. `_state_registry` assigns each state a `u(t)`
variable and builds two lookup maps: `by_rep` (keyed by conjugate representative, so a leaf
resolves to either the variable or its `conj`) and `by_canon`. `System(eqs; name)`
substitutes every drift leaf through these maps, emitting the symbolic `conj` node via
`SymbolicUtils.term(conj, var; type=Number)` so it survives `mtkcompile`. The LHS becomes
`Differential(iv)(u(t))` only here; the equation struct itself stores the raw `Average`.

[`initial_values`](@ref) computes a `u0` from a numeric `Ket`/density operator (via SQA's
`numeric_average`), [`parameter_map`](@ref) expands indexed/array parameters to the
compiled system's shapes, and [`get_solution`](@ref) reverses the registry to evaluate any
operator-average trajectory, including products that never appeared as a stored state.

## Equation containers

`equations.jl` defines the user-facing structs. [`MeanfieldEquations`](@ref) and
[`NoiseMeanfieldEquations`](@ref) share the supertype `AbstractMeanfieldEquations` and store
the averaged equations, the operator-level equations, the states/operators, the
Hamiltonian/jumps/rates, the `iv`, the `order`, the `direction`, and the per-subspace
`treatments` map. They are array-backed snapshots; the graph is the working representation
that produces them and that every pass rebuilds.
