# Internal architecture

This page is for contributors and for users who want to understand what happens inside QuantumCumulants.jl between [`meanfield`](@ref) and a numerical solution. It explains the internal architecture and, more importantly, *why* it is built the way it is. For the user-facing walkthrough of the same pipeline see [The mean-field pipeline](@ref); for the cumulant-expansion mathematics see the [theory page](@ref theory).

## Scope: the moment layer on top of the algebra

QuantumCumulants is a thin **moment layer**. The operator algebra it manipulates (`QSym`/`QAdd` atoms and sums, normal ordering, indices, symbolic sums, the diagonal split, projector completeness, sum metadata) lives in [SecondQuantizedAlgebra.jl](https://github.com/qojulia/SecondQuantizedAlgebra.jl) (SQA), which QuantumCumulants re-exports. Read SQA's [developer documentation](https://qojulia.github.io/SecondQuantizedAlgebra.jl/stable/devdocs/) for that layer; this page does not re-derive it.

The contract between the two is sharp and worth stating once:

> SQA owns every algebraic transformation. `SQA.change_index`, `SQA.Σ`, `*`, `+`, `expand_completeness`, `assume_distinct_index`, and the `SumIndices`/`SumNonEqual` metadata are the only primitives that carry canonicalisation (projector squashing, commuting reorder, the off-diagonal/diagonal split, NE propagation). QuantumCumulants never hand-rolls a `QAdd` and never sorts operators by hand. When QC needs a new index it derives the name from the user's declared vocabulary, never an invented prefix.

Everything below is the *moment* logic that sits on that algebra: turning operators into a closed set of c-number differential equations.

## One central object: the cumulant hierarchy

The whole package is organised around a single observation.

> Every step the user calls (deriving the equations of motion, closing the hierarchy, scaling, evaluating, building the numerical system) acts on **one shared object**: the cumulant hierarchy, held in memory as a graph of coupled moment equations (the **moment graph**). Each moment carries a **canonical label** that is identical for two expectation values exactly when they are physically the same, and a per-subspace **treatment** (each subspace is independently `Free`, `Scaled`, or `Concrete`) sets how strong that identification is. Because every step is just a relabelling of the same object, closing, scaling, evaluation, and hybrid systems all follow from one mechanism.

Most steps have the same shape: take the equations the user holds (a `MeanfieldEquations`), lift them into a working `MomentGraph`, transform the graph, and read a fresh equation set back out. The exception is [`meanfield`](@ref) itself, which has no prior equations to lift: it builds the first graph directly from the list of operators via `seed`.

```
  MeanfieldEquations                         MeanfieldEquations
  (equations in)                             (equations out)
        │                                          ▲
        │  _graph_from_eqs / _graph_from_stored    │  assemble_equations
        ▼                                          │
   ┌──────────────────────────────────────────────────┐
   │                   MomentGraph                    │   the working graph
   │   seed · closure! · quotient · specialize        │   each step transforms it
   └──────────────────────────────────────────────────┘
```

The whole workflow is then just a sequence of these steps:

```
 ops, H, J ──▶ meanfield ──▶ complete ──▶ scale / evaluate ──▶ System ──▶ solve
               (build)       (close)      (reduce, optional)   (numerics)
```

The source files are loaded in dependency order (`src/QuantumCumulants.jl`):

| Layer | Files | Role |
|-------|-------|------|
| Containers & identity | `equations.jl`, `tree.jl`, `canonical.jl` | equation structs, leaf traversal, moment keys |
| Algebra to moments | `operator_drift.jl`, `cumulant.jl`, `moments.jl` | Heisenberg RHS, cumulant expansion, per-moment derivation |
| The hierarchy | `graph.jl` | the cumulant hierarchy as a graph of moment equations |
| Workflow steps | `meanfield.jl`, `completion.jl`, `scaling.jl`, `evaluate.jl`, `mtk.jl`, `correlation.jl` | the operations the user calls |
| Display | `printing.jl` | plain-text and LaTeX |

## The data model: the moment graph

The moment graph (`graph.jl`) is the in-memory form of a (possibly incomplete) cumulant hierarchy.

```
MomentGraph
├─ nodes :: OrderedDict{QAdd ⟶ NodeData}    one entry per tracked moment
│    key   = canonical-representative QAdd   (the moment's identity)
│    value = NodeData
│             ├─ drift     :: Num    averaged + truncated RHS  (leaves couple to other moments)
│             ├─ op_drift  :: QAdd   operator-level RHS (inspection / latex / re-truncation)
│             ├─ noise     :: Num?   measurement-backaction drift (optional)
│             ├─ op_noise  :: Nothing always nothing (operator-level noise is deferred)
│             ├─ order     :: Int    cached get_order
│             └─ aon       :: Vector{Int}   cached acts_on
├─ sys  :: SystemSpec    frozen derivation inputs (see below)
├─ ctx  :: CanonCtx      index vocabulary + which subspaces are symmetric / selected
└─ treatments :: Dict{Int ⟶ SubspaceTreatment}   per-subspace: Free | Scaled | Concrete
```

`SystemSpec` is the immutable bundle of everything `derive` needs: `hamiltonian`, `jumps`, `jumps_dagger`, `rates`, `efficiencies` (`nothing` means deterministic), `iv`, `order`, `mix_choice`, and `direction`. It is frozen so a moment can be derived on demand at any point during a workflow step without threading a dozen arguments through every call.

Two design choices are worth the rationale:

- **Why a graph keyed by moment, not a flat list of equations.** Closing a cumulant hierarchy means repeatedly asking "do I already have an equation for this moment?". That question is a fast lookup only if the key is the moment's *physical identity*, not its written form. The `OrderedDict` keyed by a canonical representative is exactly that lookup; the closure is then just "take one moment, look at the moments it couples to, derive the ones still missing".

- **Why two representations.** The `MomentGraph` is the mutable working copy each step transforms; `MeanfieldEquations`/`NoiseMeanfieldEquations` (`equations.jl`) are the plain array of equations the user holds and inspects. Every step rebuilds that array from the graph via `assemble_equations`. Keeping them separate means the user-facing object stays a simple, stable, indexable list of equations, while the graph carries the extra bookkeeping (`ctx`, `treatments`, the cached per-moment data) that only the workflow steps need.

Every `Int` key in `treatments`, in `CanonCtx.vocab`, and in `CanonCtx.symmetric` is an SQA `space_index`: the position of a Hilbert-space factor (subspace) in the system's `ProductSpace`, as returned by `acts_on`.

## Moment identity and canonicalisation

This is the heart of the moment layer, and it is what makes the hierarchy *close on a non-redundant set*. The job of `canonical.jl` is to assign each average a **key** that is identical for two averages exactly when they are the same expectation value. Two averages can coincide three ways:

3. **Relabelling of free atom indices.** ``\langle igma^{ge}_i\rangle`` and ``\langle igma^{ge}_j\rangle`` are the same moment under a rename of the bound index.

3. **The permutation symmetry ``S_n`` of identical atoms.** In a scaled ensemble every atom is interchangeable, so ``\langle igma_i igma_j\rangle`` collapses onto one representative.

3. **Hermitian conjugation**, ``\langle A^\dagger\rangle = \langle A\rangle^*``.

`CanonCtx` holds the vocabulary used for the relabelling (always derived from the user's declared indices, never freshly minted prefixes, per SQA's naming policy) and which subspaces are symmetric or selected. Each subspace carries a `SubspaceTreatment`, and the treatment chooses which of the three identifications a key applies.

```
 subspace treatment      key function     identifies up to
 ────────────────────────────────────────────────────────────────────────────
 Free      (default)     canon_key        free-index relabelling
 Scaled    (scale)       scaled_key       + permutation symmetry of identical atoms
 Concrete  (evaluate)    concrete_key     nothing (fixed sites σ_1, σ_2, … stay distinct)

 conjugation-aware variants (used to merge conjugate pairs and build the numerics):
   canonical_rep / concrete_rep  →  (key, is_conjugate),  folding ⟨A⟩ with ⟨A†⟩

 Each subspace is independently Free / Scaled / Concrete.
 A hybrid system can be Scaled on the atoms and Concrete on a filter cavity at once.
```

`canon_key` and `scaled_key` are the same engine, `_treatment_key`, called with different treatment maps (all subspaces `Free`, and the selected subspaces `Scaled`). It puts the average into a comparable normal form (drop the summation scope, fix the order of commuting factors, drop the non-equal constraints), relabels each non-`Concrete` subspace to the vocabulary representatives, and reduces any `Scaled` subspace under the permutation symmetry of the identical atoms via `symmetric_min`. A `Concrete` subspace is simply left un-relabelled, so `_treatment_key` already handles mixed systems. `concrete_key` is a separate, simpler label for a fully materialised system: it skips the relabelling entirely so the fixed sites ``\sigma_1, \sigma_2, \ldots`` stay distinct. `canonical_rep` (and its materialised counterpart `concrete_rep`) evaluates the label on both ``\langle A\rangle`` and ``\langle A^\dagger\rangle`` and returns the smaller, together with a flag for which side it came from, so one stored moment stands in for both halves of a conjugate pair.

Because the label is chosen per subspace, **partial reductions compose**: a system can be scaled on one subspace and evaluated on another, in any order, and the labelling stays consistent. `_materialised_key` is the rule that selects the conjugation-aware label once any subspace has been made `Concrete`, and the plain treatment label otherwise.

The `MomentMap{V}` in `canonical.jl` is the single lookup that matches each moment to its stored value, shared by the state-variable registry, `System`, [`get_solution`](@ref), the spectrum routine, and the correlation steady-state lookup. It stores each value under `canonical_rep` and, on a query, reports whether the query sits on the same conjugation side as the stored representative, so a request for ``\langle A^\dagger\rangle`` against a stored ``\langle A\rangle`` returns its conjugate.

## Deriving one moment's equation of motion

`derive(op, sys, ctx)` (`moments.jl`) produces the equation of motion for a single moment.

```
derive(op, sys, ctx)
   │
   ├─ _operator_rhs(direction, op, iH, J, J†, rates)            (op_drift :: QAdd)
   │      Forward :  i[H, op] + Σ_k (r_k/2)(J†[op,J] + [J†,op]J)
   │      Backward: -i[H, op] + adjoint Lindblad + trace term      (retrodiction)
   │
   ├─ expand_completeness            σᵍᵍ fold (no-op without a ground projector)
   ├─ _assume_distinct_atom_indices  ⟨X_i X_j⟩ is the i≠j cumulant; let the diagonal split fire
   │
   ├─ average_and_truncate           cross to moments: average WITH the sum scope, then
   │      │                          cumulant_expansion to sys.order
   │      └─▶ drift :: Num           (the leaves are the moments this one couples to)
   │
   └─ _reduce_ground_in_drift        moment-level σᵍᵍ completeness fold
```

The operator RHS (`operator_drift.jl`) is the Heisenberg (quantum Langevin) equation. `Forward` uses ``+i[H,\cdot]`` plus the Lindblad recycling ``\sum_k r_k\,(J_k^\dagger\,\cdot\,J_k - \tfrac12\{J_k^\dagger J_k,\,\cdot\})``; a matrix rate selects collective decay, a `DoubleIndexedVariable` rate selects the collective cross-jump dissipator (split into an explicit diagonal self-decay plus the off-diagonal term, because completeness cannot recover the diagonal), and otherwise the per-jump term is summed over the jump's free indices. `Backward` flips the commutator sign, swaps ``J\leftrightarrow J^\dagger`` in the recycling, and adds the trace-preserving term for retrodiction.

`average_and_truncate` is where the operator RHS crosses into moments, and the ordering matters: it averages **with** the sum scope attached first, so SQA's diagonal split collapses same-index operator pairs, and only then applies the [`cumulant_expansion`](@ref). Truncating the raw product first would split the operators into separate blocks and the collapse would never fire. Index-dependent coefficients ride inside the sum so the diagonal split substitutes them; scalar coefficients stay outside.

### Cumulant expansion and truncation order

[`cumulant_expansion`](@ref) (`cumulant.jl`) rewrites any average whose order exceeds the cap as the signed sum over partitions of its factors, neglecting the joint cumulant above the cap (see the [theory page](@ref theory)). [`get_order`](@ref) measures order: a number is 0, a single operator 1, a product the number of factors. `order` is normalised to a `Vector{Int}`, one cap per Hilbert subspace; a term acting on several subspaces combines its per-subspace caps through `mix_choice` (default `maximum`).

The subtle part is the **sum-scope metadata round-trip**. A summed average such as ``\Sigma_i\,\langle a^\dagger a\,\sigma_i\rangle`` carries its scope as `SumIndices`/ `SumNonEqual` metadata on that one averaged symbol. When the cumulant expansion factorises it into a product of smaller averages, that symbol is gone, so the scope would be lost. `_stamp_sum_to_first_leaves` re-attaches each bound index onto the first factor that uses it, per additive term, rebuilding the leaf through the canonical `average(SQA.Σ(op, …))`. It must use the canonical form and not a bare `setmetadata`, because a `setmetadata` leaf is `isequal` to the un-summed leaf and the two would never cancel in a later subtraction.

## Closing the hierarchy

`seed` (`graph.jl`) derives the equations for the user's requested operators, the first entries in the graph. `closure!` then repeatedly takes one moment, looks at the moments on the right-hand side of its equation (the *leaves* of the symbolic drift), and derives an equation for any not yet present, until no new moments appear.

A moment and its conjugate are one physical unknown, and the closure exploits that *partially*. A moment whose conjugate partner is already present is covered automatically (no new equation). The `get_adjoints` flag controls only the genuinely-new case. With `get_adjoints=true` (the default for `meanfield`/`complete`) a new moment's conjugate is also added as its own moment; with `get_adjoints=false` (used by [`CorrelationFunction`](@ref)) only one representative is kept and the partner is recovered by complex conjugation when the numerical system is built. The consequences are spelled out in the worked example below and in invariant 1.

The `max_iter` guard is a runaway backstop, not a truncation limiter: hitting it raises an error rather than silently returning a non-closed system, which the numerical build would otherwise mask.

[`complete!`](@ref) is `closure!` exposed publicly: it rebuilds a graph from an existing equation set (`_graph_from_eqs`) and extends it. A `filter_func` drops unwanted moments (for example phase-invariant or ancilla terms) by zeroing them on every right-hand side. [`find_missing`](@ref) reports the still-needed moments and is representation-invariant: it folds conjugate pairs onto one key, so the answer is the same whether the system was closed with `get_adjoints` true or false.

## Reducing the system: `scale` and `evaluate`

[`scale`](@ref) and [`evaluate`](@ref) transform the graph by relabelling the moments under a new treatment. Both build the graph from the **stored** drifts (`_graph_from_stored`) rather than re-deriving, because a `complete(...; filter_func)` records filter substitutions in the stored equations that a re-derivation would discard.

- [`scale`](@ref) (`scaling.jl`) calls `quotient`, marking the selected symmetric subspaces `Scaled` and merging permutation-equivalent moments. Each surviving drift leaf ``\langle X\rangle`` becomes ``\text{prefactor}\cdot\langle X_{\text{rep}}\rangle``, where the prefactor is ``\prod_b (\text{range}_b - \#\text{ne}_b)`` over the bound indices on the scaled subspaces, and every per-atom `IndexedVariable` coupling is flattened to its scalar (all atoms share the coupling under the symmetry). The `h::Vector{Int}` keyword restricts the quotient to specific subspaces; empty means all symmetric subspaces.

- [`evaluate`](@ref) (`evaluate.jl`) calls `specialize`, unrolling each targeted symbolic range to its concrete integer size (from `limits`) and pinning indices to `Concrete` sites ``1\ldots M`` with distinct-site (``i\neq j``) semantics. A final `_arrayize_indexed_params!` rewrites any leftover per-site coupling (a callable `Sym` like `g(i_2_3)`) into `getindex` on a freshly minted Symbolics array parameter, which is what MTK can scalarise and bind; [`parameter_map`](@ref) later matches user values to that array by name.

Because both are keyed by treatment rather than a hardcoded representative, a system can be partially scaled and partially evaluated, in any order.

## Building the numerical system

`mtk.jl` bridges the moment layer to ModelingToolkitBase. `_state_registry` assigns each state a time-dependent variable `u(t)`, named deterministically by `serialize` (``\langle a^\dagger a\rangle`` becomes `avg_a_dag_a`), and builds a `MomentMap` keyed by the conjugate representative. Building the `System` then looks up every moment on a right-hand side in that map.

```
 right-hand side  (symbolic expression; each ⟨...⟩ average is a leaf)
        │  mapleaves(resolve, rhs)
        ▼
 for each leaf ⟨op⟩:
        ├─ canonical_rep(op) ─▶ (rep, side)
        ▼
   MomentMap.by_rep[rep] = (var, rep_side)
        ├─ side == rep_side ?  ─ yes ─▶  var          e.g. avg_a(t)
        └─                       no  ─▶  conj(var)     e.g. conj(avg_a(t))
```

The conjugate is emitted as `SymbolicUtils.term(conj, var; type=Number)` rather than a plain `conj(var)`, because `conj` folds to the identity on a `Real` symtype and would erase the distinction; the explicit `term` survives `mtkcompile`. The left-hand side becomes `Differential(iv)(u(t))` only here; the equation struct itself always stores the raw `Average`.

A `NoiseMeanfieldEquations` becomes a stochastic differential equation: the single Brownian term is the aggregated noise drift, with sign ``+1`` for `Forward` and ``-1`` for `Backward`, chosen by the evolution direction.

The remaining bridges: [`initial_values`](@ref) computes a `u0` from a numeric `Ket` or density operator (via SQA's `numeric_average`), [`parameter_map`](@ref) expands indexed/array parameters to the compiled system's shapes, and [`get_solution`](@ref) reverses the registry to evaluate any operator-average trajectory, including products that never appeared as a stored state. `System` uses the `MomentMap` alone; `get_solution` adds a secondary `by_canon` fallback so a query posed in a different but equivalent symbolic form still resolves.

## Correlation and spectrum

[`CorrelationFunction`](@ref) (`correlation.jl`) is another step that builds its own hierarchy. For the two-time correlation ``\langle A(t+\tau)B(t)\rangle`` it re-embeds ``B`` on a fresh ancilla subspace (`aon = max(acts_on) + 1`), derives the equations for ``A B_{\text{ancilla}}`` in the delay ``\tau``, and closes only on the averages that touch the ancilla (a closure filter), with `get_adjoints=false`. The non-ancilla averages are steady-state coefficients. The correlation inherits the parent system's treatments so its leaves resolve in the same keying. `Spectrum` is the Fourier transform of the correlation.

## A worked example: the driven Kerr cavity

The smallest model that exercises the whole pipeline, including a real cumulant truncation, is a single Kerr (anharmonic) cavity. A *linear* cavity is Gaussian and never truncates at order 2, so the Kerr term ``U\,a^\dagger a^\dagger a a`` is what makes the machinery visible.

```@example devdocs
using QuantumCumulants
h = FockSpace(:cavity)
@qnumbers a::Destroy(h)
@variables Δ U η κ
H = Δ*a'*a + U*a'*a'*a*a + η*(a + a')
nothing # hide
```

**Derivation.** Without truncation, `meanfield` produces the exact Heisenberg equation for ``\langle a\rangle``. The Kerr term injects an order-3 moment ``\langle a^\dagger a a\rangle``:

```@example devdocs
meanfield([a], H, [a]; rates=[κ])
```

**Truncation.** Passing `order=2` applies the cumulant expansion: the order-3 moment is factorised into products of moments of order ``\le 2``. This is the truncation firing.

```@example devdocs
eqs = meanfield([a], H, [a]; rates=[κ], order=2)
```

**Closure.** The right-hand side now references ``\langle a^\dagger a\rangle``, ``\langle a a\rangle`` and ``\langle a^\dagger\rangle``. The conjugate ``\langle a^\dagger\rangle`` is covered (its partner ``\langle a\rangle`` is already a state), but because the default is `get_adjoints=true`, the new moment ``\langle a a\rangle`` has its conjugate ``\langle a^\dagger a^\dagger\rangle`` tracked as its own state too, giving four states:

```@example devdocs
eqs_c = complete(eqs)
eqs_c.states
```

Asking for the minimal set with `get_adjoints=false` keeps three, recovering ``\langle a^\dagger a^\dagger\rangle`` purely by conjugation:

```@example devdocs
complete(eqs; get_adjoints=false).states
```

**Numerics.** Building the MTK system names one variable per state and looks up every moment through the `MomentMap`. The ``\langle a^\dagger\rangle`` terms resolve to `conj(avg_a(t))`:

```@example devdocs
using ModelingToolkitBase
sys = System(eqs_c; name=:kerr)
unknowns(sys)
```

```@example devdocs
equations(sys)[1]   # ⟨a⟩: note conj(avg_a(t)) standing in for ⟨a†⟩
```

Note the honest subtlety from the default `get_adjoints=true`: `avg_a_dag_a_dag(t)` is its own integrated variable, but its right-hand side resolves to `conj(avg_a_a(t))` and nothing else reads it, so it is a redundant shadow. `mtkcompile` does **not** detect the conjugate alias and keeps it, so the default integrates one redundant ODE. Use `get_adjoints=false` when the minimal state set matters.

### What changes for an indexed system

For a permutation-symmetric many-body system built from [`Index`](@ref) and [`Σ`](@ref) (see [Indexed and scaled systems](@ref)), the same graph carries indexed moments. A single-atom moment ``\langle\sigma_i\rangle`` keys to one `Free` representative regardless of the bound index. [`scale`](@ref) then marks the atom subspace `Scaled`, folds the permutation-equivalent moments onto one representative, and multiplies each by the ``(\text{range} - \#\text{ne})`` prefactor that counts how many atoms the fold stands for. [`evaluate`](@ref) instead pins the atoms to `Concrete` sites ``\sigma_1\ldots\sigma_M`` and arrayizes the per-site couplings. Treatments are per-subspace and independent, so a filter-cavity subspace can stay `Concrete` while the atoms are `Scaled` in the same system.

## Design invariants and pitfalls

These are the moment-layer rules a contributor must not break, each with its failure mode.

12. **A conjugate pair is one physical unknown, but the default does not minimise on it.** During `closure!` a moment whose conjugate is already present is covered. `get_adjoints=true` (default) tracks the conjugate of a genuinely-new moment as a second state; the numerical build keys states by `canonical_rep`/`concrete_rep`, so every occurrence of the pair on a right-hand side resolves to the first representative (or its `conj`), leaving the second variable a redundant shadow that `mtkcompile` does not eliminate. Use `get_adjoints=false` for the minimal set.

12. **`closure!` keys with `canon_key`, never `scaled_key`.** Symmetric reduction during closure would collapse the unscaled moment count. The reduction steps relabel with the treatment-matched `_materialised_key`; calling bare `canon_key` on an already-`Concrete` subspace would relabel and merge sites that must stay distinct.

12. **Build `scale`/`evaluate` from stored drifts, never re-derive.** `_graph_from_stored` reuses the recorded right-hand sides; a re-derivation discards filter substitutions a prior `complete(...; filter_func)` applied.

12. **Sum-scope metadata must round-trip.** Rebuild summed leaves through the canonical `average(SQA.Σ(...))`, not `setmetadata`, which is `isequal` to the un-summed leaf so the terms never cancel.

12. **Use `Symbolics.IM`, not Julia's `im`, in symbolic right-hand sides.** A `complex(0, …)` literal does not unify with the factored symbolic form; `_im_form`/`_rebuild` convert it.

12. **Mint indices from the user's vocabulary** (`reps[1](k)`), never invent fresh prefixes (SQA naming policy). `scale` and `evaluate` must mint the same slots so they agree.

12. **Route index operations through `SQA.change_index`/`Σ`/`*`.** Hand-rolling a `QAdd` and sorting it yourself reintroduces the bugs the algebra exists to prevent (`` igma_{ee}^2 \ne igma_{ee}``, `` igma_2 igma_1 \ne igma_1 igma_2``).

12. **Average with the sum scope before truncating**, so the diagonal split fires before factorisation.

12. **The `iszero` coefficient drop in `average_and_truncate` uses `Base.iszero` deliberately**; SQA's structural check misses uncombined zero forms like ``\lambda/2 - (1/2)\lambda``.

12. **`max_iter` in `closure!` errors, it does not truncate.** It is a runaway backstop only.

12. **The ground projector `` igma^{gg}`` stays an atom** except where the `` igma^{gg} = 1 - \Sigma\, igma^{kk}`` fold is the explicit goal (`expand_completeness`/`_reduce_ground_in_drift` in `derive`).

12. **`NodeData.op_noise` is always `nothing`** (deferred); the operator-level noise column is rebuilt in `assemble_equations` when a noise struct is assembled.
