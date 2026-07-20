"""
Flatten one level of nesting in `jumps`/`jumps_dagger`: collective-decay sources
may arrive as a vector of mode vectors, and this leaves each source a single
`QField`.
"""
_flatten_jumps(js::AbstractVector{<:QField}) = js
function _flatten_jumps(js::AbstractVector)
    isempty(js) && return QField[]
    eltype(js) <: AbstractVector || return js
    out = QField[]
    for jk in js
        append!(out, jk)
    end
    return out
end

_make_iv() = MTK.t_nounits

"""
    meanfield(ops, H, J=QField[]; kwargs...)
    meanfield(op, H, J=QField[]; kwargs...)

Equations of motion for the averages of `ops` under the Hamiltonian `H` and the
collapse operators `J`. The result is the Heisenberg (quantum Langevin) equation with
noise neglected: each jump adds a Lindblad term
``\\sum_i r_i \\left(J_i^\\dagger O J_i - \\tfrac{1}{2}\\{J_i^\\dagger J_i, O\\}\\right)``.
Returns a [`MeanfieldEquations`](@ref), or a [`NoiseMeanfieldEquations`](@ref) when
`efficiencies` is supplied (measurement backaction).

The returned right-hand sides are left unsimplified; call [`simplify!`](@ref) (or
`SymbolicUtils.simplify`) on the equations you want to inspect.

# Arguments
* `ops`: the operator(s) whose equations of motion are derived.
* `H`: the operator defining the system Hamiltonian.
* `J`: the collapse (jump) operators of the system; defaults to none.

# Keyword arguments
* `Jdagger=adjoint.(J)`: the adjoints of the jump operators.
* `rates=ones(length(J))`: decay rates for `J`. A square `Matrix` instead selects
  collective dissipation across the jump vector.
* `efficiencies=nothing`: detector efficiencies per jump. When given, the result is a
  [`NoiseMeanfieldEquations`](@ref) carrying a measurement-backaction noise drift.
* `direction=Forward()`: [`Forward`](@ref) for ordinary evolution, [`Backward`](@ref)
  for retrodiction.
* `order=nothing`: if set, a [`cumulant_expansion`](@ref) to this order is applied. An
  `Int` caps every subspace; a `Vector{Int}` caps each Hilbert subspace separately.
* `mix_choice=maximum`: for a `Vector` `order`, how to combine the per-subspace caps on a
  term that acts on several subspaces.
* `iv=ModelingToolkitBase.t_nounits`: the independent (time) variable of the system.
"""
function meanfield(
        ops::AbstractVector,
        H::QField,
        J::AbstractVector = QField[];
        Jdagger = nothing,
        rates::Union{AbstractVector, AbstractMatrix} = ones(length(J)),
        efficiencies = nothing,
        direction::EvolutionDirection = Forward(),
        order = nothing,
        mix_choice::Function = maximum,
        iv::Symbolics.Num = _make_iv(),
    )
    # Collective dissipation: wrap a flat `(J, rates::Matrix)` into the
    # nested-cluster layout the Lindblad recycling expects.
    if rates isa AbstractMatrix
        size(rates, 1) == size(rates, 2) == length(J) || throw(
            ArgumentError(
                "rates::Matrix must be square with side length(J) for collective dissipation",
            )
        )
        J = [collect(J)]
        Jdagger = Jdagger === nothing ? nothing : [collect(Jdagger)]
        rates = [rates]
    end
    Jdagger === nothing && (Jdagger = _default_jdagger(J))
    Jn, Jdn = _normalize_jumps(J, Jdagger)
    rn = _normalize_rates(rates, length(Jn))
    order_vec = _normalize_order(
        order, (; operators = ops, hamiltonian = H, jumps = Jn, jumps_dagger = Jdn)
    )
    eff = efficiencies === nothing ? nothing : _normalize_rates(efficiencies, length(Jn))

    ctx = build_ctx(ops, H, Jn, Jdn)
    sys = SystemSpec(H, Jn, Jdn, rn, eff, iv, order_vec, mix_choice, direction)
    g = seed(QAdd[op * 1 for op in ops], sys, ctx)
    return assemble_equations(g)
end

meanfield(op::QField, H::QField, args...; kw...) = meanfield([op], H, args...; kw...)

function _normalize_jumps(J, Jdagger)
    isempty(J) && return QField[], QField[]
    return collect(J), collect(Jdagger)
end

# Default `Jdagger`. For the flat form `J::Vector{QField}` this is `adjoint.(J)`;
# for collective decay `J::Vector{Vector{QField}}` adjoint each inner mode-vector
# entry-wise (the outer-level `adjoint.` would transpose the inner vector).
function _default_jdagger(J)
    isempty(J) && return QField[]
    eltype(J) <: AbstractVector && return [adjoint.(jk) for jk in J]
    return adjoint.(J)
end

function _normalize_rates(rates, n::Int)
    isempty(rates) && n == 0 && return Symbolics.Num[]
    return collect(rates)
end

"""
    simplify!(eqs::AbstractMeanfieldEquations; kwargs...)
    simplify(eqs::AbstractMeanfieldEquations; kwargs...)

Run `SymbolicUtils.simplify` on every RHS in `eqs` (and on the noise drift RHSs
of a `NoiseMeanfieldEquations`). `simplify!` mutates; `simplify` returns a fresh
struct. The derivation pipeline leaves expressions raw so the cost of
`SymbolicUtils.simplify` is paid only when explicitly requested.
"""
function simplify!(eqs::AbstractMeanfieldEquations; kwargs...)
    eqs.graph = map_drifts(eqs.graph, (_, d) -> SymbolicUtils.simplify(d; kwargs...))
    resync!(eqs; moments_unchanged = true)
    return eqs
end

SymbolicUtils.simplify(eqs::AbstractMeanfieldEquations; kwargs...) =
    simplify!(_copy(eqs); kwargs...)

# ---- measurement record: rewrite dW-parametrised SDE to dY ----

"""
    translate_W_to_Y(eqs::NoiseMeanfieldEquations; mix_choice=maximum)

Rewrite an SDE whose noise drift is parametrised by the underlying Wiener
process `dW` into one parametrised by the measurement record `dY` instead.
The substitution `dW = dY - sqrt(2η)·⟨J + J†⟩·dt` adds a deterministic
correction to the drift; the noise drift itself is unchanged.

For each equation the RHS is augmented by
`cumulant_expansion(-_dY_dS_extra_term(lhs_op, J, Jdagger, rates .* efficiencies))`,
cumulant-expanded to `eqs.order` when set. Returns a fresh
`NoiseMeanfieldEquations` of the same direction. The augmented RHS is left
unsimplified; apply `simplify` yourself for a canonical form.
"""
function translate_W_to_Y(
        eqs::NoiseMeanfieldEquations;
        mix_choice::Function = maximum,
    )
    J, Jd = collect(eqs.jumps), collect(eqs.jumps_dagger)
    rates_eff = collect(eqs.rates) .* collect(eqs.efficiencies)
    order = eqs.order
    augment = function (k, d)   # `k` is the node key, i.e. the LHS operator
        term = -_dY_dS_extra_term(k, J, Jd, rates_eff)
        order !== nothing && (term = cumulant_expansion(term, order; mix_choice))
        return d + term
    end
    g = map_drifts(eqs.graph, augment; noise = false)
    return assemble_equations(g)
end

# ---- user-supplied RHS rewrite ----

"""
    modify_equations(eqs::AbstractMeanfieldEquations, f::Function)

Return a copy of `eqs` whose RHS for each equation has been rewritten by `f`.
The function is called as `f(lhs_op, rhs)`, receiving the *operator* form of
the LHS (`undo_average(eq.lhs)`) and the symbolic RHS, and returning a new RHS.

```julia
f(lhs, rhs) = rhs + cumulant_expansion(average(commutator(1im * Hadd, lhs)), 2)
eqs_mod = modify_equations(eqs, f)
```

See also [`modify_equations!`](@ref).
"""
modify_equations(eqs::AbstractMeanfieldEquations, f::Function) =
    modify_equations!(_copy(eqs), f)

"""
    modify_equations!(eqs::AbstractMeanfieldEquations, f::Function)

In-place version of [`modify_equations`](@ref). Rewrites each drift in `eqs.graph` with
`f(lhs_op, rhs)`, where `lhs_op` is the node key (the operator form of the LHS).
"""
function modify_equations!(eqs::AbstractMeanfieldEquations, f::Function)
    eqs.graph = map_drifts(eqs.graph, (k, d) -> f(k, d); noise = false)
    resync!(eqs; moments_unchanged = true)
    return eqs
end

# ---- substitute into the RHS ----

"""
    substitute(eqs::AbstractMeanfieldEquations, dict::AbstractDict)
    substitute!(eqs::AbstractMeanfieldEquations, dict::AbstractDict)

Substitute `dict` into every RHS of `eqs` (and into the noise drift RHSs of a
`NoiseMeanfieldEquations`). Keys of `dict` are symbolic expressions to replace: parameters,
or averages known to vanish (e.g. by symmetry), so that zeroing them closes the hierarchy
onto the remaining moments. `substitute!` mutates `eqs`; `substitute` returns a fresh struct.

Only the drifts are rewritten; the moment set (the LHS states and operators) is unchanged. To
also drop the substituted moments as unknowns, complete onto the reduced subset (e.g. with a
`filter_func`) or rebuild `meanfield` on it.

```julia
# close a phase-invariant model onto ⟨σ⁺ᵢσ⁻ⱼ⟩ by asserting the single-operator averages vanish
eqs = meanfield(ops, H, J; rates, order = 2)
eqs = substitute(eqs, Dict(average(σ(2, 1, i)) => 0 for i in 1:N))
```
"""
function substitute!(eqs::AbstractMeanfieldEquations, dict::AbstractDict)
    # `_subtree_substitute` rewrites only the matched subtrees, leaving the shape
    # metadata of the untouched averages intact. `Symbolics.substitute` rebuilds
    # the whole expression and strips the scalar shape off every `⟨…⟩`, which then
    # no longer matches the (scalar) LHS averages and breaks downstream simplify.
    sub = Dict(SymbolicUtils.unwrap(k) => SymbolicUtils.unwrap(v) for (k, v) in dict)
    eqs.graph = map_drifts(eqs.graph, (_, d) -> _subtree_substitute(SymbolicUtils.unwrap(d), sub))
    resync!(eqs; moments_unchanged = true)
    return eqs
end

SymbolicUtils.substitute(eqs::AbstractMeanfieldEquations, dict::AbstractDict) =
    substitute!(_copy(eqs), dict)
