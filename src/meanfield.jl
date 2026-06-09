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
Returns a [`MeanFieldEquations`](@ref), or a [`NoiseMeanFieldEquations`](@ref) when
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
  [`NoiseMeanFieldEquations`](@ref) carrying a measurement-backaction noise drift.
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
    order_vec = _normalize_order(order, (; hamiltonian = H))
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
    simplify!(eqs::AbstractMeanFieldEquations; kwargs...)
    simplify(eqs::AbstractMeanFieldEquations; kwargs...)

Run `SymbolicUtils.simplify` on every RHS in `eqs` (and on the noise drift RHSs
of a `NoiseMeanFieldEquations`). `simplify!` mutates; `simplify` returns a fresh
struct. The derivation pipeline leaves expressions raw so the cost of
`SymbolicUtils.simplify` is paid only when explicitly requested.
"""
function simplify! end

function simplify!(eqs::MeanFieldEquations; kwargs...)
    for (i, eq) in enumerate(eqs.equations)
        eqs.equations[i] = eq.lhs ~ SymbolicUtils.simplify(eq.rhs; kwargs...)
    end
    return eqs
end

function simplify!(eqs::NoiseMeanFieldEquations; kwargs...)
    for (i, eq) in enumerate(eqs.equations)
        eqs.equations[i] = eq.lhs ~ SymbolicUtils.simplify(eq.rhs; kwargs...)
    end
    for (i, eq) in enumerate(eqs.noise_equations)
        eqs.noise_equations[i] = eq.lhs ~ SymbolicUtils.simplify(eq.rhs; kwargs...)
    end
    return eqs
end

SymbolicUtils.simplify(eqs::AbstractMeanFieldEquations; kwargs...) =
    simplify!(_copy(eqs); kwargs...)

# ---- measurement record: rewrite dW-parametrised SDE to dY ----

"""
    translate_W_to_Y(eqs::NoiseMeanFieldEquations; mix_choice=maximum)

Rewrite an SDE whose noise drift is parametrised by the underlying Wiener
process `dW` into one parametrised by the measurement record `dY` instead.
The substitution `dW = dY - sqrt(2η)·⟨J + J†⟩·dt` adds a deterministic
correction to the drift; the noise drift itself is unchanged.

For each equation the RHS is augmented by
`cumulant_expansion(-_dY_dS_extra_term(lhs_op, J, Jdagger, rates .* efficiencies))`,
cumulant-expanded to `eqs.order` when set. Returns a fresh
`NoiseMeanFieldEquations` of the same direction. The augmented RHS is left
unsimplified; apply `simplify` yourself for a canonical form.
"""
function translate_W_to_Y(
        eqs::NoiseMeanFieldEquations;
        mix_choice::Function = maximum,
    )
    out = _copy(eqs)
    J, Jd = out.jumps, out.jumps_dagger
    rates_eff = out.rates .* out.efficiencies
    for i in eachindex(out.equations)
        eq_i = out.equations[i]
        lhs_op = undo_average(eq_i.lhs)
        term = -_dY_dS_extra_term(lhs_op, J, Jd, rates_eff)
        out.order !== nothing && (term = cumulant_expansion(term, out.order; mix_choice))
        out.equations[i] = eq_i.lhs ~ eq_i.rhs + term
    end
    return out
end

# ---- user-supplied RHS rewrite ----

"""
    modify_equations(eqs::AbstractMeanFieldEquations, f::Function)

Return a copy of `eqs` whose RHS for each equation has been rewritten by `f`.
The function is called as `f(lhs_op, rhs)`, receiving the *operator* form of
the LHS (`undo_average(eq.lhs)`) and the symbolic RHS, and returning a new RHS.

```julia
f(lhs, rhs) = rhs + cumulant_expansion(average(commutator(1im * Hadd, lhs)), 2)
eqs_mod = modify_equations(eqs, f)
```

See also [`modify_equations!`](@ref).
"""
modify_equations(eqs::AbstractMeanFieldEquations, f::Function) =
    modify_equations!(_copy(eqs), f)

"""
    modify_equations!(eqs::AbstractMeanFieldEquations, f::Function)

In-place version of [`modify_equations`](@ref). Walks `eqs.equations` and replaces each
RHS with `f(undo_average(lhs), rhs)`.
"""
function modify_equations!(eqs::AbstractMeanFieldEquations, f::Function)
    for i in eachindex(eqs.equations)
        lhs = eqs.equations[i].lhs
        rhs = eqs.equations[i].rhs
        eqs.equations[i] = lhs ~ f(undo_average(lhs), rhs)
    end
    return eqs
end
