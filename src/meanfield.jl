# Public entry point. `meanfield` normalizes the jump/rate inputs, seeds a moment
# graph, and lowers it to the array container. Deterministic and noise systems
# share one path: `efficiencies` flips `SystemSpec.efficiencies`, which makes
# `derive` compute the noise drift and `lower_to_eqs` emit the noise columns. No
# bespoke per-path assembly.

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
    return lower_to_eqs(g)
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
