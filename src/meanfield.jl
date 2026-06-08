# Public entry point. `meanfield` normalizes the jump/rate inputs, seeds a moment
# graph, and lowers it to the array container. Deterministic and noise systems
# share one path: `efficiencies` flips `SystemSpec.efficiencies`, which makes
# `derive` compute the noise drift and `lower_to_eqs` emit the noise columns. No
# bespoke per-path assembly.

_make_iv() = MTK.t_nounits

"""
    meanfield(ops, H, J=QField[]; Jdagger=adjoint.(J), rates=ones(length(J)),
              efficiencies=nothing, direction=Forward(),
              order=nothing, mix_choice=maximum, iv=ModelingToolkitBase.t_nounits)

Equations of motion for the averages of `ops` under Hamiltonian `H` and collapse
operators `J` (with `rates`). Returns a `MeanFieldEquations`, or a
`NoiseMeanFieldEquations` when `efficiencies` is given.

The returned RHS is left unsimplified. Apply `SymbolicUtils.simplify` (e.g.
`Symbolics.simplify(eq.rhs; expand=true)`) to the equations you want to inspect.
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
