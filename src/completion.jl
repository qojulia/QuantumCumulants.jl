"""
    find_missing(eqs::AbstractMeanFieldEquations; filter_func=nothing,
                 get_adjoints=true)

Return the unique averages appearing on RHSes of `eqs.equations` whose
identity is not already present in `eqs.states`.

An average ⟨X⟩ and its conjugate ⟨X†⟩ are treated as the same state for
the purpose of "already covered" detection — if either is a state, neither
counts as missing. When `get_adjoints=false`, only one of each
conjugate-pair is emitted into `missing_states` (the conjugate is marked
seen so it isn't returned again on a later RHS scan).
"""
function find_missing(eqs::AbstractMeanFieldEquations; filter_func = nothing,
                      get_adjoints::Bool = true)
    seen = Set{SymbolicUtils.BasicSymbolic}()
    for s in eqs.states
        push!(seen, s)
        push!(seen, _avg_conj(s))
    end
    missing_states = SymbolicUtils.BasicSymbolic[]
    for eq in eqs.equations
        _collect_missing!(missing_states, seen, eq.rhs, get_adjoints)
    end
    if filter_func !== nothing
        filter!(filter_func, missing_states)
    end
    return missing_states
end

function _collect_missing!(missing_states, seen, x, get_adjoints::Bool = true)
    if x isa SymbolicUtils.BasicSymbolic
        if _is_leaf_average(x)
            x in seen && return
            push!(seen, x)
            get_adjoints || push!(seen, _avg_conj(x))
            push!(missing_states, x)
            return
        end
        if SymbolicUtils.iscall(x)
            for a in SymbolicUtils.arguments(x)
                _collect_missing!(missing_states, seen, a, get_adjoints)
            end
        end
    end
    return
end

# Conjugate of a leaf average: ⟨op⟩ ↦ ⟨op'⟩. We build it via SQA's `average`
# constructor rather than calling `conj` so the result lands in the same
# canonical form (sym_average head, same hashing) as states pushed by
# `meanfield`.
function _avg_conj(x::SymbolicUtils.BasicSymbolic)
    _is_leaf_average(x) || return x
    op = SecondQuantizedAlgebra.undo_average(x)
    return average(adjoint(op))
end

# `SecondQuantizedAlgebra.is_average` returns true for the whole AvgSym
# expression tree (e.g. ⟨a⟩*⟨σ₂₂⟩, which is a product of averages, still has
# symtype === AvgSym). For closure detection we want only *leaf* averages —
# the ones produced by `average(op)` directly, where the head is `sym_average`.
function _is_leaf_average(x::SymbolicUtils.BasicSymbolic)
    SecondQuantizedAlgebra.is_average(x) || return false
    return SymbolicUtils.operation(x) === SecondQuantizedAlgebra.sym_average
end

"""
    complete!(eqs::AbstractMeanFieldEquations; max_iter=200, simplify=true,
              filter_func=nothing, mix_choice=maximum)

Iteratively derive equations for missing averages until the system is closed
(or `max_iter` is reached). Mutates `eqs` in place.

The order used for cumulant expansion is the one already stored in
`eqs.order` (set at `meanfield(...; order=...)` time, or via
`cumulant_expansion(eqs, order)`). The non-mutating [`complete`](@ref)
variant additionally accepts an `order` keyword.
"""
function complete!(eqs::MeanFieldEquations; max_iter::Int = 200,
                   simplify::Bool = true, filter_func = nothing,
                   mix_choice = maximum, get_adjoints::Bool = true)
    for _ in 1:max_iter
        missing_states = find_missing(eqs; filter_func, get_adjoints)
        isempty(missing_states) && return eqs
        new_ops = QField[SecondQuantizedAlgebra.undo_average(m) for m in missing_states]
        new_eqs = _derive_for(eqs, new_ops; simplify, mix_choice)
        if filter_func !== nothing
            _filter_rhs!(new_eqs, filter_func)
        end
        _append!(eqs, new_eqs)
    end
    error("complete!: did not close within $max_iter iterations")
end

function _filter_rhs!(eqs::MeanFieldEquations, filter_func)
    for (i, eq) in enumerate(eqs.equations)
        new_rhs = _filter_expr(eq.rhs, filter_func)
        eqs.equations[i] = eq.lhs ~ new_rhs
    end
    return eqs
end

function _filter_expr(x, filter_func)
    if x isa SymbolicUtils.BasicSymbolic
        if SecondQuantizedAlgebra.is_average(x)
            return filter_func(x) ? x : 0
        end
        if SymbolicUtils.iscall(x) && _has_average(x)
            op = SymbolicUtils.operation(x)
            args = SymbolicUtils.arguments(x)
            new_args = Any[_filter_expr(a, filter_func) for a in args]
            return op(new_args...)
        end
    end
    return x
end

"""
    complete(eqs::AbstractMeanFieldEquations; order=nothing, kw...)

Non-mutating variant. Returns a new `MeanFieldEquations`. If `order` is
given, the returned system is first cumulant-expanded to that order before
the completion loop runs.
"""
function complete(eqs::MeanFieldEquations; order = nothing,
                  mix_choice = maximum, simplify::Bool = true, kw...)
    eqs_copy = order === nothing ?
        _copy(eqs) :
        cumulant_expansion(_copy(eqs), order; simplify, mix_choice)
    return complete!(eqs_copy; simplify, mix_choice, kw...)
end

function _derive_for(eqs::MeanFieldEquations, new_ops;
                     simplify::Bool, mix_choice = maximum)
    return _meanfield_forward(new_ops, eqs.hamiltonian, eqs.jumps,
                              eqs.jumps_dagger, eqs.rates, eqs.order,
                              simplify, mix_choice, eqs.iv)
end

function _append!(a::MeanFieldEquations, b::MeanFieldEquations)
    append!(a.equations, b.equations)
    append!(a.operator_equations, b.operator_equations)
    append!(a.states, b.states)
    append!(a.operators, b.operators)
    return a
end

function _copy(eqs::MeanFieldEquations)
    return MeanFieldEquations(copy(eqs.equations),
                              copy(eqs.operator_equations),
                              copy(eqs.states),
                              copy(eqs.operators),
                              eqs.hamiltonian,
                              copy(eqs.jumps),
                              copy(eqs.jumps_dagger),
                              copy(eqs.rates),
                              eqs.iv, eqs.order)
end
