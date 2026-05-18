"""
    scale!(eqs::MeanFieldEquations)

Collapse states equivalent under permutation symmetry of their indices and drop
redundant equations.

For now this is a minimal pass-through: states are inspected and equations
that are identical under index renaming are collapsed. Most production
analyses already get index canonicalization from `complete!` via SQA's
eager normal-ordering, so this function is conservative.
"""
function scale!(eqs::MeanFieldEquations)
    n = length(eqs.states)
    n <= 1 && return eqs
    canon = IdDict{SymbolicUtils.BasicSymbolic, SymbolicUtils.BasicSymbolic}()
    seen_keys = Dict{Any, Int}()
    keep_idx = Int[]
    for (k, s) in enumerate(eqs.states)
        key = _state_canonical_key(s)
        if haskey(seen_keys, key)
            canon[s] = eqs.states[seen_keys[key]]
        else
            seen_keys[key] = k
            canon[s] = s
            push!(keep_idx, k)
        end
    end
    sub = Dict(s => canon[s] for s in eqs.states if canon[s] !== s)
    new_eqs = Symbolics.Equation[]
    new_states = SymbolicUtils.BasicSymbolic[]
    new_ops = QAdd[]
    new_op_eqs = Symbolics.Equation[]
    for i in keep_idx
        new_rhs = isempty(sub) ? eqs.equations[i].rhs :
            SymbolicUtils.substitute(eqs.equations[i].rhs, sub)
        push!(new_eqs, eqs.states[i] ~ new_rhs)
        push!(new_states, eqs.states[i])
        push!(new_ops, eqs.operators[i])
        push!(new_op_eqs, eqs.operator_equations[i])
    end
    empty!(eqs.equations);          append!(eqs.equations, new_eqs)
    empty!(eqs.operator_equations); append!(eqs.operator_equations, new_op_eqs)
    empty!(eqs.states);             append!(eqs.states, new_states)
    empty!(eqs.operators);          append!(eqs.operators, new_ops)
    return eqs
end

"""
    scale(eqs::MeanFieldEquations)

Non-mutating variant.
"""
scale(eqs::MeanFieldEquations) = scale!(_copy(eqs))

function _state_canonical_key(avg::SymbolicUtils.BasicSymbolic)
    op = SQA.undo_average(avg)
    return _qfield_key(op)
end

function _qfield_key(op::QAdd)
    out = Any[]
    for (term, coeff) in op.arguments
        shape_keys = sort!([_op_shape_key(o) for o in term.ops])
        push!(out, (shape_keys, length(term.ne), coeff))
    end
    sort!(out; by = repr)
    return out
end
_qfield_key(op::SQA.QSym) = _op_shape_key(op)

_op_shape_key(op::SQA.QSym) =
    (string(typeof(op)), op.name, op.space_index)
