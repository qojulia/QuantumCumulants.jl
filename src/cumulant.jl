using Combinatorics: partitions

"""
    get_order(expr)

Highest moment order appearing in `expr`. Numbers have order 0; bare operators
have order 1; products have order equal to the number of operator factors.
"""
get_order(::Number) = 0
get_order(::SQA.QSym) = 1
function get_order(q::QAdd)
    isempty(q.arguments) && return 0
    return maximum(length(t.ops) for (t, _) in q.arguments)
end
function get_order(x::SymbolicUtils.BasicSymbolic)
    # Treat only *leaf* averages ⟨op⟩ as moments. A product of averages
    # ⟨a⟩·⟨b⟩ still has symtype === AvgSym (so `is_average` returns true)
    # but its "order" for cumulant-truncation purposes is the max order of
    # its constituent leaf averages, not the order of `undo_average`'s
    # operator-product re-assembly.
    if SQA.is_average(x) &&
            SymbolicUtils.iscall(x) &&
            SymbolicUtils.operation(x) === SQA.sym_average
        return get_order(SQA.undo_average(x))
    end
    if SymbolicUtils.iscall(x)
        return maximum(get_order, SymbolicUtils.arguments(x); init = 0)
    end
    return 0
end
get_order(x::Symbolics.Num) = get_order(SymbolicUtils.unwrap(x))

# Build the operator product for a block of QSym operators.
_prod_ops(block::AbstractVector) = isempty(block) ? 1 : reduce(*, block)

"""
    cumulant(op, n=get_order(op))

The `n`-th joint cumulant of `op` (signed sum over partitions of length n).
The result is returned in its raw, unsimplified form. Apply
`SymbolicUtils.simplify` yourself if you want a canonical representation.
"""
cumulant(op::SQA.QSym, n::Int = 1) = n == 1 ? average(op) : 0

function cumulant(op::QAdd, n::Int = get_order(op))
    n > get_order(op) && return 0
    out = 0
    for (term, coeff) in op.arguments
        out = out + coeff * _term_cumulant(term.ops, n)
    end
    return out
end

function cumulant(avg::SymbolicUtils.BasicSymbolic, args...)
    if SQA.is_average(avg)
        return cumulant(SQA.undo_average(avg), args...)
    end
    throw(ArgumentError("cumulant expects an Average or QField"))
end

function _term_cumulant(ops::Vector, n::Int)
    n > length(ops) && return 0
    acc = 0
    # n-th joint cumulant: sum over partitions of [ops] into 1, 2, ..., n blocks,
    # each block weighted by (m-1)! (-1)^(m-1).
    for k in 1:n
        for p in partitions(ops, k)
            m = length(p)
            coeff = factorial(m - 1) * (-1)^(m - 1)
            prod = 1
            for block in p
                prod = prod * average(_prod_ops(block))
            end
            acc = acc + coeff * prod
        end
    end
    return acc
end

"""
    cumulant_expansion(expr, order; mix_choice=maximum)

Expand averages with order > `order` in terms of lower-order moments using the
joint-cumulant identity (truncating at the supplied order). The output is the
raw expression. Apply `SymbolicUtils.simplify` yourself if you need a canonical
form.
"""
cumulant_expansion(x::Number, order; kw...) = x
cumulant_expansion(x::Number, ::Nothing; kw...) = x
cumulant_expansion(x::SymbolicUtils.BasicSymbolic, ::Nothing; kw...) = x
cumulant_expansion(x::Symbolics.Num, ::Nothing; kw...) = x
cumulant_expansion(eqs::MeanFieldEquations, ::Nothing; kw...) = eqs

# Returns true if the expression tree contains any averaged subexpression.
function _has_average(x::SymbolicUtils.BasicSymbolic)
    SQA.is_average(x) && return true
    SymbolicUtils.iscall(x) || return false
    for a in SymbolicUtils.arguments(x)
        a isa SymbolicUtils.BasicSymbolic && _has_average(a) && return true
    end
    return false
end
_has_average(::Any) = false

function cumulant_expansion(
        x::SymbolicUtils.BasicSymbolic, order::Int;
        mix_choice = maximum
    )
    if _is_leaf_average(x)
        get_order(x) <= order && return x
        return _expand_average(SQA.undo_average(x), order)
    end
    if SymbolicUtils.iscall(x) && _has_average(x)
        op = SymbolicUtils.operation(x)
        args = SymbolicUtils.arguments(x)
        new_args = Any[cumulant_expansion(a, order; mix_choice) for a in args]
        return op(new_args...)
    end
    return x
end

function cumulant_expansion(
        x::SymbolicUtils.BasicSymbolic, order::Vector{Int};
        mix_choice = maximum
    )
    if _is_leaf_average(x)
        ops = SQA.undo_average(x)
        aons = SQA.acts_on(ops)
        ord = isempty(aons) ? maximum(order) : mix_choice(order[k] for k in aons)
        get_order(x) <= ord && return x
        return _expand_average(ops, order; mix_choice)
    end
    if SymbolicUtils.iscall(x) && _has_average(x)
        op = SymbolicUtils.operation(x)
        args = SymbolicUtils.arguments(x)
        new_args = Any[cumulant_expansion(a, order; mix_choice) for a in args]
        return op(new_args...)
    end
    return x
end

cumulant_expansion(x::Symbolics.Num, order; kw...) =
    cumulant_expansion(SymbolicUtils.unwrap(x), order; kw...)

function _expand_average(ops::QField, order::Int)
    out = 0
    if ops isa QAdd
        for (term, coeff) in ops.arguments
            piece = coeff * _expand_product(term.ops, order)
            piece = _stamp_sum_to_first_leaves(piece, ops.indices, term.ne)
            out = out + piece
        end
    else
        out = average(ops)
    end
    return out
end

# Per-Hilbert-space variant: each emitted sub-block picks its own per-subspace
# `ord` from `acts_on(block)`, not the outer aggregate. Required for mixed
# `order = [m, n]` truncation where an atom-only length-2 block emitted from a
# mixed expansion still needs to expand under the atom cap.
function _expand_average(ops::QField, order::Vector{Int}; mix_choice = maximum)
    out = 0
    if ops isa QAdd
        for (term, coeff) in ops.arguments
            piece = coeff * _expand_product(term.ops, order; mix_choice)
            piece = _stamp_sum_to_first_leaves(piece, ops.indices, term.ne)
            out = out + piece
        end
    else
        out = average(ops)
    end
    return out
end

# After cumulant factorization, `Σ_{i_1,...,i_k} ⟨A_1⟩⟨A_2⟩⋯` would lose the
# outer sum scope (Σ is encoded as `SumIndices`/`SumNonEqual` metadata on the
# average sym, and the factored-out product is no longer a single average).
# Re-attach the sum scope by stamping each bound index onto the first averaged
# leaf in the product that references it. Bound indices that do not appear in
# any leaf are treated as "spurious" (factor 1 in `scale`'s prefactor logic),
# matching SQA's pre-existing convention for `_accumulate_with_diag!`-produced
# leftover scope.
function _stamp_sum_to_first_leaves(piece, indices::Vector{SQA.Index}, ne)
    isempty(indices) && return piece
    piece isa SymbolicUtils.BasicSymbolic || return piece
    # The factorised `piece` is a SUM of Wick products. The sum scope distributes
    # over the additive terms (`Σ_i (A + B) = Σ_i A + Σ_i B`), and within each
    # product `Σ_i` binds that product's own indexed leaf. Recurse per additive
    # term so EVERY Wick term referencing a bound index carries the scope on its
    # leaf, not just the globally-first one (which left e.g. `⟨a'σ_i⟩⟨a⟩` without
    # its `Σ_i` while `⟨σ_i⟩⟨a'a⟩` kept it).
    if SymbolicUtils.iscall(piece) && SymbolicUtils.operation(piece) === (+)
        terms = SymbolicUtils.arguments(piece)
        stamped = Any[_stamp_sum_to_first_leaves(t, indices, ne) for t in terms]
        return reduce(+, stamped)
    end
    leaves = SymbolicUtils.BasicSymbolic[]
    _collect_avg_leaves!(leaves, piece)
    isempty(leaves) && return piece
    leaf_idx_assign = [SQA.Index[] for _ in leaves]
    leaf_ne_assign = [Tuple{SQA.Index, SQA.Index}[] for _ in leaves]
    for idx in indices
        slot = _first_leaf_using(leaves, idx)
        slot === nothing && continue
        push!(leaf_idx_assign[slot], idx)
        # Carry along NE pairs that involve only indices already assigned to
        # this leaf (so the metadata stays self-consistent per leaf).
        for pair in ne
            (pair[1] == idx || pair[2] == idx) || continue
            (pair[1] in leaf_idx_assign[slot] && pair[2] in leaf_idx_assign[slot]) ||
                continue
            push!(leaf_ne_assign[slot], pair)
        end
    end
    sub = IdDict{Any, Any}()
    for (i, leaf) in enumerate(leaves)
        isempty(leaf_idx_assign[i]) && continue
        # Reconstruct the summed average through `average(SQA.Σ(op, idxs...))`, the
        # SAME canonical form a freshly-derived sum carries. A bare
        # `setmetadata(average(op), SumIndices, …)` is display- and isequal-equal to
        # the un-summed leaf but is NOT isequal to the canonical summed average
        # (their inner structure differs), so the two never cancel in arithmetic
        # (`Σ_i⟨σ_i⟩⟨a†a⟩ - ⟨a†a⟩Σ_i⟨σ_i⟩ ≠ 0`). Building it canonically makes the
        # cumulant-factored sum identical to the natural one.
        op = undo_average(leaf)
        idxs = leaf_idx_assign[i]
        summed = SQA.Σ(op, idxs[1], idxs[2:end]...)
        new_leaf = SymbolicUtils.unwrap(average(summed))
        # `SQA.Σ` over a plain operator does not carry an NE constraint; re-attach
        # the per-leaf NE pairs as metadata when the factorisation produced any.
        isempty(leaf_ne_assign[i]) ||
            (new_leaf = SymbolicUtils.setmetadata(new_leaf, SQA.SumNonEqual, copy(leaf_ne_assign[i])))
        sub[leaf] = new_leaf
    end
    isempty(sub) && return piece
    return _substitute_by_identity(piece, sub)
end

function _collect_avg_leaves!(out, x)
    x isa SymbolicUtils.BasicSymbolic || return out
    if _is_leaf_average(x)
        push!(out, x)
        return out
    end
    SymbolicUtils.iscall(x) || return out
    for a in SymbolicUtils.arguments(x)
        _collect_avg_leaves!(out, a)
    end
    return out
end

function _first_leaf_using(leaves, idx::SQA.Index)
    for (i, leaf) in enumerate(leaves)
        op = SQA.undo_average(leaf)
        if op isa QAdd
            for (term, _) in op.arguments, o in term.ops
                if SQA.has_index(o.index) && o.index == idx
                    return i
                end
            end
        end
    end
    return nothing
end

# Substitute by object identity (the leaves were obtained from `piece` itself).
function _substitute_by_identity(x, sub::IdDict)
    haskey(sub, x) && return sub[x]
    x isa SymbolicUtils.BasicSymbolic || return x
    SymbolicUtils.iscall(x) || return x
    op = SymbolicUtils.operation(x)
    args = SymbolicUtils.arguments(x)
    new_args = Any[_substitute_by_identity(a, sub) for a in args]
    all(((a, b),) -> a === b, zip(args, new_args)) && return x
    return op(new_args...)
end

function _expand_product(args::Vector, order::Int)
    n = length(args)
    n <= order && return average(_prod_ops(args))
    acc = 0
    for k in 2:n
        for p in partitions(args, k)
            m = length(p)
            coeff = -factorial(m - 1) * (-1)^(m - 1)
            prod = 1
            for block in p
                if length(block) > order
                    prod = prod * _expand_product(block, order)
                else
                    prod = prod * average(_prod_ops(block))
                end
            end
            acc = acc + coeff * prod
        end
    end
    return acc
end

function _expand_product(args::Vector, order::Vector{Int}; mix_choice = maximum)
    p = _prod_ops(args)
    aons = SQA.acts_on(p)
    ord = isempty(aons) ? maximum(order) : mix_choice(order[k] for k in aons)
    n = length(args)
    n <= ord && return average(p)
    acc = 0
    for k in 2:n
        for part in partitions(args, k)
            m = length(part)
            coeff = -factorial(m - 1) * (-1)^(m - 1)
            prod = 1
            for block in part
                prod = if length(block) == 1
                    prod * average(block[1])
                else
                    prod * _expand_product(block, order; mix_choice)
                end
            end
            acc = acc + coeff * prod
        end
    end
    return acc
end

"""
    cumulant_expansion(eqs::AbstractMeanFieldEquations, order; ...)

Apply `cumulant_expansion` to every RHS in `eqs.equations`, returning a new
struct with the same shape and the supplied order stored in `eqs.order`.
"""
function cumulant_expansion(
        eqs::MeanFieldEquations, order;
        mix_choice = maximum
    )
    order_vec = _normalize_order(order, eqs)
    new_eqs = [
        eq.lhs ~ cumulant_expansion(eq.rhs, order_vec; mix_choice)
            for eq in eqs.equations
    ]
    return MeanFieldEquations(
        new_eqs, eqs.operator_equations, eqs.states,
        eqs.operators, eqs.hamiltonian, eqs.jumps,
        eqs.jumps_dagger, eqs.rates, eqs.iv, order_vec;
        initial_operators = copy(eqs.initial_operators),
    )
end

function _normalize_order(order::Int, eqs)
    nspaces = _nspaces(eqs.hamiltonian)
    return fill(order, nspaces)
end
_normalize_order(order::Vector{Int}, _) = order
_normalize_order(::Nothing, _) = nothing

function _nspaces(op::QField)
    aons = SQA.acts_on(op)
    return isempty(aons) ? 1 : maximum(aons)
end
