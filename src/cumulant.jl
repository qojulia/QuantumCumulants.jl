"""
    get_order(expr)

Highest moment order appearing in `expr`, the order used to decide whether a term is
expanded by [`cumulant_expansion`](@ref). Numbers have order 0, a bare operator order 1,
and a product the number of its operator factors.

# Examples
```jldoctest
julia> h = FockSpace(:a) ⊗ FockSpace(:b);

julia> @qnumbers a::Destroy(h, 1) b::Destroy(h, 2);

julia> get_order(a)
1

julia> get_order(a * b)
2

julia> get_order(1)
0
```
"""
get_order(::Number) = 0
get_order(::SQA.QSym) = 1
function get_order(q::QAdd)
    isempty(q.arguments) && return 0
    return maximum(length(t.ops) for (t, _) in q.arguments)
end
function get_order(x::SymbolicUtils.BasicSymbolic)
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

_prod_ops(block::AbstractVector) = isempty(block) ? 1 : reduce(*, block)

"""
    cumulant(op, n=get_order(op))

The `n`-th joint cumulant of `op` (either an operator or an average), the signed sum
over partitions of its operator factors into `n` blocks. The result is raw and
unsimplified; apply `SymbolicUtils.simplify` for a canonical form.

# Examples
```jldoctest
julia> h = FockSpace(:a) ⊗ FockSpace(:b);

julia> @qnumbers a::Destroy(h, 1) b::Destroy(h, 2);

julia> cumulant(a * b)
⟨a * b⟩ - ⟨b⟩*⟨a⟩

julia> cumulant(a * b, 1)
⟨a * b⟩

julia> cumulant(a * b, 3)
0
```
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

Expand every average in `expr` whose order exceeds `order` into lower-order moments via
the joint-cumulant identity, neglecting the joint cumulant above `order`. `order` is an
`Int` (one cap for all subspaces) or a `Vector{Int}` (a cap per Hilbert subspace,
combined by `mix_choice` on mixed terms). The output is raw; apply `SymbolicUtils.simplify`
for a canonical form.

See also: https://en.wikipedia.org/wiki/Cumulant#Joint_cumulants

# Examples
```jldoctest
julia> h = FockSpace(:a) ⊗ FockSpace(:b);

julia> @qnumbers a::Destroy(h, 1) b::Destroy(h, 2);

julia> cumulant_expansion(average(a * b), 1)
⟨b⟩*⟨a⟩
```
"""
cumulant_expansion(x::Number, order; kw...) = x
cumulant_expansion(x::SymbolicUtils.BasicSymbolic, ::Nothing; kw...) = x
cumulant_expansion(x::Symbolics.Num, ::Nothing; kw...) = x
cumulant_expansion(eqs::MeanFieldEquations, ::Nothing; kw...) = eqs

"""
True when the expression tree contains any averaged subexpression.
"""
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
    if _is_avg_leaf(x)
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
    if _is_avg_leaf(x)
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

"""
Per-Hilbert-space variant of `_expand_average`: each emitted sub-block picks its
own per-subspace order from `acts_on(block)` rather than the outer aggregate.
Required for mixed `order = [m, n]` truncation, where an atom-only length-2 block
emitted from a mixed expansion still needs to expand under the atom cap.
"""
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

"""
Re-attach the outer sum scope lost during cumulant factorisation. After
`Σ_{i_1,…,i_k} ⟨A_1⟩⟨A_2⟩⋯` factorises, the `Σ` (encoded as
`SumIndices`/`SumNonEqual` metadata on a single average sym) is dropped because
the factored-out product is no longer one average. Stamp each bound index onto
the first averaged leaf in the product that references it. Bound indices that
appear in no leaf are treated as spurious (factor 1 in `scale`'s prefactor
logic), matching SQA's convention for `_accumulate_with_diag!`-produced leftover
scope.
"""
function _stamp_sum_to_first_leaves(piece, indices::Vector{SQA.Index}, non_equal)
    isempty(indices) && return piece
    piece isa SymbolicUtils.BasicSymbolic || return piece
    # Recurse per additive term: the sum scope distributes over `+`, so every Wick term
    # that references a bound index must carry it on its own leaf, not just the first
    # term's (else e.g. `⟨a'σ_i⟩⟨a⟩` loses its `Σ_i` while `⟨σ_i⟩⟨a'a⟩` keeps it).
    if SymbolicUtils.iscall(piece) && SymbolicUtils.operation(piece) === (+)
        terms = SymbolicUtils.arguments(piece)
        stamped = Any[_stamp_sum_to_first_leaves(t, indices, non_equal) for t in terms]
        return reduce(+, stamped)
    end
    leaves = SymbolicUtils.BasicSymbolic[]
    _collect_avg_leaves!(leaves, piece)
    isempty(leaves) && return piece
    leaf_idx_assign = [SQA.Index[] for _ in leaves]
    leaf_non_equal_assign = [Tuple{SQA.Index, SQA.Index}[] for _ in leaves]
    for idx in indices
        slot = _first_leaf_using(leaves, idx)
        slot === nothing && continue
        push!(leaf_idx_assign[slot], idx)
    end
    # Keep a non-equal pair if both indices land on this leaf, or one does and its
    # partner is external (a `Σ_{j≠ext}` constraint); a pair split across two bound
    # leaves is dropped (the bound-vs-bound truncation case).
    for (slot, bidxs) in enumerate(leaf_idx_assign)
        isempty(bidxs) && continue
        for pair in non_equal
            a_in = pair[1] in bidxs
            b_in = pair[2] in bidxs
            keep = (a_in && b_in) ||
                (a_in && !(pair[2] in indices)) ||
                (b_in && !(pair[1] in indices))
            keep && push!(leaf_non_equal_assign[slot], pair)
        end
    end
    sub = IdDict{Any, Any}()
    for (i, leaf) in enumerate(leaves)
        isempty(leaf_idx_assign[i]) && continue
        # Rebuild via the canonical `average(SQA.Σ(op, idxs...))`: a plain
        # `setmetadata(average(op), SumIndices, …)` is isequal to the un-summed leaf, not
        # to the canonical summed form, so the two never cancel
        # (`Σ_i⟨σ_i⟩⟨a†a⟩ - ⟨a†a⟩Σ_i⟨σ_i⟩ ≠ 0`).
        op = undo_average(leaf)
        idxs = leaf_idx_assign[i]
        summed = SQA.Σ(op, idxs[1], idxs[2:end]...)
        new_leaf = SymbolicUtils.unwrap(average(summed))
        # `SQA.Σ` over a plain operator does not carry a non-equal index constraint; re-attach
        # the per-leaf non-equal index pairs as metadata when the factorisation produced any.
        isempty(leaf_non_equal_assign[i]) ||
            (new_leaf = SymbolicUtils.setmetadata(new_leaf, SQA.SumNonEqual, copy(leaf_non_equal_assign[i])))
        sub[leaf] = new_leaf
    end
    isempty(sub) && return piece
    return _substitute_by_identity(piece, sub)
end

function _collect_avg_leaves!(out, x)
    x isa SymbolicUtils.BasicSymbolic || return out
    if _is_avg_leaf(x)
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

"""
The cumulant-truncation kernel: rewrite ⟨a product of `args`⟩ that exceeds `order` as
the signed sum over partitions into 2..n blocks, dropping the joint cumulant above
`order`. Each partition of `m` blocks carries the moment-cumulant inversion weight
`-(m-1)!·(-1)^(m-1)`; a block still longer than `order` is expanded recursively, a block
within the cap is averaged directly. A product already within `order` is left as one
average.
"""
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

"""
Per-Hilbert-space variant of `_expand_product`: the truncation cap for each block is its
own `mix_choice(order[acts_on])` rather than a single global `order`, so a mixed product
is capped by the subspaces it actually touches.
"""
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
    cumulant_expansion(eqs::MeanFieldEquations, order; mix_choice=maximum)

Cumulant-expand every RHS of `eqs` to `order`, returning a new `MeanFieldEquations` of
the same shape with that order stored. Noise systems take their order through
[`meanfield`](@ref), not here.
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
        eqs.jumps_dagger, eqs.rates, eqs.iv, order_vec, eqs.direction;
        treatments = eqs.treatments,
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
