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
    _is_avg_leaf(x) && return get_order(SQA.undo_average(x))
    SymbolicUtils.iscall(x) && return maximum(get_order, SymbolicUtils.arguments(x); init = 0)
    return 0
end
get_order(x::Symbolics.Num) = get_order(SymbolicUtils.unwrap(x))

_prod_ops(block::AbstractVector) = isempty(block) ? 1 : reduce(*, block)

# Bulk-build a sum/product from a Vector of terms in one SymbolicUtils pass. We call the
# add/mul workers on the Vector directly rather than splatting `+(terms...)`: splatting
# materialises an `NTuple{N}` whose type depends on `N`, so a fresh expansion size forces
# SymbolicUtils to recompile the worker (catastrophic at large `N`); a Vector keeps one
# compiled method for every size. The workers accept "an indexable list" by their own API.
function _vartype(xs)
    for x in xs
        x isa SymbolicUtils.BasicSymbolic && return SymbolicUtils.vartype(typeof(x))
    end
    return SymbolicUtils.SymReal
end
_bulk_add(terms::Vector) =
    isempty(terms) ? 0 :
    length(terms) == 1 ? terms[1] : SymbolicUtils.add_worker(_vartype(terms), terms)
_bulk_mul(factors::Vector) = SymbolicUtils.mul_worker(_vartype(factors), factors)

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
    terms = Any[]
    for (term, coeff) in op.arguments
        push!(terms, _im_form(_coeff_num(coeff)) * _term_cumulant(term.ops, n))
    end
    return _bulk_add(terms)
end

function cumulant(avg::SymbolicUtils.BasicSymbolic, args...)
    if SQA.is_average(avg)
        return cumulant(SQA.undo_average(avg), args...)
    end
    throw(ArgumentError("cumulant expects an Average or QField"))
end

function _term_cumulant(ops::Vector, n::Int)
    n > length(ops) && return 0
    terms = Any[]
    leaves = Dict{Any, Any}()   # the same block recurs across partitions; build its moment once
    # n-th joint cumulant: sum over partitions of [ops] into 1, 2, ..., n blocks,
    # each block weighted by (m-1)! (-1)^(m-1).
    for k in 1:n
        for p in partitions(ops, k)
            m = length(p)
            coeff = factorial(m - 1) * (-1)^(m - 1)
            factors = Any[coeff]
            for block in p
                push!(factors, get!(() -> average(_prod_ops(block)), leaves, Tuple(block)))
            end
            push!(terms, _bulk_mul(factors))
        end
    end
    return _bulk_add(terms)
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
cumulant_expansion(eqs::MeanfieldEquations, ::Nothing; kw...) = eqs

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
    # Whole-expression early-out: if no moment exceeds the cap, nothing expands.
    get_order(x) <= order && return x
    if _is_moment_unit(x)
        return _expand_average(SQA.undo_average(x), order)
    end
    if SymbolicUtils.iscall(x) && _has_average(x)
        op = SymbolicUtils.operation(x)
        args = SymbolicUtils.arguments(x)
        new_args = Any[cumulant_expansion(a, order; mix_choice) for a in args]
        all(i -> new_args[i] === args[i], eachindex(args)) && return x
        return op(new_args...)
    end
    return x
end

function cumulant_expansion(
        x::SymbolicUtils.BasicSymbolic, order::Vector{Int};
        mix_choice = maximum
    )
    # Conservative whole-expression early-out: total order <= the smallest per-subspace
    # cap guarantees no leaf can exceed its own cap.
    get_order(x) <= minimum(order) && return x
    if _is_moment_unit(x)
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
        all(i -> new_args[i] === args[i], eachindex(args)) && return x
        return op(new_args...)
    end
    return x
end

cumulant_expansion(x::Symbolics.Num, order; kw...) =
    cumulant_expansion(SymbolicUtils.unwrap(x), order; kw...)

"""
Expand each Wick term of `ops` to `order` and re-attach its sum scope. For a vector
`order = [m, n]`, each emitted sub-block picks its own per-subspace cap from
`acts_on(block)` rather than the outer aggregate, so an atom-only length-2 block emitted
from a mixed expansion still expands under the atom cap. `mix_choice` is ignored for a
scalar `order`.
"""
function _expand_average(ops::QField, order; mix_choice = maximum)
    ops isa QAdd || return average(ops)
    terms = Any[]
    for (term, coeff) in ops.arguments
        piece = _im_form(_coeff_num(coeff)) * _expand_product(term.ops, order; mix_choice)
        piece = _stamp_sum_to_first_leaves(piece, ops.indices, term.ne)
        push!(terms, piece)
    end
    return _bulk_add(terms)
end

"""
Re-attach the outer sum scope lost during cumulant factorisation. After
`Σ_{i_1,…,i_k} ⟨A_1⟩⟨A_2⟩⋯` factorises, the indexed sum is dropped because the
factored-out product is no longer one average. Stamp each bound index onto the first
averaged leaf in the product that references it, rebuilding that leaf as a canonical
`average(Σ(op, idx…))` indexed-sum node. Bound indices that appear in no leaf are treated
as spurious (factor 1 in `scale`'s prefactor logic), matching SQA's convention for
leftover scope.
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
        return _bulk_add(stamped)
    end
    leaves = eachleaf(piece)
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
        # Rebuild as the canonical `average(SQA.Σ(op, idx…))` indexed-sum node so it shares
        # structure with sums emitted elsewhere and cancels with them
        # (`Σ_i⟨σ_i⟩⟨a†a⟩ - ⟨a†a⟩Σ_i⟨σ_i⟩ = 0`).
        op = undo_average(leaf)
        idxs = leaf_idx_assign[i]
        summed = SQA.Σ(op, idxs[1], idxs[2:end]...)
        # `SQA.Σ` over a plain operator carries no non-equal index constraint; merge the
        # per-leaf pairs into the summed term so `average` keeps them in the sum's scope.
        isempty(leaf_non_equal_assign[i]) ||
            (summed = _carry_non_equal(summed, leaf_non_equal_assign[i], idxs))
        new_leaf = SymbolicUtils.unwrap(average(summed))
        sub[leaf] = new_leaf
    end
    isempty(sub) && return piece
    return _substitute_by_identity(piece, sub)
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
_substitute_by_identity(x, sub::IdDict) = _subtree_substitute(x, sub)

"""
The cumulant-truncation kernel: rewrite ⟨a product of `args`⟩ that exceeds `order` as
the signed sum over partitions into 2..n blocks, dropping the joint cumulant above
`order`. Each partition of `m` blocks carries the moment-cumulant inversion weight
`-(m-1)!·(-1)^(m-1)`; a block still longer than `order` is expanded recursively, a block
within the cap is averaged directly. A product already within `order` is left as one
average. `mix_choice` is accepted for signature parity with the vector method and ignored:
a global cap needs no per-subspace choice.
"""
_expand_product(args::Vector, order::Int; mix_choice = maximum) =
    _expand_product!(args, order, Dict{Any, Any}())

# Memoised worker. `_expand_product` of an ordered block at a fixed `order` is a pure
# function of the two, but the same sub-block recurs across thousands of partitions: a
# length-7 order-2 product reaches only 99 distinct blocks across 7022 recursive calls and
# 28 distinct moment leaves across ~98k constructions. Caching by block content collapses
# the redundant SymbolicUtils construction, the dominant cost (~17x on length-7 order-2).
function _expand_product!(args::Vector, order::Int, memo::Dict)
    key = Tuple(args)
    cached = get(memo, key, nothing)
    cached === nothing || return cached
    n = length(args)
    if n <= order
        res = average(_prod_ops(args))
    else
        terms = Any[]
        for k in 2:n
            for p in partitions(args, k)
                m = length(p)
                coeff = -factorial(m - 1) * (-1)^(m - 1)
                factors = Any[coeff]
                for block in p
                    push!(factors, _expand_product!(block, order, memo))
                end
                push!(terms, _bulk_mul(factors))   # one Mul per partition, built once
            end
        end
        res = _bulk_add(terms)   # one Add over all partitions
    end
    memo[key] = res
    return res
end

"""
Per-Hilbert-space variant of `_expand_product`: the truncation cap for each block is its
own `mix_choice(order[acts_on])` rather than a single global `order`, so a mixed product
is capped by the subspaces it actually touches.
"""
_expand_product(args::Vector, order::Vector{Int}; mix_choice = maximum) =
    _expand_product!(args, order, mix_choice, Dict{Any, Any}())

function _expand_product!(args::Vector, order::Vector{Int}, mix_choice, memo::Dict)
    key = Tuple(args)
    cached = get(memo, key, nothing)
    cached === nothing || return cached
    p = _prod_ops(args)
    aons = SQA.acts_on(p)
    ord = isempty(aons) ? maximum(order) : mix_choice(order[k] for k in aons)
    n = length(args)
    if n <= ord
        res = average(p)
    else
        terms = Any[]
        for k in 2:n
            for part in partitions(args, k)
                m = length(part)
                coeff = -factorial(m - 1) * (-1)^(m - 1)
                factors = Any[coeff]
                for block in part
                    push!(factors, _expand_product!(block, order, mix_choice, memo))
                end
                push!(terms, _bulk_mul(factors))
            end
        end
        res = _bulk_add(terms)
    end
    memo[key] = res
    return res
end

"""
    cumulant_expansion(eqs::MeanfieldEquations, order)

Cumulant-expand every RHS of `eqs` to `order`, returning a new `MeanfieldEquations` of
the same shape with that order stored. The system's `mix_choice` (held in `eqs.graph.sys`)
drives the mixed-order truncation. Noise systems take their order through
[`meanfield`](@ref), not here.
"""
function cumulant_expansion(eqs::MeanfieldEquations, order)
    order_vec = _normalize_order(order, eqs)
    sys = eqs.graph.sys
    # Idempotency: equations already expanded to exactly this order need no rebuild.
    sys.order == order_vec && return eqs
    g = map_drifts(
        eqs.graph, (_, d) -> cumulant_expansion(d, order_vec; mix_choice = sys.mix_choice);
        noise = false,
    )
    new_sys = SystemSpec(
        sys.hamiltonian, sys.jumps, sys.jumps_dagger, sys.rates, sys.efficiencies,
        sys.iv, order_vec, sys.mix_choice, sys.direction,
    )
    return MeanfieldEquations(MomentGraph(g.nodes, new_sys, g.ctx, g.treatments))
end

_normalize_order(order::Int, eqs) = fill(order, _nspaces(eqs))
_normalize_order(order::Vector{Int}, _) = order
_normalize_order(::Nothing, _) = nothing

function _nspaces(op::QField)
    aons = SQA.acts_on(op)
    return isempty(aons) ? 1 : maximum(aons)
end
function _nspaces(ops::AbstractVector)
    isempty(ops) && return 1
    return maximum(_nspaces, ops; init = 1)
end
# `spec` carries the full system inputs (observables, Hamiltonian, jumps and
# adjoint jumps); count subspaces across all of them so a scalar `order` sizes
# the order vector for subspaces no operator in `H` alone touches. Works for both
# the `meanfield` input NamedTuple and a stored `MeanfieldEquations`, which share
# these field names.
function _nspaces(spec)
    return max(
        _nspaces(spec.operators),
        _nspaces(spec.hamiltonian),
        _nspaces(spec.jumps),
        _nspaces(spec.jumps_dagger),
    )
end
