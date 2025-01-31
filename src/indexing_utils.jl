# TODO: move this function
function SymbolicUtils.sorted_arguments(t::SymbolicUtils.BasicSymbolic{<:CNumber})
    f = SymbolicUtils.operation(t)
    args = SymbolicUtils.arguments(t)

    if f in (*, +)
        # commutative for numbers, we may sort
        sort!(args, by=x -> nameof(SymbolicUtils.symtype(x)))
    end

    return args
end

index_invariant_hash(t::SymbolicUtils.Symbolic, h0::UInt) = _index_invariant_hash(t, h0)
index_invariant_hash(t::QNumber, h0::UInt) = _index_invariant_hash(t, h0)

function _index_invariant_hash(t, h0::UInt)
    if !TermInterface.iscall(t)
        return hash(t, h0)
    end

    inds = get_indices(t)
    if isempty(inds)
        return hash(t, h0)
    end

    f = SymbolicUtils.operation(t)
    args = SymbolicUtils.sorted_arguments(t)

    h = index_invariant_hash(f, h0)
    for arg in args
        h = index_invariant_hash(arg, h)
    end
   
    return h
end

function _index_invariant_hash(q::QMul, h0::UInt)
    inds = get_indices(q)
    if isempty(inds)
        return hash(q, h0)
    end

    h = hash((*), h0)
    h = index_invariant_hash(q.arg_c, h)
    args_nc = copy(q.args_nc)
    if was_all_merged(args_nc)
        # TODO: may have to change the checks here, e.g. for a'*s_i*s_j or so
        # could probably be done by splitting args by acts_on and sorting only each subset
        # order doesn't matter; to get consistent one sort by transition levels
        sort!(args_nc, by=_transition_levels)
    end

    for arg in args_nc
        h = index_invariant_hash(arg, h)
    end

    return h
end

function was_all_merged(args)
    for arg1 in args
        for arg2 in args
            _was_merged(arg1, arg2) || return false
        end
    end
    return true
end
_was_merged(a, b) = false
_was_merged(a::IndexedOperator{<:Transition}, b::IndexedOperator{<:Transition}) = was_merged(a, b)

const INDEX_HASH = 0x7ccb9f972ab655dc  # random UInt, but always the same
index_invariant_hash(::Index, h::UInt) = hash(INDEX_HASH, h)

function index_invariant_hash(op::IndexedOperator, h::UInt)
    op_hash = hash(op.op, h)
    return index_invariant_hash(op.ind, op_hash)
end

function _transition_levels(op::IndexedOperator{<:Transition})
    i = op.op.i
    j = op.op.j
    return i, j
end


function switch_to_extra_indices!(vs, indices_in_use, extra_indices)
    # change all indices that are already in use in e.g. the Hamiltonian
    # in an expression so that we can derive equations for it
    # the conditions for a valid replacement are:
    # * the index is not in use in the Hamiltonian (indices_in_use)
    # * the index is not in use in the current expression
    # * the index is not used as a replacement for another index
    for i=1:length(vs)
        inds = get_indices(vs[i])
        inds_to_switch = filter(in(indices_in_use), inds)
        isempty(inds_to_switch) && continue

        # build mapping for all indices
        # need to make sure indices are not re-used
        # the conditions for an index 
        to_index = Index[]
        for from in inds_to_switch
            to_position = findfirst(k -> !(k in indices_in_use) && !(k in to_index) && !(k in inds), extra_indices)
            
            (to_position === nothing) && throw("Not enough extra indices provided!")

            to = extra_indices[to_position]
            push!(to_index, to)
        end

        for (from, to) in zip(inds_to_switch, to_index)
            vs[i] = change_index(vs[i], from, to)
        end
    end

    return vs
end



function evaluate(de::MeanfieldEquations; limits=Dict{Any,Int}())
    indices = get_indices(de)
    isempty(indices) && return de

    # First, we generate indices that have actual integer upper bounds
    new_indices = _substitute_limits(indices, limits)
    index_map = Dict(indices .=> new_indices)

    # Next, we substitute the original indices in the equations with the ones that now have integer limits
    equations_with_actual_limits = MTK.Equation[]
    for eq in de.equations
        lhs = deepcopy(eq.lhs)
        rhs = deepcopy(eq.rhs)
        for (from, to) in index_map
            lhs = change_index(lhs, from, to)
            rhs = change_index(rhs, from, to)
        end

        push!(equations_with_actual_limits, MTK.Equation(lhs, rhs))
    end

    # At this point, all indices that occur have an actual upper bound that is an integer
    # So, now we only need to replace each index by 1:N in all expressions and expand the sums as well
    expanded_equations = MTK.Equation[]
    expanded_lhs = []
    for eq in equations_with_actual_limits
        lhs_original_indices = get_indices(eq.lhs)
        lhs = _expand_lhs_to_index_limits(eq.lhs, lhs_original_indices)
        for l in lhs
            if any(isequal(l), expanded_lhs)
                # dirty if here -- skip duplicates
                # TODO: check for adjoints as well
                continue
            end

            lhs_integer_indices = get_indices(l)
            lhs_index_map = Dict(lhs_original_indices .=> lhs_integer_indices)
            r = deepcopy(eq.rhs)
            for (from, to) in lhs_index_map
                r = change_index(r, from, to)
            end
            r = _expand_sums(r)  # TODO: this is stupid and slow; should just generate actual sum(term for i=1:N) expressions later on
            sort_by_integer_indices!(r)
            push!(expanded_equations, MTK.Equation(l, r))
        end
        append!(expanded_lhs, lhs)
    end

    # TODO: this can probably be skipped -- also, it's probably wrong since we already did the cumulant expansion
    operators = map(undo_average, [eq.lhs for eq in expanded_equations])
    op_rhs = map(undo_average, [eq.rhs for eq in expanded_equations])
    operator_equations = [MTK.Equation(l, r) for (l, r) in zip(operators, op_rhs)]

    vs = [eq.lhs for eq in expanded_equations]
    varmap = make_varmap(vs, de.iv)

    return MeanfieldEquations(
        expanded_equations,
        operator_equations,
        vs,
        operators,
        de.hamiltonian,
        de.jumps,
        de.jumps_dagger,
        de.rates,
        de.iv,
        varmap,
        de.order
    )
end

_substitute_limits(indices::Set, limits) = [_substitute_limits(index, limits) for index in indices]
function _substitute_limits(index::Index, limits)
    N = limits[index.range]
    return Index(index.hilb, index.name, N, index.aon)
end

function _expand_lhs_to_index_limits(lhs, indices)
    isempty(indices) && return [lhs]

    new_lhs = []
    iterators = [1:(index.range) for index in indices]

    for to_indices in Iterators.product(iterators...)
        new_lhs_ = deepcopy(lhs)
        for (from, to) in zip(indices, to_indices)
            new_lhs_ = change_index(new_lhs_, from, to)
        end
        SymbolicUtils._iszero(new_lhs_) || push!(new_lhs, new_lhs_)
    end

    for l in new_lhs
        sort_by_integer_indices!(l)
    end

    # TODO: this is only needed because we loop over the full range for all indices, so end up generating
    # e.g. s_1*s_2 and s_2*s_1 from products such as s_i*s_j; need to be smarter here
    unique_ops!(new_lhs)

    return new_lhs
end

function sort_by_integer_indices!(t)
    if !TermInterface.iscall(t)
        return nothing
    end

    for arg in SymbolicUtils.arguments(t)
        sort_by_integer_indices!(arg)
    end

    return nothing
end

function sort_by_integer_indices!(t::QMul)
    args_nc = split_by_aon(t.args_nc)

    for arg_vec in args_nc
        sort_by_integer_indices!(arg_vec)
    end

    new_args_nc = vcat(args_nc...)
    for i=1:length(t.args_nc)
        t.args_nc[i] = new_args_nc[i]
    end

    return nothing
end

function split_by_aon(args::Vector)
    all_args = []
    aon = acts_on(args[1])

    chunk_index = findfirst(x -> acts_on(x) > aon, args)
    previous_chunk_index = 1
    while chunk_index !== nothing
        chunk = args[previous_chunk_index:chunk_index-1]
        push!(all_args, chunk)

        previous_chunk_index = copy(chunk_index)
        aon = acts_on(args[chunk_index])
        chunk_index = findfirst(x -> acts_on(x) > aon, args)
    end

    # add the last chunk
    push!(all_args, args[previous_chunk_index:end])

    return all_args
end
    
function sort_by_integer_indices!(args::Vector)
    _get_index(arg) = 0
    _get_index(arg::IndexedOperator{T,I}) where {T, I<:Integer} = arg.ind
    sort!(args, by=_get_index)
    return nothing
end



function _expand_sums(t)
    if !TermInterface.iscall(t)
        return t
    end

    f = SymbolicUtils.operation(t)
    args = [_expand_sums(arg) for arg in SymbolicUtils.arguments(t)]
    return f(args...)
end

function _expand_sums(s::SymbolicUtils.BasicSymbolic{<:CSumSym})
    term, index = SymbolicUtils.arguments(s)
    args = [change_index(term, index, i) for i=1:index.range]
    return +(args...)
end