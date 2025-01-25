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

    args = SymbolicUtils.sorted_arguments(t)
    f = SymbolicUtils.operation(t)
    h = index_invariant_hash(f, h0)
    for arg in args
        h = index_invariant_hash(arg, h)
    end
   
    return h
end

function _index_invariant_hash(q::QMul, h0::UInt)
    inds = get_indices(q)
    if isempty(inds)
        return hash(t, h0)
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