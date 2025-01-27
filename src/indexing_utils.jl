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