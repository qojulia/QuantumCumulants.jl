function build_duplicates(de::DifferentialEquation{<:AbstractOperator,<:AbstractOperator})
    idx_ops = find_ncopy(de.lhs,hilbert(de.lhs[1]))
    isempty(idx_ops) && return de
    lhs, rhs = build_duplicates(de.lhs[idx_ops], de.rhs[idx_ops])
    not_idx_ops = filter(x->!(x in idx_ops), 1:length(idx_ops))
    lhs_out = [de.lhs[not_idx_ops]; lhs]
    rhs_out = [de.rhs[not_idx_ops]; rhs]
    return DifferentialEquation(lhs_out,rhs_out)
end

# Find all indexes of operators with a Hilbert space storing multiple copies
function find_ncopy(lhs, hilbert)
    idx = Int[]
    for i=1:length(lhs)
        aon = acts_on(lhs[i])
        if length(aon) < length(hilbert,aon)
            push!(idx, i)
        end
    end
    return idx
end

function build_duplicates(lhs,rhs)
    lhs_ = copy(lhs)
    rhs_ = copy(rhs)
    for i=1:length(lhs)
        idx = get_index(lhs[i])
        aon = acts_on(lhs[i])
        ns = get_sizes(hilbert(lhs[i]), aon)
        combs = Iterators.product([1:n for n in ns]...)
        for js in combs
            for jj=1:length(js)
                l = permute_idx(lhs[i], js[jj], ns[jj], aon[jj])
                (l in lhs_) && continue
                r = permute_idx(rhs[i], js[jj], ns[jj], aon[jj])
                push!(lhs_, l)
                push!(rhs_, r)
            end
        end
    end
    return lhs_, rhs_
end

get_sizes(hilbert::HilbertSpace, aon) = [length(h)]
get_sizes(hilbert::ProductSpace, aon::Int) = [length(hilbert.spaces[aon])]
get_sizes(hilbert::ProductSpace, aon::Vector) = length.(hilbert.spaces[aon])


function permute_idx(op::T,j,n) where T<:BasicOperator
    idx = get_index(op)
    k = idx + j
    while k > n
        k -= n
    end
    fields = [getfield(op, name) for name in fieldnames(T)]
    fields[end] = k
    return T(fields...)
end
function permute_idx(op::OperatorTerm,j,n,aon)
    args = Any[]
    for arg in op.arguments
        push!(args, permute_idx(arg,j,n,acts_on(arg)))
    end
    return op.f(args...)
end
function permute_idx(op::BasicOperator,j,n,aon)
    if acts_on(op) in aon
        return permute_idx(op,j,n)
    else
        return op
    end
end

permute_idx(x::Number,args...) = x
