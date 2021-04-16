## ClusterSpace methods
for f in [:levels,:ground_state]
    @eval $(f)(c::ClusterSpace{<:NLevelSpace}, args...) = $(f)(c.original_space, args...)
end

has_cluster(::HilbertSpace,args...) = false
has_cluster(::ClusterSpace,args...) = true
function has_cluster(h::ProductSpace)
    for space in h.spaces
        has_cluster(space) && return true
    end
    return false
end
has_cluster(h::ProductSpace,aon) = has_cluster(h.spaces[get_i(aon)])
has_cluster(op::QNumber,args...) = has_cluster(hilbert(op),args...)
has_cluster(avg::Average,args...) = has_cluster(undo_average(avg),args...)

# ClusterAon methods
Base.hash(c::T, h::UInt) where T<:ClusterAon = hash(T, hash(c.i, hash(c.j, h)))
Base.getindex(v::Vector{<:HilbertSpace}, c::ClusterAon) = v[c.i]
get_i(x::Integer) = x
get_j(x::Integer) = x
get_i(x::ClusterAon) = x.i
get_j(x::ClusterAon) = x.j
Base.length(::ClusterAon) = 1

extract_names(names::Vector, i::Int) = names[i]
extract_names(names::Vector, c::ClusterAon) = names[c.i][c.j]
function extract_names(names::Vector, v::Vector)
    [extract_names(names, v_) for v_ in v]
end

Base.isequal(c1::T,c2::T) where T<:ClusterAon = (c1.i==c2.i && c1.j==c2.j)
Base.isless(i::Int,c::ClusterAon) = isless(i,c.i)
Base.isless(c::ClusterAon,i::Int) = isless(c.i,i)
function Base.isless(c1::ClusterAon,c2::ClusterAon)
    if isless(c1.i, c2.i)
        return true
    elseif isequal(c1.i, c2.i)
        return isless(c1.j, c2.j)
    else
        return false
    end
end
Base.iterate(c::ClusterAon, state=1) = isone(state) ? (c,state+1) : nothing

function _cluster(h::ProductSpace, op::QSym, aon::Int)
    order = h.spaces[aon].order
    return _cluster(h, op, aon, order)
end
function _cluster(h::ClusterSpace, op::QSym, aon::Int)
    return _cluster(h, op, aon, h.order)
end
function _cluster(h, op, aon, order)
    ops = QSym[]
    for i=1:order
        name = Symbol(op.name, :_, i)
        aon_i = ClusterAon(aon[1],i)
        op_ = _remake_op(op, h, name, aon_i)
        push!(ops, op_)
    end
    return ops
end

_remake_op(op::Transition, h, name, aon) = Transition(h, name, op.i, op.j, aon)
_remake_op(op::Destroy, h, name, aon) = Destroy(h, name, aon)
_remake_op(op::Create, h, name, aon) = Create(h, name, aon)
