# TODO: error on ProductSpace
struct ClusterSpace{H<:HilbertSpace,NType,M} <: HilbertSpace
    original_space::H
    N::NType
    order::M
end
Base.:(==)(h1::T,h2::T) where T<:ClusterSpace = (h1.original_space==h2.original_space && h1.N==h2.N && h1.order==h2.order)
Base.hash(c::ClusterSpace, h::UInt) = hash(c.original_space, hash(c.N, hash(c.order, h)))

has_hilbert(::Type{T},::ClusterSpace{<:T},args...) where T<:HilbertSpace = true

for f in [:levels,:ground_state]
    @eval $(f)(c::ClusterSpace{<:NLevelSpace}, args...) = $(f)(c.original_space, args...)
end

has_cluster(::HilbertSpace) = false
has_cluster(::ClusterSpace) = true
function has_cluster(h::ProductSpace)
    for space in h.spaces
        has_cluster(space) && return true
    end
    return false
end
has_cluster(op::AbstractOperator) = has_cluster(hilbert(op))

struct ClusterAon{T<:Int}
    i::T
    j::T
end
Base.hash(c::T, h::UInt) where T<:ClusterAon = hash(T, hash(c.i, hash(c.j, h)))
Base.getindex(v::Vector{<:HilbertSpace}, i::ClusterAon) = v[i.i]
Base.length(::ClusterAon) = 1

Base.isequal(c1::T,c2::T) where T<:ClusterAon = (c1.i==c2.i && c1.j==c2.j)
Base.isless(i::Int,c::ClusterAon) = isless(i,c.i)
Base.isless(c::ClusterAon,i::Int) = isless(c.i,i)
function Base.isless(c1::ClusterAon,c2::ClusterAon)
    if isless(c1.i, c2.i)
        return true
    else
        return isless(c1.j, c2.j)
    end
end
Base.iterate(c::ClusterAon, state=1) = isone(state) ? (c,state+1) : nothing

const AonType = Union{Int,ClusterAon{Int}}

function Transition(h::ClusterSpace{<:NLevelSpace}, name, i, j, aon::Int=1)
    op = Transition(hilbert.original_space,name,i,j,1)
    return _cluster(h,op,aon)
end
function Transition(hilbert::ProductSpace,name,i,j)
    inds = findall(x->isa(x,NLevelSpace) || isa(x,ClusterSpace{<:NLevelSpace}),hilbert.spaces)
    if length(inds)==1
        if hilbert.spaces[inds] isa ClusterSpace
            op = Transition(hilbert.original_space,name,i,j,1)
            return _cluster(hilbert, op, inds[1])
        else
            return Transition(hilbert,name,i,j,inds[1])
        end
    else
        isempty(inds) && error("Can only create Transition on NLevelSpace! Not included in $(hilbert)")
        length(inds)>1 && error("More than one NLevelSpace in $(hilbert)! Specify on which Hilbert space Transition should be created with Transition(hilbert,name,i,j,acts_on)!")
    end
end
function Transition(hilbert::H,name::S,i::I,j::I,aon::Int) where {H<:ProductSpace,S,I}
    if hilbert.spaces[aon] isa ClusterSpace
        op = Transition(hilbert.spaces[aon].original_space,name,i,j,1)
        return _cluster(hilbert, op, aon)
    else
        return Transition(hilbert,name,i,j,aon)
    end
end

# TODO: Destroy, Create

function _cluster(h::ProductSpace, op::BasicOperator, aon::Int)
    order = h.spaces[aon].order
    return _cluster(h, op, aon, order)
end
function _cluster(h::ClusterSpace, op::BasicOperator, aon::Int)
    return _cluster(h, op, aon, h.order)
end
function _cluster(h, op, aon, order)
    ops = BasicOperator[]
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

# expand_clusters(c::Cluster) = sum(c.ops)
# expand_clusters(x) = x
# function expand_clusters(t::OperatorTerm)
#     args = AbstractOperator[]
#     for arg in t.arguments
#         push!(args, expand_clusters(t.arguments))
#     end
#     return t.f(args...)
# end
