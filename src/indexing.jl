import SymPy
import SymPy: sympify
import Base: ==, *

# const KroneckerDelta = SymPy.sympy.functions.special.tensor_functions.KroneckerDelta

mutable struct Index{I,L,U,ID,NID}
    label::I
    lower::L
    upper::U
    id::ID
    nid::NID
    full_id::ID
end

const IndexOrder = Index[]
function Index(label,lower,upper;neq::Vector{<:Index}=Index[])
    nid = UInt[i.id for i=unique(neq)]
    return Index(label,lower,upper,nid)
end
function Index(label,lower,upper,nid_::Vector{<:UInt})
    nid = sort(nid_)
    id = hash((label,lower,upper))
    full_id = hash((label,lower,upper,nid))
    idx = Index(label,lower,upper,id,nid,full_id)
    order_ind = findfirst(isequal(idx),IndexOrder)
    if isa(order_ind,Nothing)
        push!(IndexOrder,idx)
    end

    # Sort by equal ids -- needed to later combine indices such as i and i!=j
    sort!(IndexOrder, by=(x->x.id))

    return idx
end

==(i::T,j::T) where T<:Index = (i.full_id==j.full_id)
==(i::Index,j::Index) = false
const SymbolicIndex{L,U,ID,NID} = Index{<:Symbol,L,U,ID,NID}
Base.copy(i::Index) = Index(i.label,i.lower,i.upper,copy(i.id),copy(i.nid),copy(i.full_id))

# Indexing of symbolic variables
function Base.getindex(s::SymPy.Sym, inds...) #TODO: Type constraints
    inds_ = sympify.(inds)
    b = SymPy.sympy.IndexedBase(s, real=isreal(s))
    return b[inds_...]
end
function Base.getindex(s::SymPy.Sym, i::Int) #TODO: Type constraints
    i_ = sympify.(i)
    b = SymPy.sympy.IndexedBase(s, real=isreal(s))
    return b[i_]
end
function sympify(i::SymbolicIndex)
    u = i.upper
    N = isa(u,Symbol) ? SymPy.symbols(string(u),integer=true) : u
    label = string(i.label)
    if length(i.nid)==1
        j = IndexOrder[findfirst(x->x.id∈i.nid,IndexOrder)]
        label *= "≠"*string(j.label)
    elseif length(i.nid) > 1
        inds = IndexOrder[findall(x->x.id∈i.nid,IndexOrder)]
        unique!(x->x.id, inds)
        label *= "≠{"
        for j=inds
            label *= string(j.label)
        end
        label *= "}"
    end
    return SymPy.sympy.Idx(label,(i.lower,N))
end


mutable struct IndexedOperator{OP<:BasicOperator,I} <: BasicOperator
    operator::OP
    index::I
end
Base.getindex(op::BasicOperator,i::Union{Index,Int}) = IndexedOperator(op,i)
Base.getindex(a::Union{Identity,Zero},i::Union{Index,Int}) = a
Base.adjoint(a::IndexedOperator) = IndexedOperator(adjoint(a.operator),a.index)
==(a::IndexedOperator,b::IndexedOperator) = (a.operator==b.operator && a.index==b.index) #TODO: should this be a strict equality check for index?
Base.copy(a::T) where T<:IndexedOperator = T(copy(a.operator),copy(a.index))

ishermitian(a::IndexedOperator) = ishermitian(a.operator)
Base.iszero(x::IndexedOperator) = iszero(x.operator)
Base.isone(x::IndexedOperator) = isone(x.operator)


replace_commutator(a::IndexedOperator,::IndexedOperator) = (false,a)
function replace_commutator(a::IndexedOperator{<:Destroy},b::IndexedOperator{<:Create})
    check, op = replace_commutator(a.operator,b.operator)
    if a.index==b.index
        return (check,expression2index(op,a.index))
    else
        δ = KroneckerDelta(sympify(a.index),sympify(b.index))
        ex = expression2index(op, a.index)
        return (check, b*a+δ*ex)
    end
end

function *(a::IndexedOperator{<:Transition},b::IndexedOperator{<:Transition})
    op = a.operator*b.operator
    if a.index==b.index
        return expression2index(op,a.index)
    elseif a.index.id==b.index.id # Don't share the same nid
        nid = [a.index.nid;b.index.nid]
        unique!(nid)
        index = Index(a.label,a.lower,a.upper,nid)
        return expression2index(op,index)
    elseif a.index.id ∈ b.index.nid || b.index.id ∈ a.index.nid
        args = sort_by_inds([a,b])
        return Expression(*,args)
    else
        δ = KroneckerDelta(sympify(a.index),sympify(b.index))

        # Generate new index with non-equal index a.index for b
        i_nid = [a.index.nid;b.index.id]
        unique!(i_nid)
        i = Index(a.index.label,a.index.lower,a.index.upper,i_nid)
        j_nid = [b.index.nid;a.index.id]
        unique!(j_nid)
        j = Index(b.index.label,b.index.lower,b.index.upper,j_nid)

        # Generate new (a,b)
        a_ = a.operator[i]
        b_ = b.operator[j]
        args = sort_by_inds([a_,b_])
        ex = Expression(*,args)

        return op[a.index]*δ + ex
    end
end

function combine_prod(a::IndexedOperator{<:Transition},b::IndexedOperator{<:Transition})
    if a.index==b.index
        op = a.operator*b.operator
        return (true,expression2index(op,a.index))
    elseif a.index.id ∈ b.index.nid || b.index.id ∈ a.index.nid
        args = sort_by_inds([a,b])
        if args==[a,b]
            return (false,a)
        else
            return (true,Expression(*,args))
        end
    elseif a.index.id==b.index.id # Don't share the same nid
        index = if length(a.index.nid) > length(b.index.nid)
            a.index
        elseif length(a.index.nid) < length(b.index.nid)
            b.index
        else
            nid = [a.index.nid;b.index.nid]
            unique!(nid)
            Index(a.index.label,a.index.lower,a.index.upper,nid)
        end
        op = a.operator*b.operator
        return (true,expression2index(op,index))
    else
        op = a.operator*b.operator
        δ = KroneckerDelta(sympify(a.index),sympify(b.index))

        # Generate new index with non-equal index a.index for b
        i_nid = [a.index.nid;b.index.id]
        unique!(i_nid)
        i = Index(a.index.label,a.index.lower,a.index.upper,i_nid)
        j_nid = [b.index.nid;a.index.id]
        unique!(j_nid)
        j = Index(b.index.label,b.index.lower,b.index.upper,j_nid)

        # Generate new (a,b)
        a_ = a.operator[i]
        b_ = b.operator[j]
        args = sort_by_inds([a_,b_])
        ex = Expression(*,args)

        return (true,op[a.index]*δ + ex)
    end
end


function combine_prod(a::IndexedOperator,b::IndexedOperator)
    if a.index==b.index
        check, op = combine_prod(a.operator,b.operator)
        return (check,expression2index(op,a.index))
    elseif a.index.id==b.index.id # Don't share same .nid
        nid = [a.index.nid;b.index.nid]
        unique!(nid)
        i = Index(a.index.label,a.index.lower,a.index.upper,nid)
        j = Index(b.index.label,b.index.lower,b.index.upper,nid)
        return combine_prod(a.operator[i],b.operator[j])
    else
        return (false,a)
    end
end

function combine_add(a::IndexedOperator,b::IndexedOperator)
    if a==b
        return (true,2*a)
    elseif a.operator==b.operator && a.index.id==b.index.id
        nid = [a.index.nid;b.index.nid]
        unique!(nid)
        i = Index(a.index.label,a.index.lower,a.index.upper,nid)
        ex = 2*IndexedOperator(a.operator,i)

        a_nid = filter(x->!(x∈a.index.nid), nid)
        if !isempty(a_nid)
            js = IndexOrder[findall(x->x.id∈a_nid,IndexOrder)]
            filter!(x->isempty(x.nid),js)
            for j=js
                δ = KroneckerDelta(a.index,j)
                iszero(δ) && continue
                a_ = δ*IndexedOperator(a.operator, j)
                ex += a_
            end
        end

        b_nid = filter(x->!(x∈b.index.nid), nid)
        if !isempty(b_nid)
            ks = IndexOrder[findall(x->x.id∈b_nid,IndexOrder)]
            filter!(x->isempty(x.nid),ks)
            for k=ks
                δ = KroneckerDelta(b.index, k)
                iszero(δ) && continue
                b_ = δ*IndexedOperator(b.operator, k)
                ex += b_
            end
        end

        return (true,ex)
    else
        return (false,a)
    end
end

function expression2index(ex::Expression,i)
    args_ = [expression2index(arg,i) for arg=ex.args]
    return ex.f(args_...)
end
expression2index(a::BasicOperator,i) = getindex(a,i)
expression2index(a::IndexedOperator,i...) = error()
expression2index(x::Number,i...) = x

# TODO: optimize for non-indexed operators
function ==(a::Prod,b::Prod)
    length(a.args) == length(b.args) || return false
    a_args = sort_by_inds(a.args)
    b_args = sort_by_inds(b.args)
    return all(a_args .== b_args)
end
function sort_by_inds(args::Vector)
    args_num = filter(x->isa(x,Number),args)
    if !isempty(args_num)
        args_remain = filter(x->!isa(x,Number),args)
        return [args_num;sort_by_inds(args_remain)]
    elseif isa(args[1],IndexedOperator)
        inds = Int[]
        for a1=args
            if (iszero(a1)||isone(a1))
                i = 1
            else
                i = a1.index
            end
            if isa(i,Int)
                push!(inds,i)
            else
                push!(inds,findfirst(isequal(i),IndexOrder))
            end
        end
        p = sortperm(inds)
        return args[p]
    else
        return args
    end
end
function isequal_prod_args(as,bs)
    length(as) == length(bs) || return false
    as_ = sort_by_inds(as)
    bs_ = sort_by_inds(bs)
    return all(as_ .== bs_)
end

function acts_on(op::IndexedOperator)
    # TODO: change acts_on for TensorProd to account for IndexedOperators
    i = findfirst(isequal(op.index), IndexOrder)
    isa(i, Nothing) && error("Something went wrong here!")
    return [i]
end

function simplify_operators(op::IndexedOperator)
    op_ = simplify_operators(op.operator)
    out = expression2index(op_,op.index)
    if out==op
        return op
    else
        return simplify_operators(out)
    end
end


function KroneckerDelta(i::Index,j::Index)
    (i.id∈j.nid || j.id∈i.nid) && return 0
    return SymPy.sympy.functions.special.tensor_functions.KroneckerDelta(sympify(i),sympify(j))
end
function KroneckerDelta(i::SymPy.Sym,j::SymPy.Sym)
    @assert classname(i)==classname(j)=="Idx"
    i_ = IndexOrder[findfirst(x->sympify(x)==i,IndexOrder)]
    j_ = IndexOrder[findfirst(x->sympify(x)==j,IndexOrder)]
    (i_.id∈j_.nid || j_.id∈i_.nid) && return 0
    return SymPy.sympy.functions.special.tensor_functions.KroneckerDelta(i,j)
end
