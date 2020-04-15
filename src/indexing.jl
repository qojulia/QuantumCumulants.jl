import SymPy
import SymPy: sympify
import Base: ==, *


mutable struct Index{I,L,U,ID}
    label::I
    lower::L
    upper::U
    id::ID
end

const IndexOrder = Index[]
function Index(label,lower,upper)
    id = hash((label,lower,upper))
    idx = Index(label,lower,upper,id)
    order_ind = findfirst(isequal(idx),IndexOrder)
    if isa(order_ind,Nothing)
        push!(IndexOrder,idx)
    end

    # Sort by equal ids -- needed to later combine indices such as i and i!=j
    sort!(IndexOrder, by=(x->x.id))

    return idx
end

==(i::T,j::T) where T<:Index = (i.id==j.id)
==(i::Index,j::Index) = false
const SymbolicIndex{L,U,ID} = Index{<:Symbol,L,U,ID}
Base.copy(i::Index) = Index(i.label,i.lower,i.upper,copy(i.id))

# Indexing of symbolic variables
function Base.getindex(s::SymPy.Sym, inds::Index...)
    inds_ = sympify.(inds)
    b = SymPy.sympy.IndexedBase(s, real=isreal(s))
    return b[inds_...]
end
# function Base.getindex(s::SymPy.Sym, i::Int)
#     i_ = sympify.(i)
#     b = SymPy.sympy.IndexedBase(s, real=isreal(s))
#     return b[i_]
# end
function sympify(i::SymbolicIndex)
    u = i.upper
    N = isa(u,Symbol) ? SymPy.symbols(string(u),integer=true) : u
    label = string(i.label)
    return SymPy.sympy.Idx(label,(i.lower,N))
end


mutable struct IndexedOperator{OP<:BasicOperator,I} <: BasicOperator
    operator::OP
    index::I
end
Base.getindex(op::BasicOperator,i::Union{Index,Int}) = IndexedOperator(op,i)
Base.getindex(a::Union{Identity,Zero},i::Union{Index,Int}) = a
Base.adjoint(a::IndexedOperator) = IndexedOperator(adjoint(a.operator),a.index)
==(a::IndexedOperator,b::IndexedOperator) = (a.operator==b.operator && a.index==b.index)
Base.copy(a::T) where T<:IndexedOperator = T(copy(a.operator),copy(a.index))

ishermitian(a::IndexedOperator) = ishermitian(a.operator)
Base.iszero(x::IndexedOperator) = iszero(x.operator)
Base.isone(x::IndexedOperator) = isone(x.operator)

# Base.getindex(ex::DontSimplify, args...) = dont_simplify(getindex(ex.args[1], args...))

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

function combine_prod(a::IndexedOperator{<:Transition},b::IndexedOperator{<:Transition})
    if a.index==b.index
        op = a.operator*b.operator
        return (true,expression2index(op,a.index))
    else
        op = a.operator*b.operator
        δ = KroneckerDelta(sympify(a.index),sympify(b.index))

        args = sort_by_inds([a,b])
        ex = Expression(*,args)

        return (true,op[a.index]*δ + (1-δ)*dont_simplify(ex))
    end
end


function combine_prod(a::IndexedOperator,b::IndexedOperator)
    (a.index==b.index) || (false,a)
    check, op = combine_prod(a.operator,b.operator)
    return (check,expression2index(op,a.index))
end

function combine_add(a::IndexedOperator,b::IndexedOperator)
    if a==b
        return (true,2*a)
    # elseif a.operator==b.operator
    #     δ = KroneckerDelta(sympify(a.index),sympify(b.index))
    #     return (true,2*a.operator[a.index]*δ + (1-δ)*a*b) # TODO: avoid infinite recursion
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
    a_ = remove_dontsimplify(a)
    b_ = remove_dontsimplify(b)
    length(a_.args) == length(b_.args) || return false
    a_args = sort_by_inds(a_.args)
    b_args = sort_by_inds(b_.args)
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
    return SymPy.sympy.functions.special.tensor_functions.KroneckerDelta(sympify(i),sympify(j))
end
function KroneckerDelta(i::SymPy.Sym,j::SymPy.Sym)
    @assert classname(i)==classname(j)=="Idx"
    return SymPy.sympy.functions.special.tensor_functions.KroneckerDelta(i,j)
end
