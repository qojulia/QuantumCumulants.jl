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

        return (true,op[a.index]*δ + dont_simplify(ex))
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
    elseif all(isa.(args,IndexedOperator))
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



# Sometimes it is necessary to avoid simplification
dont_simplify(ex::Prod) = Expression(dont_simplify, ex.args)
dont_simplify(x::Number) = x
const DontSimplify{ARGS} = Expression{typeof(dont_simplify),ARGS}
dont_simplify(ex::DontSimplify) = ex

remove_dontsimplify(x::Number) = x
remove_dontsimplify(a::DontSimplify) = prod(a.args)
remove_dontsimplify(a::AbstractOperator) = a
remove_dontsimplify(ex::Expression) = ex.f(remove_dontsimplify.(ex.args)...)

==(a::DontSimplify,b::DontSimplify) = (a.args==b.args)

Base.one(ex::DontSimplify) = one(ex.args[1])
Base.zero(ex::DontSimplify) = zero(ex.args[1])
Base.iszero(ex::DontSimplify) = iszero(ex.args[1])
Base.isone(ex::DontSimplify) = isone(ex.args[1])

*(x::Number,ex::DontSimplify) = Expression(*, [x,ex])
*(ex::DontSimplify,x::Number) = x*ex
*(ex::DontSimplify,a::BasicOperator) = Expression(*,[ex,a])
*(a::BasicOperator,ex::DontSimplify) = Expression(*,[a,ex])
*(a::DontSimplify,b::DontSimplify) = Expression(*,[a,b])
*(ex::DontSimplify,a::Prod) = Expression(*,[ex,a])
*(a::Prod,ex::DontSimplify) = Expression(*,[a,ex])

function combine_add(a::DontSimplify,b::DontSimplify)
    a_args = sort_by_inds(a.args)
    b_args = sort_by_inds(b.args)
    for i=1:length(a.args)
        check, arg = combine_add(a_args[i],b_args[i])
        check || break
    end
    if check
        return (true,2*a)
    else
        return (false,a)
    end
end
function combine_prod(a::DontSimplify,b::DontSimplify)
    @assert length(a.args)==length(b.args)==2
    for i=1:length(a.args), j=1:length(b.args)
        if a.args[i].index==b.args[j].index
            check, arg = combine_prod(a.args[i],b.args[j])
            if check
                iszero(arg) && return (true,zero(a))
                return (true,arg*dont_simplify(a.args[1+mod(i,2)]*b.args[1+mod(j,2)]))
            end
        end
    end
    return (false,a)
end
function combine_prod(a::DontSimplify,b::IndexedOperator{<:Transition})
    i = findfirst(x->x.index==b.index,a.args)
    if isa(i,Nothing)
        return (false,a)
    else
        check, arg = combine_prod(a.args[i],b)
        if check
            iszero(arg) && return (true,zero(a))
            return (true, dont_simplify(a.args[1+mod(i,2)]*arg))
        end
    end
    return (false,a)
end
function combine_prod(a::IndexedOperator{<:Transition},b::DontSimplify)
    i = findfirst(x->x.index==a.index,b.args)
    if isa(i,Nothing)
        return (false,a)
    else
        check, arg = combine_prod(a,b.args[i])
        if check
            iszero(arg) && return (true,zero(a))
            return (true, dont_simplify(arg*b.args[1+mod(i,2)]))
        end
    end
    return (false,a)
end
