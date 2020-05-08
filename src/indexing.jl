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
function Base.getindex(s::SymPy.Sym, inds...)
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

        return (true,op[a.index]*δ + neq_inds_prod(ex))
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
    return all(as .== bs)
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
neq_inds_prod(x::Number) = x
const NeqIndsProd{ARGS} = Expression{typeof(neq_inds_prod),ARGS}
neq_inds_prod(ex::NeqIndsProd) = ex

function neq_inds_prod(args::Vector)
    args_ = []
    cs = Number[]
    for a=args
        if isa(a,NeqIndsProd)
            append!(args_,a.args)
        elseif isa(a,Number)
            push!(cs,a)
        elseif isa(a,BasicOperator)
            push!(args_,a)
        else
            error("Something went wrong here!")
        end
    end
    if isempty(cs)
        return Expression(neq_inds_prod, sort_by_inds(args_))
    else
        return prod(cs)*Expression(neq_inds_prod, sort_by_inds(args_))
    end
end
neq_inds_prod(ex::Prod) = neq_inds_prod(ex.args)
neq_inds_prod(ex::Add) = sum(neq_inds_prod.(ex.args))

remove_NeqIndsProd(x::Number) = x
remove_NeqIndsProd(a::NeqIndsProd) = prod(a.args)
remove_NeqIndsProd(a::AbstractOperator) = a
remove_NeqIndsProd(ex::Expression) = ex.f(remove_NeqIndsProd.(ex.args)...)

==(a::NeqIndsProd,b::NeqIndsProd) = (a.args==b.args)

Base.one(ex::NeqIndsProd) = one(ex.args[1])
Base.zero(ex::NeqIndsProd) = zero(ex.args[1])
Base.iszero(ex::NeqIndsProd) = any(iszero.(ex.args))
Base.isone(ex::NeqIndsProd) = all(isone.(ex.args))

*(x::Number,ex::NeqIndsProd) = Expression(*, [x,ex])
*(ex::NeqIndsProd,x::Number) = x*ex
*(ex::NeqIndsProd,a::BasicOperator) = Expression(*,[ex,a])
*(a::BasicOperator,ex::NeqIndsProd) = Expression(*,[a,ex])
*(a::NeqIndsProd,b::NeqIndsProd) = Expression(*,[a,b])
*(ex::NeqIndsProd,a::Prod) = Expression(*,[ex;a.args])
*(a::Prod,ex::NeqIndsProd) = Expression(*,[a.args;ex])

function combine_add(a::NeqIndsProd,b::NeqIndsProd)
    length(a.args)==length(b.args) || return (false,a)
    a_args = sort_by_inds(a.args)
    b_args = sort_by_inds(b.args)
    check = false
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
function combine_prod(a::NeqIndsProd,b::NeqIndsProd)
    # Try to combine any equal indices
    args_out = []
    checks_a = ones(Bool, length(a.args))
    checks_b = ones(Bool, length(b.args))
    for i=1:length(a.args), j=1:length(b.args)
        if a.args[i].index==b.args[j].index
            check, arg = combine_prod(a.args[i],b.args[j])
            if check
                checks_a[i] = false
                checks_b[j] = false
                iszero(arg) && return (true,zero(a))
                push!(args_out, arg)
            end
        end
    end

    if !all(checks_a) || !all(checks_b)
        append!(args_out, a.args[checks_a])
        append!(args_out, b.args[checks_b])
        return (true,neq_inds_prod(args_out))
    end

    # TODO: Combine elements if no equal indices occur
    # arg_a = copy(a.args[end])
    # arg_b = copy(b.args[1])
    # check = false
    # args_out = [AbstractOperator[]]
    # for i=0:length(a.args)-1
    #     check, arg_a = combine_prod(a.args[end-i],arg_b)
    #     check || break
    #     iszero(arg) && return (true,zero(a))
    #     for j=1:length(b.args)
    #         check, arg_ = combine_prod(arg_a,b.args[j])
    #         check || break
    #         iszero(arg_) && return (true,zero(b))
    #         push!(args_out[i+1], arg_)
    #     end
    #     check || break
    #     arg_b = neq_inds_prod((args_out[i+1]))
    # end
    # if check
    #     reverse!(args_out)
    #     out = neq_inds_prod(prod(prod.(args_out)))
    #     return (true,out)
    # end

    return (false,a)
end
function combine_prod(a::NeqIndsProd,b::IndexedOperator{<:Transition})
    # Try to combine equal indices
    inds = [a1.index for a1=a.args]
    i = findfirst(isequal(b.index), inds)
    if !isa(i,Nothing)
        check, arg = combine_prod(a.args[i],b)
        check || return (false,a)
        iszero(arg) && return (true,zero(b))
        args = sort_by_inds([a.args[1:i-1];arg;a.args[i+1:end]])
        out = prod(args)
        return (true,neq_inds_prod(out))
    end

    #Otherwise step through and combine other indices
    check, arg = combine_prod(a.args[end],b)
    if check
        iszero(arg) && return (true,zero(b))
        for i=1:length(a.args)-1
            check, arg = _combine_neq_inds_prod_rtl(a.args[end-i],arg,inds[end-1+i])
            check || break
        end
        check && return (true,arg)
    end

    return (false,a)
end
function combine_prod(a::IndexedOperator{<:Transition},b::NeqIndsProd)
    # Try to combine equal indices
    inds = [b1.index for b1=b.args]
    i = findfirst(isequal(a.index),inds)
    if !isa(i,Nothing)
        check, arg = combine_prod(a,b.args[i])
        check || return (false,a)
        iszero(arg) && return (true,zero(a))
        args = sort_by_inds([b.args[1:i-1];arg;b.args[i+1:end]])
        out = prod(args)
        return (true,neq_inds_prod(out))
    end

    #Otherwise step through and combine other indices
    check, arg = combine_prod(a,b.args[1])
    if check
        iszero(arg) && return (true,zero(a))
        for i=2:length(b.args)
            check, arg = _combine_neq_inds_prod_ltr(arg,b.args[i],inds[i-1])
            check || break
        end
        check && return (true,arg)
    end

    return (false,a)
end
function _combine_neq_inds_prod_rtl(a_arg,arg,index)
    @assert isa(arg,Add)

    # Combine the second argument
    args_ = filter(x->x.index!=index, arg.args[2].args)
    @assert length(args_)==1
    check, arg_ = combine_prod(a_arg,args_[1])
    @assert check
    args_remain = filter(x->x.index==index, arg.args[2].args)
    arg2 = neq_inds_prod(prod([arg_; args_remain]))

    # Add the first argument multiplied with a_arg if non-zero
    if iszero(arg.args[1])
        return (true, arg2)
    else
        arg1 = neq_inds_prod(prod([a_arg,arg.args[1]]))
        return (true, arg1+arg2)
    end
end
function _combine_neq_inds_prod_ltr(arg,b_arg,index)
    @assert isa(arg,Add)

    # Combine the second argument
    args_ = filter(x->x.index!=index, arg.args[2].args)
    @assert length(args_)==1
    check, arg_ = combine_prod(args_[1],b_arg)
    @assert check
    args_remain = filter(x->x.index==index, arg.args[2].args)
    arg2 = neq_inds_prod(prod([arg_; args_remain]))

    # Add the first argument multiplied with a_arg if non-zero
    if iszero(arg.args[1])
        return (true, arg2)
    else
        arg1 = neq_inds_prod(prod([arg.args[1],b_arg]))
        return (true, arg1+arg2)
    end
end
