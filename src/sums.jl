import SymPy.Introspection: classname

mutable struct Sum{I} <: Function
    index::I
end
==(s1::Sum,s2::Sum) = (s1.index==s2.index)
const SumType{argType} = Expression{<:Sum,argType}

function Sum(op::IndexedOperator,i::Index)
    if op.index==i
        return Expression(Sum(i),[op])
    elseif op.index.id==i.id # Don't have the same .nid
        nid = [op.index.nid;i.nid]
        unique!(nid)
        i_ = Index(i.label,i.lower,i.upper,nid)
        op_ = IndexedOperator(op.operator,i_)
        return Sum(op_, i_)
    else
        u = i.upper
        N = isa(u,Symbol) ? SymPy.symbols(string(u),integer=true) : u
        return (N-length(i.nid))*op
    end
end
Sum(op::IndexedOperator) = Expression(Sum(op.index),[op])
function Sum(op::BasicOperator,i::Index)
    u = i.upper
    N = isa(u,Symbol) ? SymPy.symbols(string(u),integer=true) : u
    return (N-length(i.nid))*op
end
function Sum(ex::Prod,i::Index)
    sym_inds, op_inds = find_index(ex)
    if i∈op_inds || i∈sym_inds
        return Expression(Sum(i),[ex])
    elseif id_in(i,op_inds) || id_in(i,sym_inds)
        inds = filter(x->x.id==i.id,op_inds)
        append!(inds,filter(x->x.id==i.id,sym_inds))
        unique!(inds)
        nid = UInt[]
        for j=inds
            append!(nid, j.nid)
        end
        unique!(nid)
        i_ = Index(i.label,i.lower,i.upper,nid)
        @assert (i_∈op_inds || i_∈sym_inds) # Otherwise causes infinite loop
        return Sum(ex, i_)
    else
        u = i.upper
        N = isa(u,Symbol) ? SymPy.symbols(string(u),integer=true) : u
        return (N-length(i.nid))*ex
    end
end
function Sum(ex::TensorProd,i::Index)
    # ex = simplify_operators(ex_)
    sym_inds, op_inds = find_index(ex)
    if i∈op_inds || i∈sym_inds
        return Expression(Sum(i),[ex])
    elseif id_in(i,op_inds) || id_in(i,sym_inds)
        inds = filter(x->x.id==i.id,op_inds)
        append!(inds,filter(x->x.id==i.id,sym_inds))
        unique!(inds)
        nid = UInt[]
        for j=inds
            append!(nid, j.nid)
        end
        unique!(nid)
        i_ = Index(i.label,i.lower,i.upper,nid)
        @assert (i_∈op_inds || i_∈sym_inds) # Otherwise causes infinite loop
        return Sum(ex, i_)
    else
        u = i.upper
        N = isa(u,Symbol) ? SymPy.symbols(string(u),integer=true) : u
        return (N-length(i.nid))*ex
    end
end
Sum(ex::Add,i::Index) = sum([Sum(a1,i) for a1=ex.args])

function id_in(i::Index,inds::Vector)
    ids = [j.id for j=inds]
    return i.id∈ids
end

function Sum(s::SumType,i::Index)
    sym_inds, op_inds = find_index(s)
    if i∈op_inds || i∈sym_inds
        return Expression(Sum(i),[s])
    elseif id_in(i,op_inds) || id_in(i,sym_inds)
        inds = filter(x->x.id==i.id,op_inds)
        append!(inds,filter(x->x.id==i.id,sym_inds))
        unique!(inds)
        nid = UInt[]
        for j=inds
            append!(nid, j.nid)
        end
        unique!(nid)
        i_ = Index(i.label,i.lower,i.upper,nid)
        @assert (i_∈op_inds || i_∈sym_inds) # Otherwise causes infinite loop
        return Sum(s, i_)
    else
        u = i.upper
        N = isa(u,Symbol) ? SymPy.symbols(string(u),integer=true) : u
        return (N-length(i.nid))*s
    end
end

function find_index(ex::Expression)
    c_inds = Index[]
    op_inds = Index[]
    for arg=ex.args
        c_, op_ = find_index(arg)
        append!(c_inds, c_)
        append!(op_inds, op_)
    end
    unique!(x->x.id,c_inds)
    unique!(x->x.id,op_inds)
    return c_inds, op_inds
end
find_index(op::BasicOperator) = (Index[],Index[])
find_index(op::IndexedOperator) = (Index[],[op.index])
find_index(::Number) = (Index[],Index[])
function find_index(s::SymPy.Sym)
    cs = classname(s)
    if cs == "Indexed" || cs == "KroneckerDelta"
        sympy_inds = [s.__pyobject__.indices...]
        inds = [findfirst(x->sympify(x)==i,IndexOrder) for i=sympy_inds]
        c_inds = IndexOrder[inds]
    elseif cs=="Mul"
        c, args = s.__pyobject__.as_coeff_mul()
        c_inds = Index[]
        for arg=args
            c_, _ = find_index(arg)
            append!(c_inds, c_)
        end
    elseif cs=="Add"
        c, args = s.__pyobject__.as_coeff_add()
        c_inds = Index[]
        for arg=args
            c_, _ = find_index(arg)
            append!(c_inds, c_)
        end
    elseif cs=="adjoint"
        c_inds = find_index(s')[1]
    else
        c_inds = Index[]
    end
    unique!(x->x.id,c_inds)
    return c_inds,Index[]
end

function gen_index(i::SymPy.Sym)
    return IndexOrder[findfirst(x->sympify(x)==i,IndexOrder)]
end

Base.iszero(s::SumType) = iszero(s.args[1])
Base.zero(s::SumType) = zero(s.args[1])
Base.one(s::SumType) = one(s.args[1])

# Algebra
# +
+(s1::SumType,s2::SumType) = Expression(+,[s1,s2])
-(s::SumType) = Sum(-s.args[1],s.f.index)

# *
*(s1::SumType,s2::SumType) = Sum(Sum(s1.args[1]*s2.args[1],s1.f.index),s2.f.index)
function *(op::IndexedOperator,s::SumType)
    (op.index.lower==s.f.index.lower && op.index.upper==s.f.index.upper) || error("Something went wrong here!")
    if op.index==s.f.index
        i = gen_index(op.index.lower,op.index.upper,op.index.nid)
        sum_arg = swap_index(s.args[1],s.f.index,i)
        return op*Sum(sum_arg,i)
    elseif op.index.id==s.f.index.id # Don't have same .nid
        nid = [op.index.nid;s.f.index.nid]
        unique!(nid)
        i = gen_index(op.index.lower,op.index.upper,nid)
        sum_arg = swap_index(s.args[1],s.f.index,i)
        return op*Sum(sum_arg,i)
    else
        return Sum(op*s.args[1],s.f.index)
    end
end
function *(s::SumType,op::IndexedOperator)
    (op.index.lower==s.f.index.lower && op.index.upper==s.f.index.upper) || error("Something went wrong here!")
    if op.index==s.f.index
        i = gen_index(op.index.lower,op.index.upper,op.index.nid)
        sum_arg = swap_index(s.args[1],s.f.index,i)
        return Sum(sum_arg,i)*op
    elseif op.index.id==s.f.index.id # Don't have same .nid
        nid = [op.index.nid;s.f.index.nid]
        unique!(nid)
        i = gen_index(op.index.lower,op.index.upper,nid)
        sum_arg = swap_index(s.args[1],s.f.index,i)
        return Sum(sum_arg,i)*op
    else
        return Sum(s.args[1]*op,s.f.index)
    end
end
*(s::SumType,op::BasicOperator) = Sum(s.args[1]*op,s.f.index)
*(op::BasicOperator,s::SumType) = Sum(op*s.args[1],s.f.index)
function *(s::SumType,op::Union{Prod,TensorProd})
    sym_inds, ops_inds = find_index(op)
    inds = [sym_inds;ops_inds]
    if s.f.index ∈ inds
        i = gen_index(s.f.index.lower,s.f.index.upper,s.f.index.nid)
        sum_arg = swap_index(s.args[1],s.f.index,i)
        return Sum(sum_arg,i)*op
    elseif id_in(s.f.index, inds)
        nid = [s.f.index.nid;i.nid]
        i = gen_index(s.f.index.lower,s.f.index.upper,nid)
        sum_arg = swap_index(s.args[1],s.f.index,i)
        return Sum(sum_arg,i)*op
    else
        return Sum(s.args[1]*op,s.f.index)
    end
end
function *(op::Union{Prod,TensorProd},s::SumType)
    sym_inds, ops_inds = find_index(op)
    inds = [sym_inds;ops_inds]
    if s.f.index ∈ inds
        i = gen_index(s.f.index.lower,s.f.index.upper,s.f.index.nid)
        sum_arg = swap_index(s.args[1],s.f.index,i)
        return op*Sum(sum_arg,i)
    elseif id_in(s.f.index, inds)
        nid = UInt[]
        for i=inds
            append!(nid, i.nid)
        end
        unique!(nid)
        i = gen_index(s.f.index.lower,s.f.index.upper,nid)
        sum_arg = swap_index(s.args[1],s.f.index,i)
        return op*Sum(sum_arg,i)
    else
        return Sum(op*s.args[1],s.f.index)
    end
end
*(s::SumType,op::Prod) = Sum(s.args[1]*op,s.f.index)
*(op::Prod,s::SumType) = Sum(op*s.args[1],s.f.index)
# *(s::SumType,op::Add) = sum(s.args[1]*op,s.f.index)
# *(op::TensorProd,s::SumType) = Sum(op*s.args[1],s.f.index)


*(x::Number,s::SumType) = Sum(x*s.args[1],s.f.index)
gen_index(s::Symbol,lower,upper) = Index(s,lower,upper)
# TODO: generate more legible symbols for new indices
gen_index(lower,upper) = Index(gensym(:Idx),lower,upper)
gen_index(lower,upper,nid::Vector) = Index(gensym(:Idx),lower,upper,nid)

# ⊗
⊗(s1::SumType,s2::SumType) = Sum(Sum(s1.args[1]⊗s2.args[1],s1.f.index),s2.f.index)
⊗(s::SumType,op::AbstractOperator) = Sum(s.args[1]⊗op,s.f.index)
⊗(op::AbstractOperator,s::SumType) = Sum(op⊗s.args[1],s.f.index)

# Simplification
function simplify_operators(s::SumType)
    s_ = _max_nid(s)
    !isa(s_,SumType) && return simplify_operators(s_)
    s_index = s_.f.index
    op = simplify_operators(s_.args[1])
    if isa(op,Add)
        return simplify_operators(Sum(op,s_index))
    end
    iszero(op) && return zero(s)
    c_inds, op_inds = find_index(op)
    ids = [i.id for i=[c_inds;op_inds]]
    if !(s_index.id∈ids)
        u = s_index.upper
        N = isa(u,Symbol) ? SymPy.symbols(string(u),integer=true) : u
        return N*op
    else
        check, arg = resolve_kdelta(op,s_index)
        #TODO: if check - remove s_index from all .nid fields in arg
        if check
            return _remove_nid(arg)
        else
            return _remove_nid(Sum(arg,s_index))
        end
    end
end

function _max_nid(arg::AbstractOperator)
    sym_inds, op_inds = find_index(arg)
    inds = [sym_inds;op_inds]
    swaps = Tuple{Index,Index}[]
    for i=1:length(inds)
        inds_ = [inds[1:i-1];inds[i+1:end]]
        ids = [j.id for j=inds_]
        if inds[i].id ∈ ids
            nid = copy(inds[i].nid)
            for j=inds_
                if j.id==inds[i].id
                    append!(nid, j.nid)
                end
            end
            unique!(nid)
            sort!(nid)
            if length(nid) != length(inds[i].nid)
                i_ = Index(inds[i].label,inds[i].lower,inds[i].upper,nid)
                push!(swaps, (inds[i], i_))
            end
        end
    end
    arg_ = copy(arg)
    for s=swaps
        arg_ = swap_index(arg_, s...)
    end
    return arg_
end

function _remove_nid(arg::AbstractOperator)
    sym_inds, op_inds = find_index(arg)
    inds = [sym_inds;op_inds]
    ids = [i.id for i=inds]
    nid = [i.nid for i=inds]
    swaps = Tuple{Index,Index}[]
    for i=1:length(inds)
        nid_ = UInt[]
        for j=1:length(nid[i])
            if (nid[i][j] ∈ ids)
                push!(nid_, nid[i][j])
            end
        end
        if !(length(nid[i]) == length(nid_))
            i_ = Index(inds[i].label,inds[i].lower,inds[i].upper,nid_)
            push!(swaps, (inds[i], i_))
        end
    end
    arg_ = copy(arg)
    for s=swaps
        arg_ = swap_index(arg_, s...)
    end
    return arg_
end


function resolve_kdelta(ex::SumType,index)
    check, arg = resolve_kdelta(ex.args[1],index)
    !check && return (false, Sum(arg, ex.f.index))
    return (check, Sum(arg, ex.f.index))
end
function resolve_kdelta(ex::Prod,index)
    if isa(ex.args[1],SymPy.Sym)
        arg = SymPy.expand(ex.args[1])
        iszero(arg) && return (true,Zero())
        return _resolve_kdelta(ex,index,arg)
    end
    return (false,ex)
end
resolve_kdelta(op::BasicOperator,args...) = (false,op)
function resolve_kdelta(ex::TensorProd,index)
    if isa(ex.args[1],Prod)
        arg_ = ex.args[1].args[1]
        if isa(arg_,SymPy.Sym)
            arg = SymPy.expand(arg_)
            iszero(arg) && return (true,zero(ex.args[1]))
            return _resolve_kdelta(ex,index,arg)
        end
    end
    return (false,ex)
end
function _resolve_kdelta(ex::Prod,index,arg::SymPy.Sym)
    cs = classname(arg)
    if cs=="KroneckerDelta"
        i,j = arg.__pyobject__.indices
        i_ = IndexOrder[findfirst(x->sympify(x)==i,IndexOrder)]
        j_ = IndexOrder[findfirst(x->sympify(x)==j,IndexOrder)]
        if i_.id==index.id
            out = swap_index(ex,index,j_)
            return (true, simplify_operators(out))
        elseif j_.id==index.id
            out = swap_index(ex,index,i_)
            return (true, simplify_operators(out))
        else
            return (false,ex)
        end
    elseif cs=="Mul"
        c, margs = arg.__pyobject__.as_coeff_mul()
        check = false
        for m=margs
            check, ex_ = _resolve_kdelta(ex,index,m)
            check && break
        end
        return (check, ex_)
    elseif cs=="Add"
        c, aargs = arg.__pyobject__.as_coeff_add()
        ex_ = ex.f(ex.args[2:end]...)
        add_delta_args = []
        add_others = []
        for m=(c,aargs...)
            check, _ex_ = _resolve_kdelta(m*ex_,index,m)
            if check
                push!(add_delta_args, _ex_)
            else
                push!(add_others, _ex_)
            end
        end
        check = !isempty(add_delta_args)
        if check
            out = sum(add_delta_args)
            if !isempty(add_others)
                out += Sum(simplify_operators(sum(add_others)),index)
            end
        else
            out = simplify_operators(sum(add_others))
        end
        return (check, out)
    end
    return (false,ex)
end
function _resolve_kdelta(ex::TensorProd,index,arg::SymPy.Sym)
    cs = classname(arg)
    if cs=="KroneckerDelta"
        i,j = arg.__pyobject__.indices
        i_ = IndexOrder[findfirst(x->sympify(x)==i,IndexOrder)]
        j_ = IndexOrder[findfirst(x->sympify(x)==j,IndexOrder)]
        if i_.id==index.id
            out = swap_index(ex,index,j_)
            return (true, simplify_operators(out))
        elseif j_.id==index.id
            out = swap_index(ex,index,i_)
            return (true, simplify_operators(out))
        else
            return (false,ex)
        end
    elseif cs=="Mul"
        c, margs = arg.__pyobject__.as_coeff_mul()
        check = false
        for m=margs
            check, ex_ = _resolve_kdelta(ex,index,m)
            check && break
        end
        return (check, ex_)
    elseif cs=="Add"
        c, aargs = arg.__pyobject__.as_coeff_add()
        ex_args_ = [prod(ex.args[1].args[2:end]);ex.args[2:end]]
        ex_ = ex.f(ex_args_...)
        add_delta_args = []
        add_others = []
        for m=aargs
            check, _ex_ = _resolve_kdelta(m*ex_,index,m)
            if check
                push!(add_delta_args, _ex_)
            else
                push!(add_others, _ex_)
            end
        end
        check = !isempty(add_delta_args)
        if check
            out = sum(add_delta_args)
            if !isempty(add_others)
                out += Sum(simplify_operators(sum(add_others)),index)
            end
        else
            out = simplify_operators(sum(add_others))
        end
        return (check, out)
    end
    return (false,ex)
end
function _resolve_kdelta(ex::Add,index,arg::SymPy.Sym)
    checks = Bool[]
    args = []
    for a=ex.args
        check_, arg_ = _resolve_kdelta(a,index,arg)
        push!(checks, check_)
        push!(args, arg_)
    end
    return any(checks), sum(args)
end

function combine_add(s1::SumType,s2::SumType)
    i = s1.f.index
    j = s2.f.index
    if i.lower==j.lower && i.upper==j.upper
        arg1 = s1.args[1]
        arg2 = (i.label==j.label) ? s2.args[1] : swap_index(s2.args[1], j, i)
        check, op = combine_add(arg1,arg2)
        check || return (false,s1)
        return (true,Sum(op,i))
    else
        return (false,s1)
    end
end

# function _combine_add_indexed(a::Prod,b::IndexedOperator)
#     if isa(a.args[1],Number) && length(a.args)==2
#         if a.args[2]==b
#             return (a.args[1]+1)*b
#         elseif isa(a.args[2],IndexedOperator) && a.args[2].operator==b.operator && a.args[2].index.id==b.index.id
#             nid = [a.args[2].index.nid;b.index.nid]
#             unique!(nid)
#             i = Index(a.args[2].index.label,a.args[2].index.lower,a.args[2].index.upper,nid)
#             ex = (a.args[1]+1)*IndexedOperator(a.args[2].operator,i)
#
#             a_nid = filter(x->!(x∈a.args[2].index.nid), nid)
#             if !isempty(a_nid)
#                 js = IndexOrder[findall(x->x.id∈a_nid,IndexOrder)]
#                 filter!(x->isempty(x.nid),js)
#                 for j=js
#                     δ = KroneckerDelta(a.args[2].index,j)
#                     iszero(δ) && continue
#                     a_ = (a.args[1]*δ)*IndexedOperator(a.operator, j)
#                     ex += a_
#                 end
#             end
#
#             b_nid = filter(x->!(x∈b.index.nid), nid)
#             if !isempty(b_nid)
#                 ks = IndexOrder[findall(x->x.id∈b_nid,IndexOrder)]
#                 filter!(x->isempty(x.nid),ks)
#                 for k=ks
#                     δ = KroneckerDelta(b.index, k)
#                     iszero(δ) && continue
#                     b_ = δ*IndexedOperator(b.operator, k)
#                     ex += b_
#                 end
#             end
#
#             return (true,ex)
#         end
#     end
#     return (false,a)
# end
# function _combine_add_indexed(a::Prod,b::Prod)
#     if a == b
#         return true, 2*a
#     elseif isa(a.args[1],Number) && isa(b.args[1],Number)
#         if isequal_prod_args(a.args[2:end],b.args[2:end])
#             arg_ = copy(a)
#             arg_.args[1] += b.args[1]
#             return true, arg_
#         elseif length(a.args)==length(b.args)
#             a_inds = [isa(a1,IndexedOperator) for a1=a.args[2:end]]
#             b_inds =
#         end
#     elseif isa(a.args[1],Number) && isequal_prod_args(a.args[2:end], b.args)
#         arg_ = copy(a)
#         arg_.args[1] += 1
#         return true, arg_
#     elseif isa(b.args[1],Number) && isequal_prod_args(a.args, b.args[2:end])
#         arg_ = copy(b)
#         arg_.args[1] += 1
#         return true, arg_
#     else
#         return false, a
#     end
# end


function swap_index(op::IndexedOperator, i::Index, j::Index)
    @assert i.lower==j.lower && i.upper==j.upper
    if op.index==i
        _is_neq_index(i,j) && return Zero()
        # nid = [i.nid;j.nid]
        # unique!(nid)
        # j_ = Index(j.label,j.lower,j.upper,nid)
        return IndexedOperator(op.operator, j)
    elseif op.index.id==i.id # Don't have same .nid
        _is_neq_index(i,j) && return Zero()
        nid = [op.index.nid;i.nid]
        unique!(nid)
        i_ = Index(i.label,i.lower,i.upper,nid)
        return swap_index(op.operator[i_], i_,j)
    else
        return op
    end
end
swap_index(op::BasicOperator, args...) = op
function swap_index(ex::Prod, i::Index, j::Index)
    @assert i.lower==j.lower && i.upper==j.upper
    sym_inds, op_inds = find_index(ex)
    ops = filter(x->isa(x,BasicOperator),ex.args)
    cs = filter(x->isa(x,Number),ex.args)
    if id_in(i,op_inds)
        ops = [swap_index(op, i, j) for op=ops]
    end
    if id_in(i,sym_inds)
        cs = [swap_index(c, i, j) for c=cs]
    end
    args = [cs;ops]
    return Expression(ex.f,args)
end
function swap_index(ex::TensorProd, i::Index, j::Index)
    args = [swap_index(a, i, j) for a=ex.args]
    return ex.f(args...)
end
swap_index(ex::Add, i::Index, j::Index) = sum([swap_index(a, i, j) for a=ex.args])
function swap_index(ex::SumType, i::Index, j::Index)
    if ex.f.index==i
        arg = swap_index(ex.args[1], i, j)
        return Sum(arg, j)
    elseif ex.f.index.id==i.id
        nid = [ex.f.index.nid;i.nid]
        unique!(nid)
        i_ = Index(i.label,i.lower,i.upper,nid)
        arg = swap_index(ex.args[1], i, i_)
        s = Sum(arg,i_)
        @assert i_==s.f.index
        return swap_index(s, i_, j)
    else
        arg = swap_index(ex.args[1], i, j)
        return Sum(arg, ex.f.index)
    end
end

function swap_index(x::SymPy.Sym, i::Index, j::Index)
    if classname(x) == "Indexed"
        inds = x.__pyobject__.indices
        inds_ = [sympify(_find_sympy_ind_id(ind, i, j)) for ind=inds]
        b = x.__pyobject__.base
        return b[inds_...]
    elseif classname(x) == "KroneckerDelta"
        inds = x.__pyobject__.indices
        k1 = _find_sympy_ind_id(inds[1],i,j)
        k2 = _find_sympy_ind_id(inds[2],i,j)
        isa(k1,Index) && isa(k2,Index) && _is_neq_index(k1,k2) && return 0
        return KroneckerDelta(sympify(k1),sympify(k2))
    elseif classname(x) == "adjoint"
        return swap_index(x', i, j)'
    elseif classname(x) == "Mul"
        c, args = x.__pyobject__.as_coeff_mul()
        args_ = [swap_index(a, i, j) for a=args]
        return c*prod(args_)
    elseif classname(x) == "Add"
        c, args = x.__pyobject__.as_coeff_add()
        args_ = [swap_index(a, i, j) for a=args]
        return c+sum(args_)
    else
        return x
    end
end
swap_index(x::Number, args...) = x

function _find_sympy_ind_id(ind, i, j)
    is = filter(x->x.id==i.id,IndexOrder)
    is_sym = sympify.(is)
    k = ind ∈ is_sym ? j : ind
    return k
end

# function find_sympy_inds(inds, i)
#     is = IndexOrder[filter(x->x.id==i.id, IndexOrder)]
#     inds_ = SymPy.Sym[]
#     for ii=is
#         if sympify(ii)∈inds
#             push!(inds_, ii)
#         else
#             push!(inds_, )


_is_neq_index(i::Index,j::Index) = (i.id∈j.nid || j.id∈i.nid)
function _is_neq_index(i::SymPy.Sym, j::SymPy.Sym)
    i_ = IndexOrder[findall(x->sympify(x)==i, IndexOrder)]
    j_ = IndexOrder[findall(x->sympify(x)==j, IndexOrder)]
    nid = UInt[]
    for ii=i_
        append!(nid, ii.nid)
    end
    for jj=j_
        append!(nid, jj.nid)
    end
    return (i_[1].id ∈ nid || j_[1].id ∈ nid)
end

function apply_comms(s::SumType)
    arg = simplify_operators(s.args[1])
    if isa(arg, Expression)
        arg_ = simplify_operators(apply_comms(arg))
        return Sum(arg_,s.f.index)
    else
        return Sum(arg,s.f.index)
    end
end
acts_on(s::SumType) = acts_on(s.args[1])
