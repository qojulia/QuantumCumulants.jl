import SymPy.Introspection: classname

mutable struct Sum{I} <: Function
    index::I
end
==(s1::Sum,s2::Sum) = (s1.index==s2.index)
const SumType{argType} = Expression{<:Sum,argType}

function Sum(op::IndexedOperator,i::Index)
    if op.index==i
        return Expression(Sum(i),[op])
    else
        u = i.upper
        N = isa(u,Symbol) ? SymPy.symbols(string(u),integer=true) : u
        return N*op
    end
end
Sum(op::IndexedOperator) = Expression(Sum(op.index),[op])
function Sum(op::BasicOperator,i::Index)
    u = i.upper
    N = isa(u,Symbol) ? SymPy.symbols(string(u),integer=true) : u
    return N*op
end
function Sum(ex::Prod,i::Index)
    sym_inds, op_inds = find_index(ex)
    if i∈op_inds || i∈sym_inds
        return Expression(Sum(i),[ex])
    else
        u = i.upper
        N = isa(u,Symbol) ? SymPy.symbols(string(u),integer=true) : u
        return N*ex
    end
end
function Sum(ex::TensorProd,i::Index)
    # ex = simplify_operators(ex_)
    sym_inds, op_inds = find_index(ex)
    if i∈op_inds || i∈sym_inds
        return Expression(Sum(i),[ex])
    else
        u = i.upper
        N = isa(u,Symbol) ? SymPy.symbols(string(u),integer=true) : u
        return N*ex
    end
end
Sum(ex::Add,i::Index) = sum([Sum(a1,i) for a1=ex.args])

function Sum(s::SumType,i::Index)
    sym_inds, op_inds = find_index(s)
    if i∈op_inds || i∈sym_inds
        inds = [s.f.index,i]
        p = sortperm(indexin(inds, IndexOrder)) # ensure same order in double sums
        inds_ = inds[p]
        s_ = Expression(Sum(inds_[2]), s.args)
        return Expression(Sum(inds_[1]),[s_])
    else
        u = i.upper
        N = isa(u,Symbol) ? SymPy.symbols(string(u),integer=true) : u
        return N*s
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

remove_NeqIndsProd(s::SumType) = Sum(remove_NeqIndsProd(s.args[1]),s.f.index)

# Algebra
# +
+(s1::SumType,s2::SumType) = Expression(+,[s1,s2])
-(s::SumType) = Sum(-s.args[1],s.f.index)

# *
*(s1::SumType,s2::SumType) = Sum(Sum(s1.args[1]*s2.args[1],s1.f.index),s2.f.index)
function *(op::IndexedOperator,s::SumType)
    (op.index.lower==s.f.index.lower && op.index.upper==s.f.index.upper) || error("Something went wrong here!")
    if op.index==s.f.index
        i = gen_index(op.index.lower,op.index.upper)
        sum_arg = swap_index(s.args[1],s.f.index,i)
        return op*Sum(sum_arg,i)
    else
        return Sum(op*s.args[1],s.f.index)
    end
end
function *(s::SumType,op::IndexedOperator)
    (op.index.lower==s.f.index.lower && op.index.upper==s.f.index.upper) || error("Something went wrong here!")
    if op.index==s.f.index
        i = gen_index(op.index.lower,op.index.upper)
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
        i = gen_index(s.f.index.lower,s.f.index.upper)
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
        i = gen_index(s.f.index.lower,s.f.index.upper)
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

# ⊗
⊗(s1::SumType,s2::SumType) = Sum(Sum(s1.args[1]⊗s2.args[1],s1.f.index),s2.f.index)
⊗(s::SumType,op::AbstractOperator) = Sum(s.args[1]⊗op,s.f.index)
⊗(op::AbstractOperator,s::SumType) = Sum(op⊗s.args[1],s.f.index)

# Simplification
function simplify_operators(s::SumType)
    s_index = s.f.index
    op = simplify_operators(s.args[1])
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
        check, arg = _resolve_kdelta(op,s_index)
        if check
            return arg
        else
            return Sum(arg,s_index)
        end
    end
end

function _resolve_kdelta(ex::SumType,index)
    check, arg = _resolve_kdelta(ex.args[1],index)
    !check && return (false, Sum(arg, ex.f.index))
    return (check, Sum(arg, ex.f.index))
end
function _resolve_kdelta(ex::Prod,index)
    if isa(ex.args[1],SymPy.Sym)
        arg = SymPy.expand(ex.args[1])
        iszero(arg) && return (true,Zero())
        return _resolve_kdelta(ex,index,arg)
    end
    return (false,ex)
end
_resolve_kdelta(op::BasicOperator,args...) = (false,op)
function _resolve_kdelta(ex::TensorProd,index)
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


function swap_index(op::IndexedOperator, i::Index, j::Index)
    @assert i.lower==j.lower && i.upper==j.upper
    if op.index==i
        return IndexedOperator(op.operator, j)
    else
        return op
    end
end
swap_index(op::BasicOperator, args...) = op
function swap_index(ex::Prod, i::Index, j::Index)
    @assert i.lower==j.lower && i.upper==j.upper
    sym_inds, op_inds = find_index(ex)
    ops = filter(x->isa(x,AbstractOperator),ex.args)
    @assert !isempty(ops)
    cs = filter(x->isa(x,Number),ex.args)
    if i∈op_inds
        ops = [swap_index(op, i, j) for op=ops]
    end
    if i∈sym_inds
        cs = [swap_index(c, i, j) for c=cs]
    end
    if !isempty(cs)
        c = prod(cs)
        return c*Expression(ex.f,ops)
    else
        return Expression(ex.f,ops)
    end
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
    else
        arg = swap_index(ex.args[1], i, j)
        return Sum(arg, ex.f.index)
    end
end
swap_index(ex::NeqIndsProd, i::Index, j::Index) = Expression(neq_inds_prod, sort_by_inds([swap_index(a,i,j) for a=ex.args]))

function swap_index(x::SymPy.Sym, i::Index, j::Index)
    if classname(x) == "Indexed"
        inds = x.__pyobject__.indices
        isym, jsym = sympify.((i,j))
        inds_ = [k==isym ? jsym : k for k=inds]
        b = x.__pyobject__.base
        return b[inds_...]
    elseif classname(x) == "KroneckerDelta"
        inds = x.__pyobject__.indices
        isym, jsym = sympify.((i,j))
        k1 = isym==inds[1] ? jsym : inds[1]
        k2 = isym==inds[2] ? jsym : inds[2]
        return KroneckerDelta(k1,k2)
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
function swap_index(x::SymPy.Sym, i::SymPy.Sym, j::SymPy.Sym)
    i_ = IndexOrder[findfirst(y->sympify(y)==i, IndexOrder)]
    j_ = IndexOrder[findfirst(y->sympify(y)==j, IndexOrder)]
    return swap_index(x, i_, j_)
end
function swap_index(x::SymPy.Sym, i::SymPy.Sym, j::Index)
    i_ = IndexOrder[findfirst(y->sympify(y)==i, IndexOrder)]
    return swap_index(x, i_, j)
end
function swap_index(x::SymPy.Sym, i::Index, j::SymPy.Sym)
    j_ = IndexOrder[findfirst(y->sympify(y)==j, IndexOrder)]
    return swap_index(x, i, j_)
end
function swap_index(x::Number, i, j)
    isa(x,SymPy.Sym) && error("swap_index(::Number, args..) has been called with a SymPy.Sym! Report this issue.")
    return x
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
