### Conversion to SymbolicUtils

_to_symbolic(op::BasicOperator) = OPERATORS_TO_SYMS[op](acts_on(op))
_to_symbolic(t::T) where T<:OperatorTerm = SymbolicUtils.Term{AbstractOperator}(t.f, _to_symbolic.(t.arguments))
_to_symbolic(x::Number) = x
_to_symbolic(x::SymbolicUtils.Symbolic) = x

# _to_qumulants(s::SymbolicUtils.Sym{SymbolicUtils.FnType{A,T}}) where {A,T<:AbstractOperator} = SYMS_TO_OPERATORS[s]
function _to_qumulants(t::SymbolicUtils.Term{T}) where T<:AbstractOperator
    if haskey(SYMS_TO_OPERATORS, t.f)
        return SYMS_TO_OPERATORS[t.f]
    else
        return OperatorTerm(t.f, _to_qumulants.(t.arguments))
    end
end
_to_qumulants(x::Number) = x

# Interfacing with SymbolicUtils
SymbolicUtils.promote_symtype(f, Ts::Type{<:AbstractOperator}...) = promote_type(AbstractOperator,Ts...)
SymbolicUtils.promote_symtype(f, T::Type{<:AbstractOperator}, Ts...) = promote_type(AbstractOperator,T)
# SymbolicUtils._isone(x::SymbolicUtils.Sym{T}) where T<:AbstractOperator = T <: Identity
# SymbolicUtils._iszero(x::SymbolicUtils.Sym{T}) where T<:AbstractOperator = T <: Zero

# Symbolic type promotion
for f in [+,-,*,/,^]
    @eval SymbolicUtils.promote_symtype(::$(typeof(f)),
                   T::Type{<:AbstractOperator},
                   S::Type{<:Number}) = AbstractOperator
    @eval SymbolicUtils.promote_symtype(::$(typeof(f)),
                   T::Type{<:Number},
                   S::Type{<:AbstractOperator}) = AbstractOperator
    @eval SymbolicUtils.promote_symtype(::$(typeof(f)),
                   T::Type{<:AbstractOperator},
                   S::Type{<:AbstractOperator}) = AbstractOperator#promote_type(T,S)
end

Base.one(x::SymbolicUtils.Symbolic{T}) where T<:AbstractOperator = 1#_to_symbolic(one(_to_qumulants(x)))
Base.zero(x::SymbolicUtils.Symbolic{T}) where T<:AbstractOperator = 0
Base.one(x::SymbolicUtils.Sym{SymbolicUtils.FnType{A,T}}) where {A,T<:AbstractOperator} = 1#_to_symbolic(one(_to_qumulants(x)))
Base.zero(x::SymbolicUtils.Sym{SymbolicUtils.FnType{A,T}}) where {A,T<:AbstractOperator} = 0

### End of interface

"""
    simplify_operators(op::AbstractOperator; rules=default_rules(), kwargs...)

Simplify an operator through standard algebraic rewriting, as well as using
fundamental commutation relations.

# Arguments
===========
* op: The operator expression to be simplified
* rules: The rules used. Defaults to SIMPLIFY_OPERATOR_RULES
* kwargs: Further arguments passed to `SymbolicUtils.simplify`.
"""
function simplify_operators(op::AbstractOperator; rules=default_rules(),
                commutator_rules=default_commutator_rules(), full_simplify=true,
                kwargs...)

    s = _to_symbolic(op)
    s_ = _simplify_operators(s, rules, commutator_rules; kwargs...)
    if false#full_simplify
        while !isequal(s_, s)
            s = s_
            s_ = _simplify_operators(s, rules, commutator_rules; kwargs...)
        end
    end
    (SymbolicUtils.symtype(s_) == Any) && @warn "SymbolicUtils.simplify returned symtype Any; recursion failed!"
    return _to_qumulants(s_)
end
function _simplify_operators(s, rules, commutator_rules; kwargs...)
    s_ = SymbolicUtils.simplify(s; rules=commutator_rules)
    s_ = SymbolicUtils.simplify(s_; rules=rules, kwargs...)
    return s_
end


default_rules() = SIMPLIFY_OPERATOR_RULES
default_commutator_rules() = SIMPLIFY_COMMUTATOR_RULES

function substitute(op::BasicOperator, dict)
    if haskey(dict, op)
        op_ = dict[op]
        check_hilbert(op_,op)
        return op_
    else
        return op
    end
end
substitute(t::OperatorTerm, dict) = OperatorTerm(t.f, [substitute(arg, dict) for arg in t.arguments])

### Functions needed for simplification

# # Check if an expression is equal up to number arguments
# isrepeated(x1,x2) = isequal(x1,x2)
# function isrepeated(t1::SymbolicUtils.Term,t2::SymbolicUtils.Term)
#     !(t1.f === t2.f) && return false
#     !(t1.f === (*) || t1.f === (⊗)) && (@warn "Operation other than */⊗ in isrepeated.")
#     _, args1 = separate_constants(t1)
#     _, args2 = separate_constants(t2)
#     return isequal(args1,args2)
# end
#
# function hasrepeats(x)
#     length(x) <= 1 && return false
#     for i=1:length(x)-1
#         if isrepeated(x[i], x[i+1])
#             return true
#         end
#     end
#     return false
# end
#
# separate_constants(x) = [x],[]
# separate_constants(x::SymbolicUtils.Sym{<:AbstractOperator}) = [],[x]
# function separate_constants(t::SymbolicUtils.Term)
#     c, x = (Any[],Any[])
#     for a in t.arguments
#         c_, x_ = separate_constants(a)
#         append!(c, c_)
#         append!(x, x_)
#     end
#     return c, x
# end

# # Function hierarchy * < ⊗
# for f = [:*,:⊗]
#     @eval promote_function(::typeof($f)) = ($f)
#     @eval promote_function(::typeof($f),::typeof(*)) = ($f)
#     @eval promote_function(::typeof($f),::typeof(⊗)) = (⊗)
# end
# promote_function(f1,f2,fs...) = promote_function(promote_function(f1,f2),fs...)
# promote_function(f::Function) = f
# function promote_function(xs)
#     fs = get_function.(xs)
#     return promote_function(fs...)
# end
# get_function(x) = (*)
# get_function(t::SymbolicUtils.Term) = (t.f)

# # Merge repeats even when number arguments don't match; needed for the flat expression structure
# merge_repeats(xs) = merge_repeats(promote_function(xs),xs)
# function merge_repeats(f,xs)
#     merged = Any[]
#     i=1
#     while i<=length(xs)
#         l = 1
#         c1, x1 = separate_constants(xs[i])
#         c = isempty(c1) ? 1 : prod(c1)
#         for j=i+1:length(xs)
#             c2, x2 = separate_constants(xs[j])
#             if isequal(x1,x2)
#                 c = isempty(c2) ? (c+1) : +(c,prod(c2))
#                 l += 1
#             else
#                 break
#             end
#         end
#         if l > 1
#             if !(SymbolicUtils._isone(c))
#                 push!(merged, c*f(x1...))
#             else
#                 push!(merged, f(x1...))
#             end
#         else
#             push!(merged, xs[i])
#         end
#         i+=l
#     end
#     return merged
# end


# Handle noncommutative multiplication
iscommutative(::typeof(*), x::Union{SymbolicUtils.Symbolic{SymbolicUtils.FnType{A,T}},T}) where {A,T<:AbstractOperator} = false
# iscommutative(::typeof(*), x::Union{SymbolicUtils.Symbolic{SymbolicUtils.FnType{A,T}},T}) where {A,T<:Identity} = true
# iscommutative(::typeof(*), x::Union{SymbolicUtils.Symbolic{SymbolicUtils.FnType{A,T}},T}) where {A,T<:Zero} = true
iscommutative(::typeof(*), x::SymbolicUtils.Symbolic{T}) where {T<:AbstractOperator} = false
iscommutative(::typeof(*), x::Union{SymbolicUtils.Symbolic{T},T}) where {T<:Number} = true
iscommutative(f) = x -> iscommutative(f, x)

issorted_nc(f) = x -> issorted_nc(f, x)
function issorted_nc(f::typeof(*), args)
    is_c = iscommutative.(f, args)
    args_c = args[is_c]
    args_nc = args[.!is_c]
    check1 = issorted(is_c, lt=(>))
    check2 = SymbolicUtils.issortedₑ(args_c)
    check3 = issorted(args_nc, by=acts_on)
    return  check1 && check2 && check3
end
function acts_on(t::SymbolicUtils.Term{T}) where T<:AbstractOperator
    i = findfirst(isoperator, t.arguments)
    return acts_on(t.arguments[i])
end
function acts_on(t::SymbolicUtils.Term{T}) where T<:BasicOperator
    return t.arguments[1]
end

function sort_args_nc(f::typeof(*), args)
    is_c = iscommutative.(f, args)
    args_c = SymbolicUtils.sort_args(f, args[is_c]).arguments
    args_nc = sort(args[.!is_c], by=acts_on)
    return SymbolicUtils.Term(f, [args_c; args_nc])
end

# Expand products involving sums into sums of products, e.g. x*(y+z) => x*y + x*z; needed to apply commutator rules
has_inner(f) = x -> has_inner(f,x)
function has_inner(f,args)
    for t in args
        if t isa SymbolicUtils.Term && SymbolicUtils.operation(t) === (f)
            return true
        end
    end
    return false
end

# function expand_term(f_outer, f_inner, args)
#     expanded_args = []
#     for t in args
#         if t isa SymbolicUtils.Term && SymbolicUtils.operation(t) === f_inner
#             i = findfirst(isequal(t),args)
#             lhs_args = args[1:i-1]
#             rhs_args = args[i+1:end]
#             for arg_ in SymbolicUtils.arguments(t)
#                 push!(expanded_args, SymbolicUtils.Term(f_outer, vcat(lhs_args, arg_, rhs_args)))
#             end
#             return SymbolicUtils.Term(f_inner, expanded_args)
#         end
#     end
# end
function expand_term(f_outer, f_inner, args)
    expanded_args = []
    for i=1:length(args)
        if args[i] isa SymbolicUtils.Term && SymbolicUtils.operation(args[i]) === f_inner
            inner_args = SymbolicUtils.arguments(args[i])
            for j = 1:i-1
                for inner_arg in inner_args
                    push!(expanded_args, f_outer(args[j],inner_arg))
                end
            end
            for j=i+1:length(args)
                for inner_arg in inner_args
                    push!(expanded_args, f_outer(inner_arg,args[j]))
                end
            end
            return SymbolicUtils.term(f_inner, expanded_args...)
        end
    end
end

# function expand_term(f_outer, f_inner::typeof(⊗), args)
#     expanded_args = Any[]
#     cs = filter(SymbolicUtils.sym_isa(Number), args)
#     ops = filter(!SymbolicUtils.sym_isa(Number), args)
#     for i=1:length(SymbolicUtils.arguments(ops[1]))
#         inner_args = Any[]
#         for t in ops
#             push!(inner_args, SymbolicUtils.arguments(t)[i])
#         end
#         push!(expanded_args, SymbolicUtils.Term(f_outer, inner_args))
#     end
#     if !isempty(cs)
#         c = f_outer(cs...)
#         expanded_args[1] = SymbolicUtils.Term(expanded_args[1].f, [c;SymbolicUtils.arguments(expanded_args[1])])
#     end
#     return SymbolicUtils.Term(f_inner, expanded_args)
# end

# # Remove identities from products
# function has_identity(args)
#     ops = filter(isoperator,args)
#     length(ops) <= 1 && return false
#     for arg in ops
#         isidentity(arg) && return true
#     end
#     return false
# end
# function filter_identities(f,args)
#     filtered = filter(!isidentity,args)
#     if isempty(filtered) || all(SymbolicUtils.sym_isa(Number).(filtered))
#         push!(filtered, args[findfirst(isidentity,args)])
#     end
#     return SymbolicUtils.Term(f, filtered)
# end
