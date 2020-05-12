### Conversion to SymbolicUtils

_to_symbolic(op::BasicOperator) = OPERATORS_TO_SYMS[op]
_to_symbolic(t::T) where T<:OperatorTerm = SymbolicUtils.Term{AbstractOperator}(t.f, _to_symbolic.(t.arguments))
_to_symbolic(x::Number) = x
_to_symbolic(x::SymbolicUtils.Symbolic) = x

_to_operator(s::SymbolicUtils.Sym{T}) where T<:AbstractOperator = SYMS_TO_OPERATORS[s]
_to_operator(t::SymbolicUtils.Term) = OperatorTerm(t.f, _to_operator.(t.arguments))
_to_operator(x::Number) = x
_to_operator(s::SymbolicUtils.Symbolic{<:Number}) = s

# Interfacing with SymbolicUtils
SymbolicUtils.promote_symtype(f, Ts::Type{<:AbstractOperator}...) = promote_type(Ts...)
SymbolicUtils.promote_symtype(f, T::Type{<:AbstractOperator}, Ts...) = T
# SymbolicUtils._isone(x::SymbolicUtils.Sym{T}) where T<:AbstractOperator = T <: Identity
SymbolicUtils._iszero(x::SymbolicUtils.Sym{T}) where T<:AbstractOperator = T <: Zero

# Symbolic type promotion
for f in [+,-,*,/,^]
    @eval SymbolicUtils.promote_symtype(::$(typeof(f)),
                   T::Type{<:AbstractOperator},
                   S::Type{<:Number}) = T
    @eval SymbolicUtils.promote_symtype(::$(typeof(f)),
                   T::Type{<:Number},
                   S::Type{<:AbstractOperator}) = S
    @eval SymbolicUtils.promote_symtype(::$(typeof(f)),
                   T::Type{<:AbstractOperator},
                   S::Type{<:AbstractOperator}) = promote_type(T,S)
end

Base.one(x::SymbolicUtils.Sym{T}) where T<:BasicOperator = _to_symbolic(one(_to_operator(x)))

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
function simplify_operators(op::AbstractOperator; rules=default_rules(), kwargs...)
    s = _to_symbolic(op)
    s_ = SymbolicUtils.simplify(s; rules=rules, kwargs...)
    (SymbolicUtils.symtype(s_) == Any) && @warn "SymbolicUtils.simplify returned symtype Any; recursion failed!"
    return _to_operator(s_)
end

default_rules() = SIMPLIFY_OPERATOR_RULES

### Functions needed for simplification

# Check if an expression is equal up to number arguments
isrepeated(x1,x2) = isequal(x1,x2)
function isrepeated(t1::SymbolicUtils.Term,t2::SymbolicUtils.Term)
    !(t1.f === t2.f) && return false
    !(t1.f === (*) || t1.f === (⊗)) && (@warn "Operation other than */⊗ in isrepeated.")
    _, args1 = separate_constants(t1)
    _, args2 = separate_constants(t2)
    return isequal(args1,args2)
end

function hasrepeats(x)
    length(x) <= 1 && return false
    for i=1:length(x)-1
        if isrepeated(x[i], x[i+1])
            return true
        end
    end
    return false
end

separate_constants(x) = [x],[]
separate_constants(x::SymbolicUtils.Sym{<:AbstractOperator}) = [],[x]
separate_constants(x::SymbolicUtils.Sym{<:Number}) = [x],[]
function separate_constants(t::SymbolicUtils.Term)
    c, x = (Any[],Any[])
    for a in t.arguments
        c_, x_ = separate_constants(a)
        append!(c, c_)
        append!(x, x_)
    end
    return c, x
end

# Function hierarchy * < ⊗
for f = [:*,:⊗]
    @eval promote_function(::typeof($f)) = ($f)
    @eval promote_function(::typeof($f),::typeof(*)) = ($f)
    @eval promote_function(::typeof($f),::typeof(⊗)) = (⊗)
end
promote_function(f1,f2,fs...) = promote_function(promote_function(f1,f2),fs...)
promote_function(f::Function) = f
function promote_function(xs)
    fs = get_function.(xs)
    return promote_function(fs...)
end
get_function(x) = (*)
get_function(t::SymbolicUtils.Term) = (t.f)

# Merge repeats even when number arguments don't match; needed for the flat expression structure
merge_repeats(xs) = merge_repeats(promote_function(xs),xs)
function merge_repeats(f,xs)
    merged = Any[]
    i=1
    while i<=length(xs)
        l = 1
        c1, x1 = separate_constants(xs[i])
        c = isempty(c1) ? 1 : prod(c1)
        for j=i+1:length(xs)
            c2, x2 = separate_constants(xs[j])
            if isequal(x1,x2)
                if !isempty(c2)
                    c = +(c,prod(c2))
                end
                l += 1
            else
                break
            end
        end
        if l > 1
            if !isone(c)
                push!(merged, c*f(x1...))
            else
                push!(merged, f(x1...))
            end
        else
            push!(merged, xs[i])
        end
        i+=l
    end
    return merged
end


# Handle noncommutative multiplication
iscommutative(::typeof(*), x::Union{SymbolicUtils.Symbolic{T},T}) where T<:AbstractOperator = false
iscommutative(::typeof(*), x::Union{SymbolicUtils.Symbolic{T},T}) where T<:Identity = true
iscommutative(::typeof(*), x::Union{SymbolicUtils.Symbolic{T},T}) where T<:Zero = true
iscommutative(::typeof(*), x::Union{SymbolicUtils.Symbolic{T},T}) where T<:Number = true
iscommutative(f) = x -> iscommutative(f, x)

issorted_nc(f) = x -> issorted_nc(f, x)
function issorted_nc(f::typeof(*), args)
    is_c = iscommutative.(f, args)
    return (issorted(is_c, lt=(>)) && SymbolicUtils.issortedₑ(args[is_c]))
end

function sort_args_nc(f::typeof(*), args)
    is_c = iscommutative.(f, args)
    args_c = SymbolicUtils.sort_args(f, args[is_c]).arguments
    return SymbolicUtils.Term(f, [args_c; args[.!is_c]])
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

function expand_term(f_outer, f_inner, args)
    expanded_args = []
    for t in args
        if t isa SymbolicUtils.Term && SymbolicUtils.operation(t) === f_inner
            i = findfirst(isequal(t),args)
            lhs_args = args[1:i-1]
            rhs_args = args[i+1:end]
            for arg_ in SymbolicUtils.arguments(t)
                push!(expanded_args, SymbolicUtils.Term(f_outer, vcat(lhs_args, arg_, rhs_args)))
            end
            return SymbolicUtils.Term(f_inner, expanded_args)
        end
    end
end
function expand_term(f_outer, f_inner::typeof(⊗), args)
    expanded_args = Any[]
    cs = filter(SymbolicUtils.sym_isa(Number), args)
    ops = filter(!SymbolicUtils.sym_isa(Number), args)
    for i=1:length(SymbolicUtils.arguments(ops[1]))
        inner_args = Any[]
        for t in ops
            push!(inner_args, SymbolicUtils.arguments(t)[i])
        end
        push!(expanded_args, SymbolicUtils.Term(f_outer, inner_args))
    end
    if !isempty(cs)
        c = f_outer(cs...)
        expanded_args[1] = SymbolicUtils.Term(expanded_args[1].f, [c;SymbolicUtils.arguments(expanded_args[1])])
    end
    return SymbolicUtils.Term(f_inner, expanded_args)
end

# Remove identities from products
function has_identity(args)
    ops = filter(isoperator,args)
    length(ops) <= 1 && return false
    for arg in ops
        isidentity(arg) && return true
    end
    return false
end
function filter_identities(f,args)
    filtered = filter(!isidentity,args)
    if isempty(filtered) || all(SymbolicUtils.sym_isa(Number).(filtered))
        push!(filtered, args[findfirst(isidentity,args)])
    end
    return SymbolicUtils.Term(f, filtered)
end
