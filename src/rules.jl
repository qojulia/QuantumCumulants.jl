
# Rules
const SIMPLIFY_OPERATOR_RULES = SymbolicUtils.RuleSet([
    # SymbolicUtils.@rule ~t::SymbolicUtils.sym_isa(Number) => SymbolicUtils.NUMBER_RULES(~t, applyall=true, recurse=true)
    SymbolicUtils.@rule ~t::SymbolicUtils.sym_isa(AbstractOperator) => OPERATOR_RULES(~t, applyall=true, recurse=true)
])

const OPERATOR_RULES = SymbolicUtils.RuleSet([
    SymbolicUtils.@rule ~t               => ASSORTED_RULES(~t, recurse=false)
    SymbolicUtils.@rule ~t::SymbolicUtils.is_operation(+) =>  PLUS_RULES(~t, recurse=false)
    SymbolicUtils.@rule ~t::SymbolicUtils.is_operation(⊗) => TENSOR_RULES(~t, recurse=false)
    # SymbolicUtils.@rule ~t::SymbolicUtils.is_operation(^) =>   POW_RULES(~t, recurse=false)
    SymbolicUtils.@rule ~t::SymbolicUtils.is_operation(*) => NC_TIMES_RULES(~t, recurse=false)
    SymbolicUtils.@rule ~t               => COMMUTATOR_RULES(~t, recurse=false)
])


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


const PLUS_RULES = SymbolicUtils.RuleSet([
    SymbolicUtils.@rule(+(~~x::SymbolicUtils.isnotflat(+)) => SymbolicUtils.flatten_term(+, ~~x))
    SymbolicUtils.@rule(+(~~x::!(SymbolicUtils.issortedₑ)) => SymbolicUtils.sort_args(+, ~~x))
    SymbolicUtils.@acrule(~a::SymbolicUtils.isnumber + ~b::SymbolicUtils.isnumber => ~a + ~b)

    # *
    # SymbolicUtils.@acrule(*(~~x) + *(~β, ~~x) => *(1 + ~β, (~~x)...))
    # SymbolicUtils.@acrule(*(~α, ~~x) + *(~β, ~~x) => *(~α +  ~β, (~~x)...))
    # SymbolicUtils.@acrule(*(~~x, ~α) + *(~~x, ~β,) => *((~~x)..., ~α +  ~β))
    #
    # SymbolicUtils.@acrule(~x + *(~β, ~x) => *(1 + ~β, ~x))
    # SymbolicUtils.@acrule(*(~α::SymbolicUtils.isnumber, ~x) + ~x => *(~α + 1, ~x))

    # ⊗
    # SymbolicUtils.@rule(⊗(~~x) + ⊗(~β, ~~x) => ⊗(1 + ~β, (~~x)...))
    # SymbolicUtils.@rule(⊗(~α, ~~x) + ⊗(~β, ~~x) => ⊗(~α +  ~β, (~~x)...))
    # SymbolicUtils.@rule(⊗(~~x, ~α) + ⊗(~~x, ~β,) => ⊗((~~x)..., ~α +  ~β))

    SymbolicUtils.@rule(+(~~x::SymbolicUtils.hasrepeats) => +(SymbolicUtils.merge_repeats(*, ~~x)...))
    SymbolicUtils.@rule(+(~~x::hasrepeats) => +(merge_repeats(~~x)...))

    SymbolicUtils.@acrule((~z::SymbolicUtils._iszero + ~x) => ~x)
    SymbolicUtils.@rule(+(~x) => ~x)
])

# Handle noncommutative operators
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



const NC_TIMES_RULES = SymbolicUtils.RuleSet([
    SymbolicUtils.@rule(*(~~x::SymbolicUtils.isnotflat(*)) => SymbolicUtils.flatten_term(*, ~~x))
    SymbolicUtils.@rule(*(~~x::!(issorted_nc(*))) => sort_args_nc(*, ~~x))
    SymbolicUtils.@rule(*(~~x::has_identity) => filter_identities(*, ~~x))

    # SymbolicUtils.@rule(*(~~x::SymbolicUtils.hasrepeats) => *(SymbolicUtils.merge_repeats(^, ~~x)...))

    SymbolicUtils.@acrule(~a::SymbolicUtils.isnumber * ~b::SymbolicUtils.isnumber => ~a * ~b)
    # SymbolicUtils.@acrule((~y)^(~n) * ~y => (~y)^(~n+1))
    # SymbolicUtils.@acrule((~x)^(~n) * (~x)^(~m) => (~x)^(~n + ~m))

    SymbolicUtils.@acrule((~z::SymbolicUtils._isone  * ~x) => ~x)
    SymbolicUtils.@acrule((~z::SymbolicUtils._iszero *  ~x) => ~z)

    SymbolicUtils.@rule((~z::isidentity * ~x::isoperator) => ~x) #TODO: 2*Identity() != 2
    SymbolicUtils.@rule((~x::isoperator * ~z::isidentity) => ~x)
    SymbolicUtils.@rule(*(~x) => ~x)

    SymbolicUtils.@rule(*(~x::SymbolicUtils.sym_isa(Number), ⊗(~y, ~z)) => ⊗(*(~x,~y), ~z))
    SymbolicUtils.@rule(*(~x ⊗ ~y, ~z ⊗ ~w) => ⊗(*(~x, ~z), *(~y, ~w)))

    SymbolicUtils.@rule(*(~x, +(~y, ~z)) => +(~x*~y, ~x*~z))
    SymbolicUtils.@rule(*(+(~x, ~y), ~z) => +(~x*~z, ~y*~z))
    SymbolicUtils.@rule(*(~~x::has_inner(+)) => expand_term(*,+,~~x))
    SymbolicUtils.@rule(*(~~x::has_inner(⊗)) => expand_term(*,⊗,~~x))
])

const POW_RULES = SymbolicUtils.RuleSet([
    SymbolicUtils.@rule(^(*(~~x), ~y) => *(map(a->pow(a, ~y), ~~x)...))
    SymbolicUtils.@rule((((~x)^(~p))^(~q)) => (~x)^((~p)*(~q)))
    SymbolicUtils.@rule(^(~x, ~z::SymbolicUtils._iszero) => 1)
    SymbolicUtils.@rule(^(~x, ~z::SymbolicUtils._isone) => ~x)

    SymbolicUtils.@rule(^(~x::SymbolicUtils._iszero, ~z) => ~x)
    SymbolicUtils.@rule(^(~x::SymbolicUtils._isone, ~z) => ~x)

    # ⊗
    SymbolicUtils.@rule(^(⊗(~~x), ~y) => ⊗(map(a->(^(a, ~y)), ~~x)...))
])

const ASSORTED_RULES = SymbolicUtils.RuleSet([
    SymbolicUtils.@rule(identity(~x) => ~x)
    SymbolicUtils.@rule(-(~x, ~y) => ~x + -1(~y))
])

const TENSOR_RULES = SymbolicUtils.RuleSet([
    SymbolicUtils.@rule(⊗(~~x::SymbolicUtils.isnotflat(⊗)) => SymbolicUtils.flatten_term(⊗, ~~x))

    SymbolicUtils.@rule(⊗(~z::SymbolicUtils._iszero, ~x) => 0)
    SymbolicUtils.@rule(⊗(~x, ~z::SymbolicUtils._iszero) => 0)

    # SymbolicUtils.@acrule((~z::SymbolicUtils._iszero *  ~x) => ~z)
    SymbolicUtils.@rule(⊗(~x, *(~c::SymbolicUtils.sym_isa(Number), ~y)) => ⊗(*(~c, ~x), ~y))
    SymbolicUtils.@rule(⊗(*(~b::SymbolicUtils.sym_isa(Number), ~x),
                    *(~c::SymbolicUtils.sym_isa(Number), ~y)) => ⊗(*(*(~b, ~c) ~x), ~y))
    SymbolicUtils.@rule(⊗(~~x::has_inner(+)) => expand_term(⊗,+,~~x))
    SymbolicUtils.@rule(⊗(~x) => ~x)
])



const COMMUTATOR_RULES = SymbolicUtils.RuleSet([

])
