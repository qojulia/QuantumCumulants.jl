
# Rules
const SIMPLIFY_OPERATOR_RULES = SymbolicUtils.RuleSet([
    SymbolicUtils.@rule ~t::SymbolicUtils.sym_isa(Number) => SymbolicUtils.NUMBER_RULES(~t, applyall=true, recurse=true)
    SymbolicUtils.@rule ~t::SymbolicUtils.sym_isa(AbstractOperator) => OPERATOR_RULES(~t, applyall=true, recurse=true)
])

const OPERATOR_RULES = SymbolicUtils.RuleSet([
    SymbolicUtils.@rule ~t               => ASSORTED_RULES(~t, recurse=false)
    SymbolicUtils.@rule ~t::SymbolicUtils.is_operation(+) =>  PLUS_RULES(~t, recurse=false)
    SymbolicUtils.@rule ~t::SymbolicUtils.is_operation(^) =>   POW_RULES(~t, recurse=false)
    SymbolicUtils.@rule ~t::SymbolicUtils.is_operation(*) => NC_TIMES_RULES(~t, recurse=false)
])


const PLUS_RULES = SymbolicUtils.RuleSet([
    SymbolicUtils.@rule(+(~~x::SymbolicUtils.isnotflat(+)) => SymbolicUtils.flatten_term(+, ~~x))
    SymbolicUtils.@rule(+(~~x::!(SymbolicUtils.issortedₑ)) => SymbolicUtils.sort_args(+, ~~x))
    SymbolicUtils.@acrule(~a::SymbolicUtils.isnumber + ~b::SymbolicUtils.isnumber => ~a + ~b)

    # *
    SymbolicUtils.@acrule(*(~~x) + *(~β, ~~x) => *(1 + ~β, (~~x)...))
    SymbolicUtils.@acrule(*(~α, ~~x) + *(~β, ~~x) => *(~α +  ~β, (~~x)...))
    SymbolicUtils.@acrule(*(~~x, ~α) + *(~~x, ~β,) => *((~~x)..., ~α +  ~β))

    SymbolicUtils.@acrule(~x + *(~β, ~x) => *(1 + ~β, ~x))
    SymbolicUtils.@acrule(*(~α::SymbolicUtils.isnumber, ~x) + ~x => *(~α + 1, ~x))
    SymbolicUtils.@rule(+(~~x::SymbolicUtils.hasrepeats) => +(SymbolicUtils.merge_repeats(*, ~~x)...))

    # ⊗
    # SymbolicUtils.@acrule(⊗(~~x) + ⊗(~β, ~~x) => ⊗(1 + ~β, (~~x)...))
    SymbolicUtils.@rule(⊗(~α, ~~x) + ⊗(~β, ~~x) => ⊗(~α +  ~β, (~~x)...))
    SymbolicUtils.@rule(⊗(~~x, ~α) + ⊗(~~x, ~β,) => ⊗((~~x)..., ~α +  ~β))

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

const NC_TIMES_RULES = SymbolicUtils.RuleSet([
    SymbolicUtils.@rule(*(~~x::SymbolicUtils.isnotflat(*)) => SymbolicUtils.flatten_term(*, ~~x))
    SymbolicUtils.@rule(*(~~x::!(issorted_nc(*))) => sort_args_nc(*, ~~x))

    SymbolicUtils.@rule(*(~~x::SymbolicUtils.hasrepeats) => *(SymbolicUtils.merge_repeats(^, ~~x)...))

    SymbolicUtils.@acrule(~a::SymbolicUtils.isnumber * ~b::SymbolicUtils.isnumber => ~a * ~b)
    SymbolicUtils.@acrule((~y)^(~n) * ~y => (~y)^(~n+1))
    SymbolicUtils.@acrule((~x)^(~n) * (~x)^(~m) => (~x)^(~n + ~m))

    SymbolicUtils.@acrule((~z::SymbolicUtils._isone  * ~x) => ~x)
    SymbolicUtils.@acrule((~z::SymbolicUtils._iszero *  ~x) => ~z)

    # SymbolicUtils.@rule((~z::isidentity * ~x::isoperator) => ~x) TODO: 2*Identity() != 2
    SymbolicUtils.@rule(*(~x) => ~x)

    SymbolicUtils.@rule(*(~x ⊗ ~y, ~z ⊗ ~w) => ⊗(~x*~z, ~y*~w))
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

    # SymbolicUtils.@acrule((~z::SymbolicUtils._iszero *  ~x) => ~z)

    SymbolicUtils.@rule(⊗(~x) => ~x)
])

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


const SIMPLIFY_COMMUTATOR_RULES = SymbolicUtils.RuleSet([
    SymbolicUtils.@rule ~t => COMMUTATOR_RULES(~t, recurse=true, applyall=true)

    # Simplify result
    SymbolicUtils.@rule ~t               => SIMPLIFY_OPERATOR_RULES(~t, recurse=true, applyall=true)
])

const COMMUTATOR_RULES = SymbolicUtils.RuleSet([
    # Expand inner sums
    SymbolicUtils.@rule(*(~~x::has_inner(+)) => expand_term(*,+,~~x))
    # SymbolicUtils.@rule(*(~~x::has_inner(⊗)) => expand_term(*,⊗,~~x))
    SymbolicUtils.@rule(⊗(~~x::has_inner(+)) => expand_term(⊗,+,~~x))
])
