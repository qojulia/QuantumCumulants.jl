const SIMPLIFY_OPERATOR_RULES = SymbolicUtils.RuleSet([
    SymbolicUtils.@rule ~t::SymbolicUtils.sym_isa(Number) => SymbolicUtils.NUMBER_RULES(~t, applyall=true, recurse=true)
    SymbolicUtils.@rule ~t::SymbolicUtils.sym_isa(AbstractOperator) => OPERATOR_RULES(~t, applyall=true, recurse=true)
])

const OPERATOR_RULES = SymbolicUtils.RuleSet([
    SymbolicUtils.@rule ~t               => SymbolicUtils.ASSORTED_RULES(~t, recurse=false)
    SymbolicUtils.@rule ~t::SymbolicUtils.is_operation(+) =>  SymbolicUtils.PLUS_RULES(~t, recurse=false)
    # SymbolicUtils.@rule ~t::SymbolicUtils.is_operation(⊗) => TENSOR_RULES(~t, recurse=false)
    # SymbolicUtils.@rule ~t::SymbolicUtils.is_operation(^) =>   POW_RULES(~t, recurse=false)
    SymbolicUtils.@rule ~t::SymbolicUtils.is_operation(*) => NC_TIMES_RULES(~t, recurse=false)
    # SymbolicUtils.@rule ~t               => COMMUTATOR_RULES(~t, recurse=false)
])

# const PLUS_RULES = SymbolicUtils.RuleSet([
#     SymbolicUtils.@rule(+(~~x::SymbolicUtils.isnotflat(+)) => SymbolicUtils.flatten_term(+, ~~x))
#     SymbolicUtils.@rule(+(~~x::!(SymbolicUtils.issortedₑ)) => SymbolicUtils.sort_args(+, ~~x))
#     SymbolicUtils.@acrule(~a::SymbolicUtils.isnumber + ~b::SymbolicUtils.isnumber => ~a + ~b)
#
#     # *
#     SymbolicUtils.@acrule(*(~~x) + *(~β, ~~x) => *(1 + ~β, (~~x)...))
#     SymbolicUtils.@acrule(*(~α, ~~x) + *(~β, ~~x) => *(~α +  ~β, (~~x)...))
#     SymbolicUtils.@acrule(*(~~x, ~α) + *(~~x, ~β,) => *((~~x)..., ~α +  ~β))
#
#     SymbolicUtils.@acrule(~x + *(~β, ~x) => *(1 + ~β, ~x))
#     # SymbolicUtils.@acrule(*(~α::SymbolicUtils.isnumber, ~x) + ~x => *(~α + 1, ~x))
#
#     # ⊗
#     # SymbolicUtils.@rule(⊗(~~x) + ⊗(~β, ~~x) => ⊗(1 + ~β, (~~x)...))
#     # SymbolicUtils.@rule(⊗(~α, ~~x) + ⊗(~β, ~~x) => ⊗(~α +  ~β, (~~x)...))
#     # SymbolicUtils.@rule(⊗(~~x, ~α) + ⊗(~~x, ~β,) => ⊗((~~x)..., ~α +  ~β))
#
#     SymbolicUtils.@rule(+(~~x::SymbolicUtils.hasrepeats) => +(SymbolicUtils.merge_repeats(*, ~~x)...))
#
#     SymbolicUtils.@acrule((~z::SymbolicUtils._iszero + ~x) => ~x)
#     SymbolicUtils.@rule(+(~x) => ~x)
# ])


const NC_TIMES_RULES = SymbolicUtils.RuleSet([
    SymbolicUtils.@rule(*(~~x::SymbolicUtils.isnotflat(*)) => SymbolicUtils.flatten_term(*, ~~x))
    SymbolicUtils.@rule(*(~~x::!(issorted_nc(*))) => sort_args_nc(*, ~~x))

    SymbolicUtils.@acrule(~a::SymbolicUtils.isnumber * ~b::SymbolicUtils.isnumber => ~a * ~b)

    # SymbolicUtils.@rule(*(~~x::SymbolicUtils.hasrepeats) => *(SymbolicUtils.merge_repeats(^, ~~x)...))

    # SymbolicUtils.@acrule((~y)^(~n) * ~y => (~y)^(~n+1))
    # SymbolicUtils.@acrule((~x)^(~n) * (~x)^(~m) => (~x)^(~n + ~m))

    SymbolicUtils.@acrule((~z::SymbolicUtils._isone  * ~x) => ~x)
    SymbolicUtils.@acrule((~z::SymbolicUtils._iszero *  ~x) => ~z)

    # SymbolicUtils.@rule((~z::isidentity * ~x::isoperator) => ~x) #TODO: 2*Identity() != 2
    # SymbolicUtils.@rule((~x::isoperator * ~z::isidentity) => ~x)
    SymbolicUtils.@rule(*(~x) => ~x)

    # SymbolicUtils.@rule(*(~x::SymbolicUtils.sym_isa(Number), ⊗(~y, ~z)) => ⊗(*(~x,~y), ~z))
    # SymbolicUtils.@rule(*(~x ⊗ ~y, ~z ⊗ ~w) => ⊗(*(~x, ~z), *(~y, ~w)))

])
#
# const POW_RULES = SymbolicUtils.RuleSet([
#     SymbolicUtils.@rule(^(*(~~x), ~y) => *(map(a->pow(a, ~y), ~~x)...))
#     SymbolicUtils.@rule((((~x)^(~p))^(~q)) => (~x)^((~p)*(~q)))
#     SymbolicUtils.@rule(^(~x, ~z::SymbolicUtils._iszero) => 1)
#     SymbolicUtils.@rule(^(~x, ~z::SymbolicUtils._isone) => ~x)
#
#     SymbolicUtils.@rule(^(~x::SymbolicUtils._iszero, ~z) => ~x)
#     SymbolicUtils.@rule(^(~x::SymbolicUtils._isone, ~z) => ~x)
#     #
#     # # ⊗
#     # SymbolicUtils.@rule(^(⊗(~~x), ~y) => ⊗(map(a->(^(a, ~y)), ~~x)...))
# ])
#
# const ASSORTED_RULES = SymbolicUtils.RuleSet([
#     SymbolicUtils.@rule(identity(~x) => ~x)
#     SymbolicUtils.@rule(-(~x, ~y) => ~x + -1(~y))
# ])

# const TENSOR_RULES = SymbolicUtils.RuleSet([
#     SymbolicUtils.@rule(⊗(~~x::SymbolicUtils.isnotflat(⊗)) => SymbolicUtils.flatten_term(⊗, ~~x))
#
#     SymbolicUtils.@rule(⊗(~z::SymbolicUtils._iszero, ~x) => 0)
#     SymbolicUtils.@rule(⊗(~x, ~z::SymbolicUtils._iszero) => 0)
#
#     # SymbolicUtils.@acrule((~z::SymbolicUtils._iszero *  ~x) => ~z)
#     SymbolicUtils.@rule(⊗(~x, *(~c::SymbolicUtils.sym_isa(Number), ~y)) => ⊗(*(~c, ~x), ~y))
#     SymbolicUtils.@rule(⊗(*(~b::SymbolicUtils.sym_isa(Number), ~x),
#                     *(~c::SymbolicUtils.sym_isa(Number), ~y)) => ⊗(*(*(~b, ~c) ~x), ~y))
#     SymbolicUtils.@rule(⊗(~~x::has_inner(+)) => expand_term(⊗,+,~~x))
#     SymbolicUtils.@rule(⊗(~x) => ~x)
# ])


# Added incrementally depending on what Hilbert spaces are defined
const SIMPLIFY_COMMUTATOR_RULES = SymbolicUtils.RuleSet([
    SymbolicUtils.@rule(~t => SymbolicUtils.ASSORTED_RULES(~t, applyall=true, recurse=true))
    SymbolicUtils.@rule(~t => EXPAND_RULES(~t, applyall=true, recurse=true))
    SymbolicUtils.@rule(~t => COMMUTATOR_RULES(~t, applyall=true, recurse=true))
])
const COMMUTATOR_RULES = SymbolicUtils.RuleSet([])

const EXPAND_RULES = SymbolicUtils.RuleSet([
    # SymbolicUtils.@acrule(~a::SymbolicUtils.isnumber + ~b::SymbolicUtils.isnumber => ~a + ~b)
    # SymbolicUtils.@acrule(~a::SymbolicUtils.isnumber * ~b::SymbolicUtils.isnumber => ~a * ~b)
    SymbolicUtils.@rule(*(~~x::SymbolicUtils.isnotflat(*)) => SymbolicUtils.flatten_term(*, ~~x))
    SymbolicUtils.@rule(*(~~x::!(issorted_nc(*))) => sort_args_nc(*, ~~x))
    SymbolicUtils.@rule(*(~~x::has_inner(+)) => expand_term(*,+,~~x))
    SymbolicUtils.@rule(+(~~x::SymbolicUtils.isnotflat(+)) => SymbolicUtils.flatten_term(+,~~x))
])
