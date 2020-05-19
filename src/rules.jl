# General simplification of noncommutative operators
const SIMPLIFY_OPERATOR_RULES = SymbolicUtils.RuleSet([
    SymbolicUtils.@rule ~t::SymbolicUtils.sym_isa(Number) => SymbolicUtils.NUMBER_RULES(~t, applyall=true, recurse=true)
    SymbolicUtils.@rule ~t::SymbolicUtils.sym_isa(AbstractOperator) => OPERATOR_RULES(~t, applyall=true, recurse=true)
])

const OPERATOR_RULES = SymbolicUtils.RuleSet([
    SymbolicUtils.@rule ~t               => SymbolicUtils.ASSORTED_RULES(~t, recurse=false)
    SymbolicUtils.@rule ~t::SymbolicUtils.is_operation(+) =>  SymbolicUtils.PLUS_RULES(~t, recurse=false)
    # SymbolicUtils.@rule ~t::SymbolicUtils.is_operation(^) =>   POW_RULES(~t, recurse=false) TODO
    SymbolicUtils.@rule ~t::SymbolicUtils.is_operation(*) => NC_TIMES_RULES(~t, recurse=false)
])

const NC_TIMES_RULES = SymbolicUtils.RuleSet([
    SymbolicUtils.@rule(*(~~x::SymbolicUtils.isnotflat(*)) => SymbolicUtils.flatten_term(*, ~~x))
    SymbolicUtils.@rule(*(~~x::!(issorted_nc(*))) => sort_args_nc(*, ~~x))

    SymbolicUtils.@acrule(~a::SymbolicUtils.isnumber * ~b::SymbolicUtils.isnumber => ~a * ~b)

    # SymbolicUtils.@rule(*(~~x::SymbolicUtils.hasrepeats) => *(SymbolicUtils.merge_repeats(^, ~~x)...))
    # SymbolicUtils.@acrule((~y)^(~n) * ~y => (~y)^(~n+1))
    # SymbolicUtils.@acrule((~x)^(~n) * (~x)^(~m) => (~x)^(~n + ~m))

    SymbolicUtils.@acrule((~z::SymbolicUtils._isone  * ~x) => ~x)
    SymbolicUtils.@acrule((~z::SymbolicUtils._iszero *  ~x) => ~z)
    SymbolicUtils.@rule(*(~x) => ~x)
])

# Commutation rules
const SIMPLIFY_COMMUTATOR_RULES = SymbolicUtils.RuleSet([
    SymbolicUtils.@rule(~t => SymbolicUtils.ASSORTED_RULES(~t, applyall=true, recurse=true))
    SymbolicUtils.@rule(~t => EXPAND_RULES(~t, applyall=true, recurse=true))
    SymbolicUtils.@rule(~t => COMMUTATOR_RULES(~t, applyall=true, recurse=true))
])

# Added incrementally depending on what Hilbert spaces are defined
const COMMUTATOR_RULES = SymbolicUtils.RuleSet([])

const EXPAND_RULES = SymbolicUtils.RuleSet([
    SymbolicUtils.@acrule(~a::SymbolicUtils.isnumber + ~b::SymbolicUtils.isnumber => ~a + ~b)
    SymbolicUtils.@acrule(~a::SymbolicUtils.isnumber * ~b::SymbolicUtils.isnumber => ~a * ~b)
    SymbolicUtils.@rule(*(~~x::SymbolicUtils.isnotflat(*)) => SymbolicUtils.flatten_term(*, ~~x))
    SymbolicUtils.@rule(*(~~x::!(issorted_nc(*))) => sort_args_nc(*, ~~x))
    SymbolicUtils.@rule(*(~~x::has_inner(+)) => expand_term(*,+,~~x))
    SymbolicUtils.@rule(+(~~x::SymbolicUtils.isnotflat(+)) => SymbolicUtils.flatten_term(+,~~x))
])
