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

    SymbolicUtils.ACRule(SymbolicUtils.@rule(~a::SymbolicUtils.isnumber * ~b::SymbolicUtils.isnumber => ~a * ~b), 2)

    # SymbolicUtils.@rule(*(~~x::SymbolicUtils.hasrepeats) => *(SymbolicUtils.merge_repeats(^, ~~x)...))
    # SymbolicUtils.ACRule(SymbolicUtils.@rule((~y)^(~n) * ~y => (~y)^(~n+1)), 3)
    # SymbolicUtils.ACRule(SymbolicUtils.@rule((~x)^(~n) * (~x)^(~m) => (~x)^(~n + ~m)), 4)

    SymbolicUtils.ACRule(SymbolicUtils.@rule((~z::SymbolicUtils._isone  * ~x) => ~x), 2)
    SymbolicUtils.ACRule(SymbolicUtils.@rule((~z::SymbolicUtils._iszero *  ~x) => ~z), 2)
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
    SymbolicUtils.ACRule(SymbolicUtils.@rule(~a::SymbolicUtils.isnumber + ~b::SymbolicUtils.isnumber => ~a + ~b), 2)
    SymbolicUtils.ACRule(SymbolicUtils.@rule(~a::SymbolicUtils.isnumber * ~b::SymbolicUtils.isnumber => ~a * ~b), 2)
    SymbolicUtils.@rule(*(~~x::SymbolicUtils.isnotflat(*)) => SymbolicUtils.flatten_term(*, ~~x))
    SymbolicUtils.@rule(*(~~x::!(issorted_nc(*))) => sort_args_nc(*, ~~x))
    SymbolicUtils.@rule(*(~~x::has_inner(+)) => expand_term(*,+,~~x))
    SymbolicUtils.@rule(+(~~x::SymbolicUtils.isnotflat(+)) => SymbolicUtils.flatten_term(+,~~x))
])
