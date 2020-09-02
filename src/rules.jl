
let
    NC_TIMES_RULES = [
        SymbolicUtils.@rule(~x::SymbolicUtils.isnotflat(*) => SymbolicUtils.flatten_term(*, ~~x))
        SymbolicUtils.@rule(~x::needs_sorting_nc => sort_args_nc(~x))

        SymbolicUtils.ACRule(combinations, SymbolicUtils.@rule(~a::SymbolicUtils.isnumber * ~b::SymbolicUtils.isnumber => ~a * ~b), 2)

        SymbolicUtils.@rule(*(~~x::SymbolicUtils.hasrepeats) => *(SymbolicUtils.merge_repeats(^, ~~x)...))
        SymbolicUtils.ACRule(combinations, SymbolicUtils.@rule((~y)^(~n) * ~y => (~y)^(~n+1)), 3)
        SymbolicUtils.ACRule(combinations, SymbolicUtils.@rule((~x)^(~n) * (~x)^(~m) => (~x)^(~n + ~m)), 3)

        SymbolicUtils.ACRule(combinations, SymbolicUtils.@rule((~z::SymbolicUtils._isone  * ~x) => ~x), 2)
        SymbolicUtils.ACRule(combinations, SymbolicUtils.@rule((~z::SymbolicUtils._iszero *  ~x) => ~z), 2)
        SymbolicUtils.@rule(*(~x) => ~x)
    ]

    EXP_RULES = [
        # SymbolicUtils.@rule() TODO collect exp(a)*exp(b) => exp(a+b) for nc variables
        # TODO rewrite as sin/cos ?
        # SymbolicUtils.@rule(exp(~x) => cos(~x) + im*sin(~x))
    ]

    COMMUTATOR_RULES = [
        # Fock space rules
        SymbolicUtils.@rule(*(~~a, ~x::SymbolicUtils.sym_isa(Destroy), ~y::SymbolicUtils.sym_isa(Create), ~~b) => apply_commutator(commute_bosonic, ~~a, ~~b, ~x, ~y))

        # NLevel rules
        SymbolicUtils.@rule(*(~~a, ~x::SymbolicUtils.sym_isa(Transition), ~y::SymbolicUtils.sym_isa(Transition), ~~b) => apply_commutator(merge_transitions, ~~a, ~~b, ~x, ~y))
        SymbolicUtils.@rule(~x::SymbolicUtils.sym_isa(Transition) => rewrite_gs(~x))
    ]

    EXPAND_TIMES_RULES = [
        SymbolicUtils.@rule(~x::SymbolicUtils.isnotflat(*) => SymbolicUtils.flatten_term(*, ~x))
        SymbolicUtils.@rule(~x::needs_sorting_nc => sort_args_nc(~x))
        SymbolicUtils.@rule(*(~~a, +(~~b), ~~c) => +(map(b -> *((~~a)..., b, (~~c)...), ~~b)...))
    ]

    EXPAND_POW_RULES = [
        SymbolicUtils.@rule(^(~x::SymbolicUtils.sym_isa(AbstractOperator),~y::SymbolicUtils.isliteral(Integer)) => *((~x for i=1:~y)...))
    ]

    EXPAND_EXP_RULES = [

        # Euler expansion
        SymbolicUtils.@rule(cos(~x) => 0.5*exp(-im*~x) + 0.5*exp(im*~x))
        SymbolicUtils.@rule(sin(~x) => 0.5im*exp(-im*~x) + (-0.5im)*exp(im*~x))

        # Exponentials
        SymbolicUtils.@rule(exp(~x::SymbolicUtils._iszero) => 1)
        SymbolicUtils.ACRule(permutations, SymbolicUtils.@rule(exp(~x + ~y::SymbolicUtils.isnumber) => exp(~x)*exp(~y)), 2)
        SymbolicUtils.@rule(exp(~x + ~y) => separate_exp(~x,~y))
        SymbolicUtils.@rule(^(exp(~x), ~y) => exp(~y*~x))

        # Normal order in products
        SymbolicUtils.@rule(*(~~a, ~x::SymbolicUtils.sym_isa(Destroy), exp(~y::SymbolicUtils.sym_isa(AbstractOperator)), ~~b) => commute_destroy_exp(~~a, ~~b, ~x, ~y))
        SymbolicUtils.@rule(*(~~a, exp(~x::SymbolicUtils.sym_isa(AbstractOperator)), ~y::SymbolicUtils.sym_isa(Create), ~~b) => commute_exp_create(~~a, ~~b, ~x, ~y))
    ]


    # Copied directly from SymbolicUtils
    PLUS_RULES = [
        SymbolicUtils.@rule(~x::SymbolicUtils.isnotflat(+) => SymbolicUtils.flatten_term(+, ~~x))
        SymbolicUtils.@rule(~x::SymbolicUtils.needs_sorting(+) => SymbolicUtils.sort_args(+, ~~x))
        SymbolicUtils.ACRule(combinations, SymbolicUtils.@rule(~a::SymbolicUtils.isnumber + ~b::SymbolicUtils.isnumber => ~a + ~b), 2)

        SymbolicUtils.ACRule(permutations, SymbolicUtils.@rule(*(~~x) + *(~β, ~~x) => *(1 + ~β, (~~x)...)), 2)
        SymbolicUtils.ACRule(permutations, SymbolicUtils.@rule(*(~α, ~~x) + *(~β, ~~x) => *(~α + ~β, (~~x)...)), 2)

        SymbolicUtils.ACRule(permutations, SymbolicUtils.@rule(~x + *(~β, ~x) => *(1 + ~β, ~x)), 2)
        SymbolicUtils.ACRule(permutations, SymbolicUtils.@rule(*(~α::SymbolicUtils.isnumber, ~x) + ~x => *(~α + 1, ~x)), 2)
        SymbolicUtils.@rule(+(~~x::SymbolicUtils.hasrepeats) => +(SymbolicUtils.merge_repeats(*, ~~x)...))

        SymbolicUtils.ACRule(combinations, SymbolicUtils.@rule((~z::SymbolicUtils._iszero + ~x) => ~x), 2)
        SymbolicUtils.@rule(+(~x) => ~x)
    ]

    POW_RULES = [
        SymbolicUtils.@rule(^(*(~~x), ~y::SymbolicUtils.isliteral(Integer)) => *(map(a->SymbolicUtils.pow(a, ~y), ~~x)...))
        SymbolicUtils.@rule((((~x)^(~p::SymbolicUtils.isliteral(Integer)))^(~q::SymbolicUtils.isliteral(Integer))) => (~x)^((~p)*(~q)))
        SymbolicUtils.@rule(^(~x, ~z::SymbolicUtils._iszero) => 1)
        SymbolicUtils.@rule(^(~x, ~z::SymbolicUtils._isone) => ~x)
    ]

    ASSORTED_RULES = [
        SymbolicUtils.@rule(identity(~x) => ~x)
        SymbolicUtils.@rule(-(~x) => -1*~x)
        SymbolicUtils.@rule(-(~x, ~y) => ~x + -1(~y))
        SymbolicUtils.@rule(~x / ~y => ~x * SymbolicUtils.pow(~y, -1))
        SymbolicUtils.@rule(one(~x) => one(SymbolicUtils.symtype(~x)))
        SymbolicUtils.@rule(zero(~x) => zero(SymbolicUtils.symtype(~x)))
        # SymbolicUtils.@rule(SymbolicUtils.cond(~x::SymbolicUtils.isnumber, ~y, ~z) => ~x ? ~y : ~z)
    ]


    # Rewriter functions
    global operator_simplifier
    global commutator_simplifier
    global default_operator_simplifier
    global default_expand_simplifier
    global noncommutative_simplifier

    function default_operator_simplifier(; kwargs...)
        SymbolicUtils.IfElse(
            SymbolicUtils.sym_isa(AbstractOperator), SymbolicUtils.Postwalk(operator_simplifier()),
            SymbolicUtils.default_simplifier(; kwargs...)
        )
    end

    function operator_simplifier()
        rw_comms = commutator_simplifier()
        rw_nc = noncommutative_simplifier()
        rw = SymbolicUtils.Chain([rw_comms,rw_nc])
        return SymbolicUtils.Postwalk(rw)
    end

    function commutator_simplifier()
        rule_tree = [SymbolicUtils.If(SymbolicUtils.istree, SymbolicUtils.Chain(ASSORTED_RULES)),
                    SymbolicUtils.If(has_exps, SymbolicUtils.Chain(EXPAND_EXP_RULES)),
                    SymbolicUtils.If(SymbolicUtils.is_operation(*), SymbolicUtils.Chain(EXPAND_TIMES_RULES)),
                    SymbolicUtils.If(SymbolicUtils.is_operation(^), SymbolicUtils.Chain(EXPAND_POW_RULES)),
                    SymbolicUtils.Chain(COMMUTATOR_RULES)
                    ] |> SymbolicUtils.Chain
        return SymbolicUtils.Fixpoint(SymbolicUtils.Postwalk(rule_tree))
    end

    function noncommutative_simplifier()
        rule_tree = [SymbolicUtils.If(SymbolicUtils.is_operation(+), SymbolicUtils.Chain(PLUS_RULES)),
                     SymbolicUtils.If(SymbolicUtils.is_operation(*), SymbolicUtils.Chain(NC_TIMES_RULES)),
                     SymbolicUtils.If(SymbolicUtils.is_operation(^), SymbolicUtils.Chain(POW_RULES)),
                     # SymbolicUtils.If(SymbolicUtils.is_operation(exp), SymbolicUtils.Chain(EXP_RULES))
                     ] |> SymbolicUtils.Chain
        return SymbolicUtils.Fixpoint(SymbolicUtils.Postwalk(rule_tree))
    end

    function default_expand_simplifier()
        rule_tree = [SymbolicUtils.If(SymbolicUtils.istree, SymbolicUtils.Chain(ASSORTED_RULES)),
                    SymbolicUtils.If(SymbolicUtils.is_operation(*), SymbolicUtils.Chain(EXPAND_TIMES_RULES)),
                    SymbolicUtils.If(SymbolicUtils.is_operation(^), SymbolicUtils.Chain(EXPAND_POW_RULES))
                    ] |> SymbolicUtils.Chain
        return SymbolicUtils.Fixpoint(SymbolicUtils.Postwalk(rule_tree))
    end
end
