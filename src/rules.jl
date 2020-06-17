isdestroy(a) = false
iscreate(a) = false
istransition(a) = false

let
    NC_TIMES_RULES = [
        SymbolicUtils.@rule(*(~~x::SymbolicUtils.isnotflat(*)) => SymbolicUtils.flatten_term(*, ~~x))
        SymbolicUtils.@rule(*(~~x::!(issorted_nc(*))) => sort_args_nc(*, ~~x))

        SymbolicUtils.ACRule(combinations, SymbolicUtils.@rule(~a::SymbolicUtils.isnumber * ~b::SymbolicUtils.isnumber => ~a * ~b), 2)

        # SymbolicUtils.@rule(*(~~x::SymbolicUtils.hasrepeats) => *(SymbolicUtils.merge_repeats(^, ~~x)...))
        # SymbolicUtils.ACRule(permutations, SymbolicUtils.@rule((~y)^(~n) * ~y => (~y)^(~n+1)), 3)
        # SymbolicUtils.ACRule(permutations, SymbolicUtils.@rule((~x)^(~n) * (~x)^(~m) => (~x)^(~n + ~m)), 4)

        SymbolicUtils.ACRule(combinations, SymbolicUtils.@rule((~z::SymbolicUtils._isone  * ~x) => ~x), 2)
        SymbolicUtils.ACRule(combinations, SymbolicUtils.@rule((~z::SymbolicUtils._iszero *  ~x) => ~z), 2)
        SymbolicUtils.@rule(*(~x) => ~x)
    ]

    COMMUTATOR_RULES = [
        # Fock space rules
        SymbolicUtils.@rule(*(~~x::has_consecutive(isdestroy,iscreate)) => commute_bosonic(*, ~~x))
        SymbolicUtils.@rule(*(~~x::(!issorted_by_inds(isdestroy))) => sort_by_inds(*, isdestroy, ~~x))
        SymbolicUtils.@rule(*(~~x::(!issorted_by_inds(iscreate))) => sort_by_inds(*, iscreate, ~~x))

        # NLevel rules
        SymbolicUtils.@rule(*(~~x::has_consecutive(istransition)) => merge_transitions(*, ~~x))
        SymbolicUtils.@rule(~x::istransition => rewrite_gs(~x))
    ]

    EXPAND_RULES = [
        SymbolicUtils.ACRule(combinations, SymbolicUtils.@rule(~a::SymbolicUtils.isnumber + ~b::SymbolicUtils.isnumber => ~a + ~b), 2)
        SymbolicUtils.ACRule(combinations, SymbolicUtils.@rule(~a::SymbolicUtils.isnumber * ~b::SymbolicUtils.isnumber => ~a * ~b), 2)
        SymbolicUtils.@rule(*(~~x::SymbolicUtils.isnotflat(*)) => SymbolicUtils.flatten_term(*, ~~x))
        SymbolicUtils.@rule(*(~~x::!(issorted_nc(*))) => sort_args_nc(*, ~~x))
        SymbolicUtils.@rule(*(~~x::has_inner(+)) => expand_term(*,+,~~x))
        SymbolicUtils.@rule(+(~~x::SymbolicUtils.isnotflat(+)) => SymbolicUtils.flatten_term(+,~~x))
    ]


    # Copied directly from SymbolicUtils
    PLUS_RULES = [
        SymbolicUtils.@rule(+(~~x::SymbolicUtils.isnotflat(+)) => SymbolicUtils.flatten_term(+, ~~x))
        SymbolicUtils.@rule(+(~~x::!(SymbolicUtils.issortedₑ)) => SymbolicUtils.sort_args(+, ~~x))
        SymbolicUtils.ACRule(combinations, SymbolicUtils.@rule(~a::SymbolicUtils.isnumber + ~b::SymbolicUtils.isnumber => ~a + ~b), 2)

        SymbolicUtils.ACRule(permutations, SymbolicUtils.@rule(*(~~x) + *(~β, ~~x) => *(1 + ~β, (~~x)...)), 2)
        SymbolicUtils.ACRule(permutations, SymbolicUtils.@rule(*(~α, ~~x) + *(~β, ~~x) => *(~α + ~β, (~~x)...)), 2)
        SymbolicUtils.ACRule(permutations, SymbolicUtils.@rule(*(~~x, ~α) + *(~~x, ~β) => *((~~x)..., ~α + ~β)), 2)

        SymbolicUtils.ACRule(permutations, SymbolicUtils.@rule(~x + *(~β, ~x) => *(1 + ~β, ~x)), 2)
        SymbolicUtils.ACRule(permutations, SymbolicUtils.@rule(*(~α::SymbolicUtils.isnumber, ~x) + ~x => *(~α + 1, ~x)), 2)
        SymbolicUtils.@rule(+(~~x::SymbolicUtils.hasrepeats) => +(SymbolicUtils.merge_repeats(*, ~~x)...))

        SymbolicUtils.ACRule(combinations, SymbolicUtils.@rule((~z::SymbolicUtils._iszero + ~x) => ~x), 2)
        SymbolicUtils.@rule(+(~x) => ~x)
    ]

    ASSORTED_RULES = [
        SymbolicUtils.@rule(identity(~x) => ~x)
        SymbolicUtils.@rule(-(~x) => -1*~x)
        SymbolicUtils.@rule(-(~x, ~y) => ~x + -1(~y))
        SymbolicUtils.@rule(~x / ~y => ~x * pow(~y, -1))
        SymbolicUtils.@rule(one(~x) => one(SymbolicUtils.symtype(~x)))
        SymbolicUtils.@rule(zero(~x) => zero(SymbolicUtils.symtype(~x)))
        # SymbolicUtils.@rule(SymbolicUtils.cond(~x::SymbolicUtils.isnumber, ~y, ~z) => ~x ? ~y : ~z)
    ]

    # Rewriter functions
    global operator_simplifier
    global default_commutator_simplifier
    global default_operator_simplifier
    global default_expand_simplifier

    function default_operator_simplifier(; kwargs...)
        SymbolicUtils.IfElse(
            SymbolicUtils.sym_isa(AbstractOperator), SymbolicUtils.Postwalk(operator_simplifier()),
            SymbolicUtils.default_simplifier(; kwargs...)
        )
    end

    function operator_simplifier()
        rule_tree = [SymbolicUtils.If(SymbolicUtils.istree, SymbolicUtils.Chain(ASSORTED_RULES)),
                    SymbolicUtils.If(SymbolicUtils.is_operation(+), SymbolicUtils.Chain(PLUS_RULES)),
                     SymbolicUtils.If(SymbolicUtils.is_operation(*), SymbolicUtils.Chain(NC_TIMES_RULES)),
                     # SymbolicUtils.If(is_operation(^), SymbolicUtils.Chain(POW_RULES)) TODO
                     ] |> SymbolicUtils.RestartedChain
        return rule_tree
    end

    function default_commutator_simplifier()
        rule_tree = [SymbolicUtils.If(SymbolicUtils.istree, SymbolicUtils.Chain(ASSORTED_RULES)),
                    SymbolicUtils.If(SymbolicUtils.istree, SymbolicUtils.Chain(EXPAND_RULES)),
                    SymbolicUtils.Chain(COMMUTATOR_RULES)
                    ] |> SymbolicUtils.RestartedChain
        SymbolicUtils.Postwalk(rule_tree)
    end

    function default_expand_simplifier()
        rule_tree = [SymbolicUtils.If(SymbolicUtils.istree, SymbolicUtils.Chain(ASSORTED_RULES)),
                    SymbolicUtils.If(SymbolicUtils.istree, SymbolicUtils.Chain(EXPAND_RULES))
                    ] |> SymbolicUtils.RestartedChain
        SymbolicUtils.Postwalk(rule_tree)
    end
end
